use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(version, about = "Build a raw skeleton motif atlas from a FASTA proteome.")]
struct Args {
    fasta: PathBuf,
    outdir: PathBuf,

    #[arg(long = "min", default_value_t = 5)]
    min_k: usize,

    #[arg(long = "max", default_value_t = 12)]
    max_k: usize,

    #[arg(long, default_value_t = 0.35)]
    specificity: f64,

    #[arg(long)]
    include_trembl: bool,

    #[arg(long)]
    include_other: bool,

    #[arg(long, default_value = "*BJOXZU")]
    exclude_aa: String,

    #[arg(long)]
    no_expand_duplicates: bool,

    #[arg(long, default_value_t = 16)]
    chunk_size: usize,

    #[arg(long, default_value_t = 2048)]
    branch_min_hits: usize,

    #[arg(long)]
    threads: Option<usize>,

    #[arg(long)]
    overwrite: bool,
}

#[derive(Clone, Debug)]
struct FastaRecord {
    accession: String,
    header: String,
    sequence: String,
    source: Source,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Source {
    SwissProt,
    Trembl,
    Other,
}

impl Source {
    fn as_str(self) -> &'static str {
        match self {
            Source::SwissProt => "SwissProt",
            Source::Trembl => "TrEMBL",
            Source::Other => "Other",
        }
    }

    fn rank(self) -> u8 {
        match self {
            Source::SwissProt => 0,
            Source::Trembl => 1,
            Source::Other => 2,
        }
    }
}

#[derive(Clone, Debug)]
struct RepresentativeGroup {
    representative_protein_id: u32,
    sequence: String,
    support_protein_ids: Vec<u32>,
    all_protein_ids: Vec<u32>,
}

#[derive(Clone, Copy, Debug)]
struct WindowHit {
    group_index: usize,
    start: usize,
    protein_id: u32,
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct EndpointKey {
    first: u8,
    last: u8,
}

#[derive(Clone, Debug)]
struct CandidateMotif {
    motif_text: String,
    postings: Vec<u32>,
}

#[derive(Clone, Debug)]
struct SplitChild {
    position: usize,
    aa: u8,
    entries: Vec<WindowHit>,
    postings: Vec<u32>,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct SplitScore {
    coverage_count: usize,
    new_child_count: usize,
    total_child_count: usize,
}

#[derive(Clone, Debug)]
struct AcceptedSplit {
    position: usize,
    children: Vec<SplitChild>,
}

struct MotifWriter {
    motifs: BufWriter<File>,
    postings: BufWriter<File>,
    next_motif_id: u64,
    posting_offset: u64,
}

impl MotifWriter {
    fn new(outdir: &Path) -> Result<Self> {
        let mut motifs = BufWriter::new(File::create(outdir.join("motifs.tsv"))?);

        writeln!(
            motifs,
            "motif_id\tmotif_text\tprotein_count\tposting_offset\tposting_bytes"
        )?;

        Ok(Self {
            motifs,
            postings: BufWriter::new(File::create(outdir.join("postings.bin"))?),
            next_motif_id: 0,
            posting_offset: 0,
        })
    }

    fn write_candidate(&mut self, candidate: &CandidateMotif) -> Result<()> {
        let encoded = encode_postings(&candidate.postings);
        let nbytes = encoded.len() as u64;

        self.postings.write_all(&encoded)?;

        writeln!(
            self.motifs,
            "{}\t{}\t{}\t{}\t{}",
            self.next_motif_id,
            candidate.motif_text,
            candidate.postings.len(),
            self.posting_offset,
            nbytes
        )?;

        self.next_motif_id += 1;
        self.posting_offset += nbytes;

        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        self.motifs.flush()?;
        self.postings.flush()?;
        Ok(())
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    validate_args(&args)?;

    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .context("failed to configure Rayon thread pool")?;
    }

    prepare_output_dir(&args.outdir, args.overwrite)?;

    let records = load_records(&args)?;
    eprintln!("kept {} proteins", records.len());

    let groups = deduplicate_records(&records, args.no_expand_duplicates);
    eprintln!("kept {} representative sequences", groups.len());

    write_proteins_tsv(&args.outdir, &records, &groups)?;
    write_build_info_tsv(&args.outdir, &args, records.len(), groups.len())?;

    let invalid = build_invalid_table(&args.exclude_aa);
    let mut writer = MotifWriter::new(&args.outdir)?;

    for k in args.min_k..=args.max_k {
        let budget = internal_budget(k, args.specificity);
        let split_cutoff = required_split_coverage(args.specificity);

        eprintln!(
            "k={k}: internal specificity budget={budget} split_coverage_cutoff={split_cutoff:.3}"
        );

        let buckets = collect_endpoint_buckets(k, &groups, &args, &invalid);
        let window_hits: usize = buckets.values().map(Vec::len).sum();

        eprintln!(
            "k={k}: endpoint groups={} window/protein hits={}",
            buckets.len(),
            window_hits
        );

        let mut bucket_items: Vec<(EndpointKey, Vec<WindowHit>)> = buckets.into_iter().collect();
        bucket_items.sort_by_key(|(key, _)| *key);

        let mut retained_for_k = 0usize;

        for batch in bucket_items.chunks(args.chunk_size) {
            let mut batch_candidates: Vec<CandidateMotif> = batch
                .par_iter()
                .flat_map_iter(|(_, hits)| build_candidates_for_endpoint_group(k, hits, &args, &groups))
                .collect();

            batch_candidates.sort_by(|a, b| {
                a.motif_text
                    .cmp(&b.motif_text)
                    .then(a.postings.len().cmp(&b.postings.len()))
                    .then(a.postings.cmp(&b.postings))
            });

            for candidate in &batch_candidates {
                writer.write_candidate(candidate)?;
            }

            retained_for_k += batch_candidates.len();
        }

        eprintln!("k={k}: retained raw skeleton motifs={retained_for_k}");
    }

    writer.flush()?;
    eprintln!("done: wrote {} motifs", writer.next_motif_id);

    Ok(())
}

fn validate_args(args: &Args) -> Result<()> {
    if args.min_k < 2 {
        bail!("--min must be >= 2");
    }
    if args.max_k < args.min_k {
        bail!("--max must be >= --min");
    }
    if !(args.specificity > 0.0 && args.specificity <= 1.0) {
        bail!("--specificity must be > 0 and <= 1");
    }
    if args.chunk_size == 0 {
        bail!("--chunk-size must be >= 1");
    }
    if args.branch_min_hits == 0 {
        bail!("--branch-min-hits must be >= 1");
    }
    Ok(())
}

fn prepare_output_dir(outdir: &Path, overwrite: bool) -> Result<()> {
    if outdir.exists() {
        if !overwrite {
            bail!("output directory already exists; use --overwrite to replace it");
        }

        fs::remove_dir_all(outdir)
            .with_context(|| format!("failed to remove {}", outdir.display()))?;
    }

    fs::create_dir_all(outdir)
        .with_context(|| format!("failed to create {}", outdir.display()))?;

    Ok(())
}

fn load_records(args: &Args) -> Result<Vec<FastaRecord>> {
    let mut records = parse_fasta(&args.fasta)?;

    records.retain(|r| match r.source {
        Source::SwissProt => true,
        Source::Trembl => args.include_trembl,
        Source::Other => args.include_other,
    });

    if records.is_empty() {
        bail!("no proteins kept after source filtering");
    }

    Ok(records)
}

fn parse_fasta(path: &Path) -> Result<Vec<FastaRecord>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open FASTA {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();

    for line in reader.lines() {
        let line = line?;

        if let Some(rest) = line.strip_prefix('>') {
            if let Some(h) = header.take() {
                push_record(&mut records, h, &seq);
                seq.clear();
            }

            header = Some(rest.trim().to_string());
        } else {
            seq.push_str(line.trim());
        }
    }

    if let Some(h) = header.take() {
        push_record(&mut records, h, &seq);
    }

    Ok(records)
}

fn push_record(records: &mut Vec<FastaRecord>, header: String, seq: &str) {
    let sequence = seq.trim().to_ascii_uppercase();

    if sequence.is_empty() {
        return;
    }

    let token = header.split_whitespace().next().unwrap_or(&header);
    let parts: Vec<&str> = token.split('|').collect();

    let (source, accession) = if parts.len() >= 2 {
        let source = match parts[0] {
            "sp" => Source::SwissProt,
            "tr" => Source::Trembl,
            _ => Source::Other,
        };

        (source, parts[1].to_string())
    } else {
        (Source::Other, token.to_string())
    };

    records.push(FastaRecord {
        accession,
        header,
        sequence,
        source,
    });
}

fn deduplicate_records(
    records: &[FastaRecord],
    no_expand_duplicates: bool,
) -> Vec<RepresentativeGroup> {
    let mut by_sequence: HashMap<&str, Vec<usize>> = HashMap::new();

    for (idx, record) in records.iter().enumerate() {
        by_sequence
            .entry(record.sequence.as_str())
            .or_default()
            .push(idx);
    }

    let mut groups = Vec::with_capacity(by_sequence.len());

    for (sequence, mut ids) in by_sequence {
        ids.sort_by(|&a, &b| {
            records[a]
                .source
                .rank()
                .cmp(&records[b].source.rank())
                .then(records[a].accession.cmp(&records[b].accession))
        });

        let rep_idx = ids[0];
        let all_protein_ids: Vec<u32> = ids.iter().map(|&idx| idx as u32).collect();

        let support_protein_ids = if no_expand_duplicates {
            vec![rep_idx as u32]
        } else {
            all_protein_ids.clone()
        };

        groups.push(RepresentativeGroup {
            representative_protein_id: rep_idx as u32,
            sequence: sequence.to_string(),
            support_protein_ids,
            all_protein_ids,
        });
    }

    groups.sort_by(|a, b| {
        a.representative_protein_id
            .cmp(&b.representative_protein_id)
            .then(a.sequence.len().cmp(&b.sequence.len()))
    });

    groups
}

fn build_invalid_table(exclude_aa: &str) -> [bool; 256] {
    let mut invalid = [false; 256];

    for b in exclude_aa.bytes() {
        invalid[b as usize] = true;
        invalid[b.to_ascii_lowercase() as usize] = true;
    }

    invalid
}

fn internal_budget(k: usize, specificity: f64) -> usize {
    if k <= 2 {
        return 0;
    }

    let internal = k - 2;
    let budget = ((internal as f64) * specificity).round() as usize;

    budget.min(internal)
}

fn required_split_coverage(specificity: f64) -> f64 {
    1.0 - specificity
}

fn required_coverage_count(parent_support: usize, specificity: f64) -> usize {
    ((parent_support as f64) * required_split_coverage(specificity)).ceil() as usize
}

fn required_child_support(
    parent_support: usize,
    fixed_internal: usize,
    budget: usize,
    specificity: f64,
) -> usize {
    if budget == 0 {
        return 1;
    }

    let next_fixed_internal = fixed_internal + 1;
    let progress = (next_fixed_internal as f64) / (budget as f64);
    let required_fraction = specificity * progress;

    ((parent_support as f64) * required_fraction)
        .ceil()
        .max(1.0) as usize
}

fn collect_endpoint_buckets(
    k: usize,
    groups: &[RepresentativeGroup],
    args: &Args,
    invalid: &[bool; 256],
) -> HashMap<EndpointKey, Vec<WindowHit>> {
    groups
        .par_chunks(args.chunk_size)
        .enumerate()
        .map(|(chunk_index, chunk)| {
            let mut local: HashMap<EndpointKey, Vec<WindowHit>> = HashMap::new();
            let base = chunk_index * args.chunk_size;

            for (within, group) in chunk.iter().enumerate() {
                collect_group_windows(base + within, group, k, invalid, &mut local);
            }

            local
        })
        .reduce(HashMap::new, |mut acc, local| {
            for (key, mut hits) in local {
                acc.entry(key).or_default().append(&mut hits);
            }

            acc
        })
}

fn collect_group_windows(
    group_index: usize,
    group: &RepresentativeGroup,
    k: usize,
    invalid: &[bool; 256],
    buckets: &mut HashMap<EndpointKey, Vec<WindowHit>>,
) {
    let bytes = group.sequence.as_bytes();

    if bytes.len() < k {
        return;
    }

    let mut bad = bytes[..k]
        .iter()
        .filter(|&&b| invalid[b as usize])
        .count();

    if bad == 0 {
        add_window_hit(group_index, group, 0, k, buckets);
    }

    for start in 1..=(bytes.len() - k) {
        let leaving = bytes[start - 1];
        let entering = bytes[start + k - 1];

        if invalid[leaving as usize] {
            bad -= 1;
        }

        if invalid[entering as usize] {
            bad += 1;
        }

        if bad == 0 {
            add_window_hit(group_index, group, start, k, buckets);
        }
    }
}

fn add_window_hit(
    group_index: usize,
    group: &RepresentativeGroup,
    start: usize,
    k: usize,
    buckets: &mut HashMap<EndpointKey, Vec<WindowHit>>,
) {
    let bytes = group.sequence.as_bytes();

    let key = EndpointKey {
        first: bytes[start],
        last: bytes[start + k - 1],
    };

    let bucket = buckets.entry(key).or_default();

    for &protein_id in &group.support_protein_ids {
        bucket.push(WindowHit {
            group_index,
            start,
            protein_id,
        });
    }
}

fn build_candidates_for_endpoint_group(
    k: usize,
    hits: &[WindowHit],
    args: &Args,
    groups: &[RepresentativeGroup],
) -> Vec<CandidateMotif> {
    let fixed_positions = if k == 2 {
        vec![0, 1]
    } else {
        vec![0, k - 1]
    };

    let mut candidates = refine_node(
        k,
        hits.to_vec(),
        fixed_positions,
        Vec::new(),
        args,
        groups,
    );

    dedup_candidates(&mut candidates);
    candidates
}

fn refine_node(
    k: usize,
    entries: Vec<WindowHit>,
    fixed_positions: Vec<usize>,
    blocked_positions: Vec<usize>,
    args: &Args,
    groups: &[RepresentativeGroup],
) -> Vec<CandidateMotif> {
    if entries.is_empty() {
        return Vec::new();
    }

    let parent_postings = postings_from_entries(&entries);
    let budget = internal_budget(k, args.specificity);

    let fixed_internal = fixed_positions
        .iter()
        .filter(|&&p| p > 0 && p + 1 < k)
        .count();

    if budget == 0 || fixed_internal >= budget {
        return vec![candidate_from_node(
            k,
            &entries[0],
            &fixed_positions,
            &parent_postings,
            groups,
        )];
    }

    let child_min = required_child_support(
        parent_postings.len(),
        fixed_internal,
        budget,
        args.specificity,
    );

    let splits = choose_best_splits(
        k,
        &entries,
        &fixed_positions,
        &blocked_positions,
        &parent_postings,
        child_min,
        args,
        groups,
    );

    if splits.is_empty() {
        return vec![candidate_from_node(
            k,
            &entries[0],
            &fixed_positions,
            &parent_postings,
            groups,
        )];
    }

    let tie_positions: Vec<usize> = splits.iter().map(|split| split.position).collect();
    let mut branches: Vec<(usize, Vec<WindowHit>, Vec<usize>)> = Vec::new();

    for split in splits {
        let mut child_blocked = blocked_positions.clone();

        for &position in &tie_positions {
            if position < split.position {
                child_blocked.push(position);
            }
        }

        child_blocked.sort_unstable();
        child_blocked.dedup();

        for child in split.children {
            branches.push((child.position, child.entries, child_blocked.clone()));
        }
    }

    let parallel = branches.len() > 1 && entries.len() >= args.branch_min_hits;

    if parallel {
        branches
            .into_par_iter()
            .flat_map_iter(|(position, child_entries, child_blocked)| {
                let mut next_fixed = fixed_positions.clone();

                if !next_fixed.contains(&position) {
                    next_fixed.push(position);
                    next_fixed.sort_unstable();
                }

                refine_node(k, child_entries, next_fixed, child_blocked, args, groups)
            })
            .collect()
    } else {
        let mut out = Vec::new();

        for (position, child_entries, child_blocked) in branches {
            let mut next_fixed = fixed_positions.clone();

            if !next_fixed.contains(&position) {
                next_fixed.push(position);
                next_fixed.sort_unstable();
            }

            out.extend(refine_node(
                k,
                child_entries,
                next_fixed,
                child_blocked,
                args,
                groups,
            ));
        }

        out
    }
}

fn choose_best_splits(
    k: usize,
    entries: &[WindowHit],
    fixed_positions: &[usize],
    blocked_positions: &[usize],
    parent_postings: &[u32],
    child_min: usize,
    args: &Args,
    groups: &[RepresentativeGroup],
) -> Vec<AcceptedSplit> {
    let min_covered = required_coverage_count(parent_postings.len(), args.specificity);
    let mut best_score: Option<SplitScore> = None;
    let mut best_splits: Vec<AcceptedSplit> = Vec::new();

    for position in 1..(k - 1) {
        if fixed_positions.contains(&position) || blocked_positions.contains(&position) {
            continue;
        }

        let mut by_aa: HashMap<u8, Vec<WindowHit>> = HashMap::new();

        for &entry in entries {
            let aa = aa_at_entry(entry, position, groups);
            by_aa.entry(aa).or_default().push(entry);
        }

        let mut children = Vec::new();

        for (aa, child_entries) in by_aa {
            let postings = postings_from_entries(&child_entries);

            if postings.len() >= child_min {
                children.push(SplitChild {
                    position,
                    aa,
                    entries: child_entries,
                    postings,
                });
            }
        }

        if children.is_empty() {
            continue;
        }

        children.sort_by(|a, b| a.aa.cmp(&b.aa));

        let mut covered: Vec<u32> = Vec::new();
        let mut same_parent_count = 0usize;

        for child in &children {
            covered.extend_from_slice(&child.postings);

            if child.postings == parent_postings {
                same_parent_count += 1;
            }
        }

        covered.sort_unstable();
        covered.dedup();

        let coverage_count = covered.len();

        if coverage_count < min_covered {
            continue;
        }

        let total_child_count = children.len();
        let new_child_count = total_child_count - same_parent_count;

        let score = SplitScore {
            coverage_count,
            new_child_count,
            total_child_count,
        };

        let split = AcceptedSplit { position, children };

        match &best_score {
            None => {
                best_score = Some(score);
                best_splits.push(split);
            }
            Some(current) if score > *current => {
                best_score = Some(score);
                best_splits.clear();
                best_splits.push(split);
            }
            Some(current) if score == *current => {
                best_splits.push(split);
            }
            _ => {}
        }
    }

    best_splits.sort_by_key(|split| split.position);
    best_splits
}

fn postings_from_entries(entries: &[WindowHit]) -> Vec<u32> {
    let mut postings: Vec<u32> = entries.iter().map(|e| e.protein_id).collect();
    postings.sort_unstable();
    postings.dedup();
    postings
}

fn candidate_from_node(
    k: usize,
    example: &WindowHit,
    fixed_positions: &[usize],
    postings: &[u32],
    groups: &[RepresentativeGroup],
) -> CandidateMotif {
    CandidateMotif {
        motif_text: motif_text_from_entry(k, *example, fixed_positions, groups),
        postings: postings.to_vec(),
    }
}

fn dedup_candidates(candidates: &mut Vec<CandidateMotif>) {
    candidates.sort_by(|a, b| {
        a.motif_text
            .cmp(&b.motif_text)
            .then(a.postings.cmp(&b.postings))
    });

    candidates.dedup_by(|a, b| a.motif_text == b.motif_text && a.postings == b.postings);
}

fn aa_at_entry(entry: WindowHit, position: usize, groups: &[RepresentativeGroup]) -> u8 {
    groups[entry.group_index].sequence.as_bytes()[entry.start + position]
}

fn motif_text_from_entry(
    k: usize,
    entry: WindowHit,
    fixed_positions: &[usize],
    groups: &[RepresentativeGroup],
) -> String {
    let bytes = groups[entry.group_index].sequence.as_bytes();
    let mut fixed = vec![false; k];

    for &p in fixed_positions {
        fixed[p] = true;
    }

    let mut out = String::with_capacity(k);

    for i in 0..k {
        if fixed[i] {
            out.push(bytes[entry.start + i] as char);
        } else {
            out.push('.');
        }
    }

    out
}

fn encode_postings(postings: &[u32]) -> Vec<u8> {
    let mut out = Vec::new();

    encode_varint(postings.len() as u64, &mut out);

    let mut prev = 0u32;

    for (i, &protein_id) in postings.iter().enumerate() {
        let delta = if i == 0 {
            protein_id
        } else {
            protein_id - prev
        };

        encode_varint(delta as u64, &mut out);
        prev = protein_id;
    }

    out
}

fn encode_varint(mut value: u64, out: &mut Vec<u8>) {
    while value >= 0x80 {
        out.push((value as u8) | 0x80);
        value >>= 7;
    }

    out.push(value as u8);
}

fn write_proteins_tsv(
    outdir: &Path,
    records: &[FastaRecord],
    groups: &[RepresentativeGroup],
) -> Result<()> {
    let mut representative_for = vec![0u32; records.len()];

    for group in groups {
        for &protein_id in &group.all_protein_ids {
            representative_for[protein_id as usize] = group.representative_protein_id;
        }
    }

    let mut writer = BufWriter::new(File::create(outdir.join("proteins.tsv"))?);

    writeln!(
        writer,
        "protein_id\taccession\tsource\trepresentative_protein_id\trepresentative_accession\theader"
    )?;

    for (idx, record) in records.iter().enumerate() {
        let rep_id = representative_for[idx];
        let rep_acc = &records[rep_id as usize].accession;

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            idx,
            clean_tsv(&record.accession),
            record.source.as_str(),
            rep_id,
            clean_tsv(rep_acc),
            clean_tsv(&record.header)
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn write_build_info_tsv(
    outdir: &Path,
    args: &Args,
    protein_count: usize,
    representative_count: usize,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(outdir.join("build_info.tsv"))?);

    writeln!(writer, "key\tvalue")?;
    writeln!(
        writer,
        "fasta\t{}",
        clean_tsv(&args.fasta.display().to_string())
    )?;
    writeln!(writer, "protein_count\t{}", protein_count)?;
    writeln!(
        writer,
        "representative_sequence_count\t{}",
        representative_count
    )?;
    writeln!(writer, "min\t{}", args.min_k)?;
    writeln!(writer, "max\t{}", args.max_k)?;
    writeln!(writer, "specificity\t{}", args.specificity)?;
    writeln!(
        writer,
        "split_coverage_cutoff\t{}",
        required_split_coverage(args.specificity)
    )?;
    writeln!(writer, "child_support_policy\tdynamic_specificity_threshold")?;
    writeln!(
        writer,
        "child_support_formula\tceil(parent_support * specificity * ((fixed_internal + 1) / internal_budget))"
    )?;
    writeln!(writer, "branch_min_hits\t{}", args.branch_min_hits)?;
    writeln!(writer, "chunk_size\t{}", args.chunk_size)?;
    writeln!(writer, "split_tie_policy\tspawn_all_exact_best_splits_with_canonical_order_blocking")?;
    writeln!(writer, "position_tiebreak_policy\tnone")?;
    writeln!(writer, "dedup_policy\texact_motif_text_and_postings_only")?;
    writeln!(writer, "include_trembl\t{}", args.include_trembl)?;
    writeln!(writer, "include_other\t{}", args.include_other)?;
    writeln!(writer, "exclude_aa\t{}", clean_tsv(&args.exclude_aa))?;
    writeln!(writer, "expand_duplicates\t{}", !args.no_expand_duplicates)?;

    writer.flush()?;
    Ok(())
}

fn clean_tsv(value: &str) -> String {
    value
        .replace('\t', " ")
        .replace('\n', " ")
        .replace('\r', " ")
}
