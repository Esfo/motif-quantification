use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug)]
#[command(version, about = "Build a skeleton motif atlas from a FASTA proteome.")]
struct Args {
    fasta: PathBuf,
    outdir: PathBuf,

    #[arg(long = "min", default_value_t = 5)]
    min_k: usize,

    #[arg(long = "max", default_value_t = 12)]
    max_k: usize,

    #[arg(long, default_value_t = 0.35)]
    specificity: f64,

    #[arg(long, default_value_t = 2)]
    min_proteins: usize,

    #[arg(long)]
    include_trembl: bool,

    #[arg(long)]
    include_other: bool,

    #[arg(long, default_value = "*BJOXZU")]
    exclude_aa: String,

    #[arg(long)]
    no_expand_duplicates: bool,

    #[arg(long, default_value_t = 64)]
    chunk_size: usize,

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

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
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
    new_child_count: usize,
    new_coverage_count: usize,
    same_parent_count: usize,
    total_supported: usize,
}

struct MotifWriter {
    motifs: BufWriter<File>,
    postings: BufWriter<File>,
    next_motif_id: u64,
    posting_offset: u64,
}

impl MotifWriter {
    fn new(outdir: &Path) -> Result<Self> {
        let motifs = BufWriter::new(File::create(outdir.join("motifs.tsv"))?);
        let postings = BufWriter::new(File::create(outdir.join("postings.bin"))?);

        let mut writer = Self {
            motifs,
            postings,
            next_motif_id: 0,
            posting_offset: 0,
        };

        writeln!(
            writer.motifs,
            "motif_id\tmotif_text\tprotein_count\tposting_offset\tposting_bytes"
        )?;

        Ok(writer)
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
        eprintln!("k={k}: internal specificity budget={budget}");

        let buckets = collect_endpoint_buckets(k, &groups, &args, &invalid);
        let window_hits: usize = buckets.values().map(Vec::len).sum();

        eprintln!(
            "k={k}: endpoint groups={} window/protein hits={}",
            buckets.len(),
            window_hits
        );

        let mut candidates: Vec<CandidateMotif> = buckets
            .par_iter()
            .flat_map_iter(|(_, hits)| build_candidates_for_endpoint_group(k, hits, &args, &groups))
            .collect();

        candidates.sort_by(|a, b| {
            a.motif_text
                .cmp(&b.motif_text)
                .then(a.postings.len().cmp(&b.postings.len()))
                .then(a.postings.cmp(&b.postings))
        });

        let n = candidates.len();

        for candidate in &candidates {
            writer.write_candidate(candidate)?;
        }

        eprintln!("k={k}: retained skeleton motifs={n}");
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
    if args.min_proteins < 1 {
        bail!("--min-proteins must be >= 1");
    }
    if args.chunk_size == 0 {
        bail!("--chunk-size must be >= 1");
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

fn deduplicate_records(records: &[FastaRecord], no_expand_duplicates: bool) -> Vec<RepresentativeGroup> {
    let mut by_sequence: HashMap<&str, Vec<usize>> = HashMap::new();

    for (idx, record) in records.iter().enumerate() {
        by_sequence.entry(record.sequence.as_str()).or_default().push(idx);
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

fn dynamic_child_min(
    parent_support: usize,
    next_internal_fixed: usize,
    budget: usize,
    args: &Args,
) -> usize {
    if budget == 0 {
        return usize::MAX;
    }

    let step = next_internal_fixed as f64 / budget as f64;
    let required_fraction = args.specificity * step;
    let dynamic_floor = (parent_support as f64 * required_fraction).ceil() as usize;

    args.min_proteins.max(dynamic_floor)
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

    let mut bad = 0usize;

    for &b in &bytes[..k] {
        if invalid[b as usize] {
            bad += 1;
        }
    }

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
    let mut out = Vec::new();
    let mut seen_postings: HashSet<Vec<u32>> = HashSet::new();

    let fixed_positions = if k == 2 {
        vec![0, 1]
    } else {
        vec![0, k - 1]
    };

    refine_node(
        k,
        hits.to_vec(),
        fixed_positions,
        args,
        groups,
        &mut seen_postings,
        &mut out,
    );

    out
}

fn refine_node(
    k: usize,
    entries: Vec<WindowHit>,
    fixed_positions: Vec<usize>,
    args: &Args,
    groups: &[RepresentativeGroup],
    seen_postings: &mut HashSet<Vec<u32>>,
    out: &mut Vec<CandidateMotif>,
) {
    if entries.is_empty() {
        return;
    }

    let parent_postings = postings_from_entries(&entries);

    if parent_postings.len() < args.min_proteins {
        return;
    }

    let budget = internal_budget(k, args.specificity);

    let fixed_internal = fixed_positions
        .iter()
        .filter(|&&p| p > 0 && p + 1 < k)
        .count();

    if fixed_internal >= budget {
        emit_candidate(
            k,
            &entries[0],
            &fixed_positions,
            &parent_postings,
            groups,
            seen_postings,
            out,
        );
        return;
    }

    let split = choose_best_split(
        k,
        &entries,
        &fixed_positions,
        &parent_postings,
        budget,
        fixed_internal,
        args,
        groups,
    );

    match split {
        None => {
            emit_candidate(
                k,
                &entries[0],
                &fixed_positions,
                &parent_postings,
                groups,
                seen_postings,
                out,
            );
        }
        Some(children) => {
            let child_same_as_parent = children.iter().any(|child| child.postings == parent_postings);

            if !child_same_as_parent {
                emit_candidate(
                    k,
                    &entries[0],
                    &fixed_positions,
                    &parent_postings,
                    groups,
                    seen_postings,
                    out,
                );
            }

            for child in children {
                let mut next_fixed = fixed_positions.clone();

                if !next_fixed.contains(&child.position) {
                    next_fixed.push(child.position);
                    next_fixed.sort_unstable();
                }

                refine_node(
                    k,
                    child.entries,
                    next_fixed,
                    args,
                    groups,
                    seen_postings,
                    out,
                );
            }
        }
    }
}

fn choose_best_split(
    k: usize,
    entries: &[WindowHit],
    fixed_positions: &[usize],
    parent_postings: &[u32],
    budget: usize,
    fixed_internal: usize,
    args: &Args,
    groups: &[RepresentativeGroup],
) -> Option<Vec<SplitChild>> {
    let mut best_score: Option<SplitScore> = None;
    let mut best_children: Option<Vec<SplitChild>> = None;

    let child_min = dynamic_child_min(
        parent_postings.len(),
        fixed_internal + 1,
        budget,
        args,
    );

    for position in 1..(k - 1) {
        if fixed_positions.contains(&position) {
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

        let mut same_parent_count = 0usize;
        let mut new_child_count = 0usize;
        let mut new_coverage: Vec<u32> = Vec::new();

        for child in &children {
            if child.postings == parent_postings {
                same_parent_count += 1;
            } else {
                new_child_count += 1;
                new_coverage.extend_from_slice(&child.postings);
            }
        }

        new_coverage.sort_unstable();
        new_coverage.dedup();

        let score = SplitScore {
            new_child_count,
            new_coverage_count: new_coverage.len(),
            same_parent_count,
            total_supported: children.len(),
        };

        if best_score.as_ref().map_or(true, |s| score > *s) {
            best_score = Some(score);
            best_children = Some(children);
        }
    }

    best_children
}

fn postings_from_entries(entries: &[WindowHit]) -> Vec<u32> {
    let mut postings: Vec<u32> = entries.iter().map(|e| e.protein_id).collect();
    postings.sort_unstable();
    postings.dedup();
    postings
}

fn emit_candidate(
    k: usize,
    example: &WindowHit,
    fixed_positions: &[usize],
    postings: &[u32],
    groups: &[RepresentativeGroup],
    seen_postings: &mut HashSet<Vec<u32>>,
    out: &mut Vec<CandidateMotif>,
) {
    if !seen_postings.insert(postings.to_vec()) {
        return;
    }

    let motif_text = motif_text_from_entry(k, *example, fixed_positions, groups);

    out.push(CandidateMotif {
        motif_text,
        postings: postings.to_vec(),
    });
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
    writeln!(writer, "fasta\t{}", clean_tsv(&args.fasta.display().to_string()))?;
    writeln!(writer, "protein_count\t{}", protein_count)?;
    writeln!(writer, "representative_sequence_count\t{}", representative_count)?;
    writeln!(writer, "min\t{}", args.min_k)?;
    writeln!(writer, "max\t{}", args.max_k)?;
    writeln!(writer, "specificity\t{}", args.specificity)?;
    writeln!(writer, "min_proteins_floor\t{}", args.min_proteins)?;
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
