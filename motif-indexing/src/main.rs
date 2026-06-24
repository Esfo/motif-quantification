use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/*
Skeleton-first motif atlas builder
==================================

This program builds a compact skeleton-first motif atlas.

It does not store exact peptide examples, exact peptide support rows, occurrence
events, or motif/protein edge rows. The central output is:

    skeleton motif -> compressed protein posting list

Core idea
---------

For each peptide/window length k, every valid k-window is assigned to an initial
endpoint skeleton group.

Example for k=5:

    A---G -> A...G

Inside each endpoint group, the builder recursively adds internal fixed positions
when those positions create useful protein-set structure.

Examples:

    A...G
    A.P.G
    A.PQG

This keeps the motif universe compact without enumerating and storing every exact
peptide/protein relationship.

Output files
------------

proteins.tsv
    protein_id
    accession
    source
    representative_protein_id
    is_representative

motifs.tsv
    motif_id
    motif_text
    posting_offset
    posting_bytes

postings.bin
    Concatenated delta-varint encoded sorted protein_id lists.

build_info.tsv
    Build settings and high-level counts.

Not stored
----------

The current design deliberately does not store:

    exact peptide examples
    exact peptide support rows
    occurrence counts
    motif length
    protein count
    total count
    protein -> motif reverse index
    DuckDB database
    Parquet exports
    raw build TSV event logs

Reasons:

    motif length is inherent in motif_text

    protein count is recoverable by decoding the posting list

    total occurrence count is not needed for presence-based motif grouping

    exact peptide examples are debugging material, not part of the atlas

    protein -> motif can be derived later if needed

    raw event logs recreate the disk-space problem this design avoids

No hard motif-length packing limit
----------------------------------

Earlier versions packed each window into a fixed-width integer. That created an
artificial max-k limit.

This version does not pack windows into u64 or u128. Each hit stores:

    group_index
    start_position
    protein_id

Amino acids are read directly from the original representative protein sequence
during skeleton refinement.

That means there is no artificial hard limit like k <= 12 or k <= 25. Runtime and
memory still grow with larger k, but the program no longer blocks longer motif
lengths because of integer packing.
*/

const DEFAULT_EXCLUDED_AA: &str = "*BJOXZU";

#[derive(Parser, Debug)]
#[command(author, version, about = "Build a compact skeleton-first protein motif atlas.")]
struct Args {
    /// Input protein FASTA.
    fasta: PathBuf,

    /// Output folder. Created by the program.
    outdir: PathBuf,

    /// Minimum peptide/window length.
    #[arg(long, default_value_t = 5)]
    min_k: usize,

    /// Maximum peptide/window length.
    #[arg(long, default_value_t = 12)]
    max_k: usize,

    /// Include tr| TrEMBL entries. Default is SwissProt-only.
    #[arg(long)]
    include_trembl: bool,

    /// Include FASTA entries that are neither sp| nor tr|.
    #[arg(long)]
    include_other: bool,

    /// Amino acid symbols that invalidate a window.
    #[arg(long, default_value = DEFAULT_EXCLUDED_AA)]
    exclude_aa: String,

    /// Minimum number of proteins required for a motif to be kept.
    #[arg(long, default_value_t = 2)]
    min_proteins: usize,

    /// Maximum fixed positions in a skeleton, including the two endpoints.
    #[arg(long, default_value_t = 4)]
    max_fixed: usize,

    /// Only index representative sequences, not duplicate accessions.
    #[arg(long)]
    no_expand_duplicates: bool,

    /// Representative protein groups per Rayon task during endpoint collection.
    #[arg(long, default_value_t = 64)]
    chunk_size: usize,

    /// Worker threads. Defaults to Rayon default.
    #[arg(long)]
    threads: Option<usize>,

    /// Overwrite output folder if it exists.
    #[arg(long)]
    overwrite: bool,
}

#[derive(Clone, Debug)]
struct FastaRecord {
    protein_id: u32,
    accession: String,
    source: Source,
    sequence: String,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Source {
    SwissProt,
    Trembl,
    Other,
}

#[derive(Clone, Debug)]
struct RepresentativeGroup {
    representative_id: u32,
    sequence: String,
    protein_ids: Vec<u32>,
}

#[derive(Clone, Copy, Debug)]
struct WindowHit {
    group_index: usize,
    start: usize,
    protein_id: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
struct EndpointKey {
    first: u8,
    last: u8,
}

#[derive(Clone, Debug)]
struct CandidateMotif {
    motif_text: String,
    postings: Vec<u32>,
}

struct MotifWriter {
    motifs: BufWriter<File>,
    postings: BufWriter<File>,
    next_motif_id: u64,
    posting_offset: u64,
}

impl Source {
    fn as_str(self) -> &'static str {
        match self {
            Source::SwissProt => "swissprot",
            Source::Trembl => "trembl",
            Source::Other => "other",
        }
    }
}

impl MotifWriter {
    fn new(outdir: &Path) -> Result<Self> {
        let motifs_path = outdir.join("motifs.tsv");
        let postings_path = outdir.join("postings.bin");

        let mut motifs = BufWriter::new(File::create(&motifs_path)?);
        let postings = BufWriter::new(File::create(&postings_path)?);

        writeln!(motifs, "motif_id\tmotif_text\tposting_offset\tposting_bytes")?;

        Ok(Self {
            motifs,
            postings,
            next_motif_id: 0,
            posting_offset: 0,
        })
    }

    fn write_candidate(&mut self, candidate: CandidateMotif) -> Result<()> {
        let motif_id = self.next_motif_id;
        self.next_motif_id += 1;

        let encoded_postings = encode_postings(&candidate.postings);
        let posting_bytes = encoded_postings.len() as u64;
        let posting_offset = self.posting_offset;
        self.posting_offset += posting_bytes;

        self.postings.write_all(&encoded_postings)?;

        writeln!(
            self.motifs,
            "{}\t{}\t{}\t{}",
            motif_id,
            clean_tsv(&candidate.motif_text),
            posting_offset,
            posting_bytes
        )?;

        Ok(())
    }

    fn finish(mut self) -> Result<u64> {
        self.motifs.flush()?;
        self.postings.flush()?;
        Ok(self.next_motif_id)
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    validate_args(&args)?;

    if let Some(thread_count) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build_global()
            .context("Could not initialize Rayon thread pool")?;
    }

    prepare_output_dir(&args)?;

    eprintln!("Loading FASTA...");
    let records = load_records(&args)?;
    let representative_groups = deduplicate_records(&records);

    eprintln!(
        "Kept {} proteins: {} SwissProt, {} TrEMBL, {} other.",
        records.len(),
        records
            .iter()
            .filter(|record| record.source == Source::SwissProt)
            .count(),
        records
            .iter()
            .filter(|record| record.source == Source::Trembl)
            .count(),
        records
            .iter()
            .filter(|record| record.source == Source::Other)
            .count()
    );

    eprintln!(
        "Representative sequences after exact deduplication: {}",
        representative_groups.len()
    );

    write_proteins_tsv(&args.outdir, &records, &representative_groups)?;
    write_build_info_tsv(&args, &records, &representative_groups)?;

    let invalid = build_invalid_table(&args.exclude_aa);
    let mut writer = MotifWriter::new(&args.outdir)?;

    for k in args.min_k..=args.max_k {
        eprintln!("[k={}] collecting endpoint skeleton groups...", k);
        let buckets = collect_endpoint_buckets(&args, &representative_groups, &invalid, k)?;

        let total_hits: usize = buckets.values().map(Vec::len).sum();
        eprintln!(
            "[k={}] endpoint groups: {}, window/protein hits: {}",
            k,
            buckets.len(),
            total_hits
        );

        eprintln!("[k={}] refining skeletons...", k);
        let mut candidate_batches: Vec<Vec<CandidateMotif>> = buckets
            .par_iter()
            .map(|(_endpoint, hits)| {
                build_candidates_for_endpoint_group(k, hits, &args, &representative_groups)
            })
            .collect();

        let candidate_count: usize = candidate_batches.iter().map(Vec::len).sum();
        eprintln!("[k={}] retained skeleton motifs: {}", k, candidate_count);

        for batch in candidate_batches.iter_mut() {
            for candidate in batch.drain(..) {
                writer.write_candidate(candidate)?;
            }
        }
    }

    let motif_count = writer.finish()?;
    eprintln!("Done. Wrote {} motifs.", motif_count);
    eprintln!("Output folder: {}", args.outdir.display());

    Ok(())
}

fn validate_args(args: &Args) -> Result<()> {
    if args.min_k < 2 {
        bail!("--min-k must be >= 2 because endpoint skeletons require two endpoints");
    }

    if args.max_k < args.min_k {
        bail!("--max-k must be >= --min-k");
    }

    if args.min_proteins < 2 {
        bail!("--min-proteins must be >= 2");
    }

    if args.max_fixed < 2 {
        bail!("--max-fixed must be >= 2 because endpoints are always fixed");
    }

    if args.chunk_size == 0 {
        bail!("--chunk-size must be >= 1");
    }

    Ok(())
}

fn prepare_output_dir(args: &Args) -> Result<()> {
    if args.outdir.exists() {
        if args.overwrite {
            fs::remove_dir_all(&args.outdir)
                .with_context(|| format!("Could not remove {}", args.outdir.display()))?;
        } else {
            bail!(
                "Output folder already exists: {}. Use --overwrite to replace it.",
                args.outdir.display()
            );
        }
    }

    fs::create_dir_all(&args.outdir)?;
    Ok(())
}

fn source_from_accession(accession: &str) -> Source {
    if accession.starts_with("sp|") {
        Source::SwissProt
    } else if accession.starts_with("tr|") {
        Source::Trembl
    } else {
        Source::Other
    }
}

fn load_records(args: &Args) -> Result<Vec<FastaRecord>> {
    let raw = parse_fasta(&args.fasta)?;

    if raw.is_empty() {
        bail!("No FASTA records found in {}", args.fasta.display());
    }

    let filtered: Vec<FastaRecord> = raw
        .into_iter()
        .filter(|record| match record.source {
            Source::SwissProt => true,
            Source::Trembl => args.include_trembl,
            Source::Other => args.include_other,
        })
        .enumerate()
        .map(|(protein_id, mut record)| {
            record.protein_id = protein_id as u32;
            record
        })
        .collect();

    if filtered.is_empty() {
        bail!(
            "No records remained after filtering. Default is SwissProt-only. Use --include-trembl or --include-other if needed."
        );
    }

    Ok(filtered)
}

fn parse_fasta(path: &Path) -> Result<Vec<FastaRecord>> {
    let file = File::open(path).with_context(|| format!("Could not open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut accession: Option<String> = None;
    let mut sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        if trimmed.starts_with('>') {
            if let Some(previous_accession) = accession.take() {
                let source = source_from_accession(&previous_accession);

                records.push(FastaRecord {
                    protein_id: records.len() as u32,
                    accession: previous_accession,
                    source,
                    sequence: sequence.to_ascii_uppercase(),
                });

                sequence.clear();
            }

            let header = trimmed.trim_start_matches('>');
            let first_token = header
                .split_whitespace()
                .next()
                .context("FASTA header had no accession token")?;

            accession = Some(first_token.to_string());
        } else {
            sequence.push_str(trimmed);
        }
    }

    if let Some(previous_accession) = accession.take() {
        let source = source_from_accession(&previous_accession);

        records.push(FastaRecord {
            protein_id: records.len() as u32,
            accession: previous_accession,
            source,
            sequence: sequence.to_ascii_uppercase(),
        });
    }

    Ok(records)
}

fn deduplicate_records(records: &[FastaRecord]) -> Vec<RepresentativeGroup> {
    let mut sequence_to_protein_ids: HashMap<&str, Vec<u32>> = HashMap::new();

    for record in records {
        sequence_to_protein_ids
            .entry(record.sequence.as_str())
            .or_default()
            .push(record.protein_id);
    }

    let mut groups = Vec::with_capacity(sequence_to_protein_ids.len());

    for (sequence, mut protein_ids) in sequence_to_protein_ids {
        protein_ids.sort_by_key(|protein_id| representative_rank(records, *protein_id));

        let representative_id = protein_ids[0];

        groups.push(RepresentativeGroup {
            representative_id,
            sequence: sequence.to_string(),
            protein_ids,
        });
    }

    groups.sort_by_key(|group| group.representative_id);
    groups
}

fn representative_rank(records: &[FastaRecord], protein_id: u32) -> (u8, u32) {
    let source_rank = match records[protein_id as usize].source {
        Source::SwissProt => 0,
        Source::Other => 1,
        Source::Trembl => 2,
    };

    (source_rank, protein_id)
}

fn build_invalid_table(exclude_aa: &str) -> [bool; 256] {
    let mut invalid = [true; 256];

    for byte in b'A'..=b'Z' {
        invalid[byte as usize] = false;
    }

    for byte in exclude_aa.bytes() {
        invalid[byte as usize] = true;
        invalid[byte.to_ascii_lowercase() as usize] = true;
    }

    invalid
}

fn collect_endpoint_buckets(
    args: &Args,
    groups: &[RepresentativeGroup],
    invalid: &[bool; 256],
    k: usize,
) -> Result<HashMap<EndpointKey, Vec<WindowHit>>> {
    let expand_duplicates = !args.no_expand_duplicates;

    let partial_maps: Vec<HashMap<EndpointKey, Vec<WindowHit>>> = groups
        .par_chunks(args.chunk_size)
        .enumerate()
        .map(|(chunk_index, chunk)| {
            let mut local: HashMap<EndpointKey, Vec<WindowHit>> = HashMap::new();
            let group_offset = chunk_index * args.chunk_size;

            for (within_chunk_index, group) in chunk.iter().enumerate() {
                let group_index = group_offset + within_chunk_index;

                collect_group_windows(
                    k,
                    group_index,
                    group,
                    invalid,
                    expand_duplicates,
                    &mut local,
                );
            }

            local
        })
        .collect();

    let mut merged: HashMap<EndpointKey, Vec<WindowHit>> = HashMap::new();

    for mut partial in partial_maps {
        for (key, mut hits) in partial.drain() {
            merged.entry(key).or_default().append(&mut hits);
        }
    }

    Ok(merged)
}

fn collect_group_windows(
    k: usize,
    group_index: usize,
    group: &RepresentativeGroup,
    invalid: &[bool; 256],
    expand_duplicates: bool,
    buckets: &mut HashMap<EndpointKey, Vec<WindowHit>>,
) {
    let bytes = group.sequence.as_bytes();

    if bytes.len() < k {
        return;
    }

    let output_protein_ids: &[u32] = if expand_duplicates {
        &group.protein_ids
    } else {
        std::slice::from_ref(&group.representative_id)
    };

    let mut bad = 0usize;

    for byte in &bytes[0..k] {
        if invalid[*byte as usize] {
            bad += 1;
        }
    }

    if bad == 0 {
        add_window_hit(bytes, 0, k, group_index, output_protein_ids, buckets);
    }

    for start in 1..=(bytes.len() - k) {
        let left = bytes[start - 1];
        let right = bytes[start + k - 1];

        if invalid[left as usize] {
            bad -= 1;
        }

        if invalid[right as usize] {
            bad += 1;
        }

        if bad == 0 {
            add_window_hit(bytes, start, k, group_index, output_protein_ids, buckets);
        }
    }
}

fn add_window_hit(
    bytes: &[u8],
    start: usize,
    k: usize,
    group_index: usize,
    protein_ids: &[u32],
    buckets: &mut HashMap<EndpointKey, Vec<WindowHit>>,
) {
    let endpoint = EndpointKey {
        first: bytes[start],
        last: bytes[start + k - 1],
    };

    let bucket = buckets.entry(endpoint).or_default();

    for protein_id in protein_ids {
        bucket.push(WindowHit {
            group_index,
            start,
            protein_id: *protein_id,
        });
    }
}

fn build_candidates_for_endpoint_group(
    k: usize,
    hits: &[WindowHit],
    args: &Args,
    groups: &[RepresentativeGroup],
) -> Vec<CandidateMotif> {
    if hits.is_empty() {
        return Vec::new();
    }

    let mut seen_postings: HashSet<Vec<u32>> = HashSet::new();
    let fixed_positions = vec![0, k - 1];
    let mut candidates = Vec::new();

    refine_node(
        k,
        hits.to_vec(),
        fixed_positions,
        args,
        groups,
        &mut seen_postings,
        &mut candidates,
    );

    candidates
}

fn refine_node(
    k: usize,
    entries: Vec<WindowHit>,
    fixed_positions: Vec<usize>,
    args: &Args,
    groups: &[RepresentativeGroup],
    seen_postings: &mut HashSet<Vec<u32>>,
    candidates: &mut Vec<CandidateMotif>,
) {
    let parent_postings = postings_from_entries(&entries);

    if parent_postings.len() < args.min_proteins {
        return;
    }

    if fixed_positions.len() >= args.max_fixed {
        emit_candidate(
            k,
            &entries,
            &fixed_positions,
            parent_postings,
            groups,
            seen_postings,
            candidates,
        );
        return;
    }

    let best_split = choose_best_split(
        k,
        &entries,
        &fixed_positions,
        &parent_postings,
        args,
        groups,
    );

    let Some(children) = best_split else {
        emit_candidate(
            k,
            &entries,
            &fixed_positions,
            parent_postings,
            groups,
            seen_postings,
            candidates,
        );
        return;
    };

    let has_redundant_child = children
        .iter()
        .any(|child| child.postings == parent_postings);

    if !has_redundant_child {
        emit_candidate(
            k,
            &entries,
            &fixed_positions,
            parent_postings,
            groups,
            seen_postings,
            candidates,
        );
    }

    for child in children {
        let mut child_fixed_positions = fixed_positions.clone();
        child_fixed_positions.push(child.position);
        child_fixed_positions.sort_unstable();

        refine_node(
            k,
            child.entries,
            child_fixed_positions,
            args,
            groups,
            seen_postings,
            candidates,
        );
    }
}

struct SplitChild {
    position: usize,
    entries: Vec<WindowHit>,
    postings: Vec<u32>,
}

fn choose_best_split(
    k: usize,
    entries: &[WindowHit],
    fixed_positions: &[usize],
    parent_postings: &[u32],
    args: &Args,
    groups: &[RepresentativeGroup],
) -> Option<Vec<SplitChild>> {
    let fixed: HashSet<usize> = fixed_positions.iter().copied().collect();
    let mut best_score: Option<(usize, usize, usize)> = None;
    let mut best_children: Option<Vec<SplitChild>> = None;

    for position in 1..(k - 1) {
        if fixed.contains(&position) {
            continue;
        }

        let mut partitions: Vec<Vec<WindowHit>> = vec![Vec::new(); 26];

        for entry in entries {
            let aa = aa_at_entry(groups, entry, position);

            if aa < b'A' || aa > b'Z' {
                continue;
            }

            let code = (aa - b'A') as usize;

            if code < partitions.len() {
                partitions[code].push(*entry);
            }
        }

        let mut children = Vec::new();
        let mut same_as_parent = 0usize;
        let mut nonredundant = 0usize;
        let mut supported = 0usize;

        for partition_entries in partitions {
            if partition_entries.is_empty() {
                continue;
            }

            let postings = postings_from_entries(&partition_entries);

            if postings.len() < args.min_proteins {
                continue;
            }

            supported += 1;

            if postings == parent_postings {
                same_as_parent += 1;
            } else {
                nonredundant += 1;
            }

            children.push(SplitChild {
                position,
                entries: partition_entries,
                postings,
            });
        }

        if children.is_empty() {
            continue;
        }

        if same_as_parent == 0 && nonredundant == 0 {
            continue;
        }

        let score = (same_as_parent, nonredundant, supported);

        if best_score.map_or(true, |current| score > current) {
            best_score = Some(score);
            best_children = Some(children);
        }
    }

    best_children
}

fn postings_from_entries(entries: &[WindowHit]) -> Vec<u32> {
    let mut postings: Vec<u32> = entries.iter().map(|entry| entry.protein_id).collect();
    postings.sort_unstable();
    postings.dedup();
    postings
}

fn emit_candidate(
    k: usize,
    entries: &[WindowHit],
    fixed_positions: &[usize],
    postings: Vec<u32>,
    groups: &[RepresentativeGroup],
    seen_postings: &mut HashSet<Vec<u32>>,
    candidates: &mut Vec<CandidateMotif>,
) {
    if !seen_postings.insert(postings.clone()) {
        return;
    }

    let motif_text = motif_text_from_entry(groups, entries[0], k, fixed_positions);

    candidates.push(CandidateMotif {
        motif_text,
        postings,
    });
}

fn aa_at_entry(groups: &[RepresentativeGroup], entry: &WindowHit, position: usize) -> u8 {
    groups[entry.group_index].sequence.as_bytes()[entry.start + position]
}

fn motif_text_from_entry(
    groups: &[RepresentativeGroup],
    entry: WindowHit,
    k: usize,
    fixed_positions: &[usize],
) -> String {
    let fixed: HashSet<usize> = fixed_positions.iter().copied().collect();
    let mut text = String::with_capacity(k);

    for position in 0..k {
        if fixed.contains(&position) {
            text.push(aa_at_entry(groups, &entry, position) as char);
        } else {
            text.push('.');
        }
    }

    text
}

fn encode_postings(postings: &[u32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(postings.len() * 2);
    let mut previous = 0u32;

    for (index, protein_id) in postings.iter().copied().enumerate() {
        let delta = if index == 0 {
            protein_id
        } else {
            protein_id - previous
        };

        encode_varint(delta, &mut out);
        previous = protein_id;
    }

    out
}

fn encode_varint(mut value: u32, out: &mut Vec<u8>) {
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
    let path = outdir.join("proteins.tsv");
    let mut representative_by_protein_id = HashMap::new();

    for group in groups {
        for protein_id in &group.protein_ids {
            representative_by_protein_id.insert(*protein_id, group.representative_id);
        }
    }

    let mut writer = BufWriter::new(File::create(&path)?);

    writeln!(
        writer,
        "protein_id\taccession\tsource\trepresentative_protein_id\tis_representative"
    )?;

    for record in records {
        let representative_id = representative_by_protein_id[&record.protein_id];
        let is_representative = if representative_id == record.protein_id {
            1
        } else {
            0
        };

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            record.protein_id,
            clean_tsv(&record.accession),
            record.source.as_str(),
            representative_id,
            is_representative
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn write_build_info_tsv(
    args: &Args,
    records: &[FastaRecord],
    groups: &[RepresentativeGroup],
) -> Result<()> {
    let path = args.outdir.join("build_info.tsv");
    let mut writer = BufWriter::new(File::create(&path)?);

    writeln!(writer, "key\tvalue")?;
    writeln!(writer, "input_fasta\t{}", args.fasta.display())?;
    writeln!(writer, "output_folder\t{}", args.outdir.display())?;
    writeln!(writer, "min_k\t{}", args.min_k)?;
    writeln!(writer, "max_k\t{}", args.max_k)?;
    writeln!(writer, "include_trembl\t{}", args.include_trembl)?;
    writeln!(writer, "include_other\t{}", args.include_other)?;
    writeln!(writer, "exclude_aa\t{}", args.exclude_aa)?;
    writeln!(writer, "min_proteins\t{}", args.min_proteins)?;
    writeln!(writer, "max_fixed\t{}", args.max_fixed)?;
    writeln!(
        writer,
        "expand_duplicates\t{}",
        !args.no_expand_duplicates
    )?;
    writeln!(writer, "kept_proteins\t{}", records.len())?;
    writeln!(writer, "representative_sequences\t{}", groups.len())?;

    writer.flush()?;
    Ok(())
}

fn clean_tsv(text: &str) -> String {
    text.replace('\t', " ").replace('\n', " ")
}

#[allow(dead_code)]
fn decode_postings(bytes: &[u8]) -> Vec<u32> {
    let mut postings = Vec::new();
    let mut index = 0usize;
    let mut previous = 0u32;

    while index < bytes.len() {
        let (delta, next_index) = decode_varint(bytes, index);
        let protein_id = if postings.is_empty() {
            delta
        } else {
            previous + delta
        };

        postings.push(protein_id);
        previous = protein_id;
        index = next_index;
    }

    postings
}

#[allow(dead_code)]
fn decode_varint(bytes: &[u8], mut index: usize) -> (u32, usize) {
    let mut shift = 0u32;
    let mut value = 0u32;

    loop {
        let byte = bytes[index];
        index += 1;
        value |= ((byte & 0x7f) as u32) << shift;

        if byte & 0x80 == 0 {
            return (value, index);
        }

        shift += 7;
    }
}
