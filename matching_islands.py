#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
matching_islands.py

Script for comparative analysis of genomic islands (GIPSy2) using
protein similarity (DIAMOND blastp) or nucleotide similarity (BLAST+ blastn).
Generates similarity matrices, networks, heatmaps, community detection
(Louvain), and, for proteins, shared/exclusive protein files.

Author: Janaíne Aparecida de Paula
Date: 17/08/2026
Version: best-hit protein averaging + reciprocal nucleotide HSP aggregation.
"""

from Bio import SeqIO
from pathlib import Path
from multiprocessing import Pool, cpu_count
import re
import csv
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
import time
from datetime import datetime
from scipy.cluster.hierarchy import dendrogram
import networkx as nx
import community as community_louvain
import subprocess
from collections import defaultdict

start_time = time.time()
analysis_time = datetime.now()

#  Biological comparison parameters:

# Protein: each query protein contributes to the average using the pident of
# its best significant hit in the target island; absence of a hit contributes
# 0%. The threshold below is used ONLY to define high-similarity homologs,
# shared proteins, and protein-content coverage.

PROTEIN_IDENTITY_THRESHOLD = 90.0
PROTEIN_CONTAINMENT_COVERAGE = 0.95

# Nucleotide: island conservation is evaluated using the average identity
# calculated from non-redundant HSPs and reciprocal island coverage.

NUCLEOTIDE_IDENTITY_THRESHOLD = 95.0
NUCLEOTIDE_CONTAINMENT_COVERAGE = 0.995
COVERAGE_EPSILON = 1e-6

# Common utility def

def read_fasta_sequences(faa_path):
    """Read a FASTA file and return a list of SeqRecord objects."""
    return list(SeqIO.parse(faa_path, "fasta"))

def read_fna_sequence(fna_path):
    """Read a .fna file containing a single sequence and return the sequence
    as a string. If multiple records are present, concatenate them."""
    records = list(SeqIO.parse(fna_path, "fasta"))
    if not records:
        return None
    seq = "".join(str(rec.seq) for rec in records)
    return seq

def prepare_all_island_fasta(faa_files, combined_fasta_path, mode):
    """
    Combine all protein or nucleotide sequences into a single FASTA file.
    For proteins, add a suffix to distinguish proteins from the same island.
    For nucleotides, use the island_id as the unique sequence ID.
    """
    with open(combined_fasta_path, 'w') as out:
        for island_id, file_path in faa_files.items():
            if mode == 'protein':
                records = SeqIO.parse(file_path, "fasta")
                for record in records:
                    record.id = f"{island_id}__{record.id}"
                    record.description = ""
                    SeqIO.write(record, out, "fasta")
            else:
                seq = read_fna_sequence(file_path)
                if seq is None:
                    continue
                from Bio.Seq import Seq
                from Bio.SeqRecord import SeqRecord
                record = SeqRecord(Seq(seq), id=island_id, description="")
                SeqIO.write(record, out, "fasta")


def run_alignment_allvsall(combined_fasta, db_name, output_file, mode, threads=8, n_targets=None):
    """Run DIAMOND blastp or BLASTn all-vs-all alignment."""
    if mode == 'protein':
        makedb_cmd = ["diamond", "makedb", "--in", str(combined_fasta), "-d", str(db_name)]
        subprocess.run(makedb_cmd, check=True)
        blast_cmd = [
            "diamond", "blastp",
            "-d", str(db_name),
            "-q", str(combined_fasta),
            "-o", str(output_file),
            "--threads", str(threads),
            "--evalue", "1e-5",
            "--max-target-seqs", "0",
            "--outfmt", "6", "qseqid", "sseqid", "pident", "nident", "length",
            "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        subprocess.run(blast_cmd, check=True)
    else:
        makeblastdb_cmd = [
            "makeblastdb", "-in", str(combined_fasta),
            "-dbtype", "nucl", "-out", str(db_name)
        ]
        subprocess.run(makeblastdb_cmd, check=True)
        if n_targets is None or n_targets < 1:
            raise ValueError("n_targets must be provided and >= 1 for BLASTn all-vs-all.")
        blast_cmd = [
            "blastn", "-db", str(db_name), "-query", str(combined_fasta),
            "-out", str(output_file), "-num_threads", str(threads),
            "-evalue", "1e-5", "-max_target_seqs", str(n_targets),
            "-outfmt",
            "6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
        ]
        subprocess.run(blast_cmd, check=True)
    return output_file


def parse_alignment_output(output_file, mode):
    """Read the tabular DIAMOND/BLASTn output."""
    hits = defaultdict(list)
    with open(output_file, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if mode == 'protein':
                (qseqid, sseqid, pident, nident, length, qlen, slen,
                 qstart, qend, sstart, send, evalue, bitscore) = parts
                hit = {
                    'target': sseqid, 'pident': float(pident), 'nident': int(nident),
                    'length': int(length), 'qlen': int(qlen), 'slen': int(slen),
                    'qstart': int(qstart), 'qend': int(qend),
                    'sstart': int(sstart), 'send': int(send),
                    'evalue': float(evalue), 'bitscore': float(bitscore)
                }
            else:
                (qseqid, sseqid, pident, nident, length, mismatch, gapopen,
                 qstart, qend, sstart, send, evalue, bitscore, qseq, sseq) = parts
                hit = {
                    'target': sseqid, 'pident': float(pident), 'nident': int(nident),
                    'length': int(length), 'mismatch': int(mismatch),
                    'gapopen': int(gapopen), 'evalue': float(evalue),
                    'bitscore': float(bitscore), 'qstart': int(qstart),
                    'qend': int(qend), 'sstart': int(sstart), 'send': int(send),
                    'qseq': qseq, 'sseq': sseq
                }
            hits[qseqid].append(hit)
    for q in hits:
        hits[q].sort(key=lambda x: x['bitscore'], reverse=True)
    return hits


def write_protein_comparison(seqs1, seqs2, high_matches, unmatched1, unmatched2,
                             output_directory, island_a, island_b):
    """Write FASTA files containing shared and exclusive proteins."""
    output_directory.mkdir(parents=True, exist_ok=True)
    seqs1_dict = {record.id: record for record in seqs1}
    seqs2_dict = {record.id: record for record in seqs2}
    shared_a, shared_b = [], []
    shared_pairs_file = output_directory / "shared_pairs.tsv"
    with open(shared_pairs_file, "w", encoding="utf-8") as out:
        out.write("Protein_A\tProtein_B\tIdentity\n")
        seen_a, seen_b = set(), set()
        for protein_a, protein_b, identity in high_matches:
            if protein_a in seqs1_dict and protein_a not in seen_a:
                shared_a.append(seqs1_dict[protein_a]); seen_a.add(protein_a)
            if protein_b in seqs2_dict and protein_b not in seen_b:
                shared_b.append(seqs2_dict[protein_b]); seen_b.add(protein_b)
            out.write(f"{protein_a}\t{protein_b}\t{identity*100:.2f}\n")
    exclusive_a = [seqs1_dict[pid] for pid in unmatched1 if pid in seqs1_dict]
    exclusive_b = [seqs2_dict[pid] for pid in unmatched2 if pid in seqs2_dict]
    SeqIO.write(shared_a, output_directory / "shared_A.faa", "fasta")
    SeqIO.write(shared_b, output_directory / "shared_B.faa", "fasta")
    SeqIO.write(exclusive_a, output_directory / "exclusive_A.faa", "fasta")
    SeqIO.write(exclusive_b, output_directory / "exclusive_B.faa", "fasta")


def generate_shared_protein_files(results, seq_files, island_output, island_type):
    """Generate shared and exclusive protein FASTA files only for island pairs containing at least one high-similarity shared protein."""

    protein_comparison_dir = island_output / "protein_comparison"
    protein_comparison_dir.mkdir(exist_ok=True)
    print(f"Generating shared and exclusive proteins for "f"{island_type} island pairs...")
    generated_pairs = 0
    for result in results: 
        (key1, key2, summary, avg_identity, cov_ab, cov_ba, n1, n2, high_matches, unmatched1, unmatched2) = result

        # Skip pairs without any shared protein meeting the identity threshold.
        if high_matches is None or not high_matches:
            continue

        island_a = simplify_matrix_id(key1)
        island_b = simplify_matrix_id(key2)
        pair_directory = (protein_comparison_dir / f"{island_a}__{island_b}")
        seqs1 = read_fasta_sequences(seq_files[key1])
        seqs2 = read_fasta_sequences(seq_files[key2])

        write_protein_comparison(seqs1, seqs2, high_matches, unmatched1, unmatched2, pair_directory, island_a, island_b)

        generated_pairs += 1

    if generated_pairs == 0:
        print("No island pairs with shared proteins above the "f"{PROTEIN_IDENTITY_THRESHOLD:.0f}% identity threshold were found.")
    else:
        print(f"Protein comparison files generated for "f"{generated_pairs} island pair(s).")
    return protein_comparison_dir


def _protein_hits_between(query_island, target_island, seqs_by_island, diamond_hits):
    """
        Compare each protein from query_island against ALL proteins from
        target_island and select exactly one best hit per query protein.

        Protein island similarity rule:
        - if any significant hit exists (already filtered by the DIAMOND
            E-value), use the pident of the BEST hit, regardless of whether
            it is >=90%, <90%, or any other reported value;
        - if no significant hit exists in the target island, the protein
            contributes 0%;
        - A->B similarity is the average of these values across ALL
            proteins in island A.

        The 90% threshold does NOT filter the average. It is applied only
        afterward to identify high-similarity homologs used for --shared
        and protein-content coverage.
        """

    seqs_query = seqs_by_island[query_island]
    best_matches = []
    high_query = set()
    high_target = set()
    identities = []

    for record in seqs_query:
        query_composed = f"{query_island}__{record.id}"
        candidate_hits = []

        for hit in diamond_hits.get(query_composed, []):
            target_composed = hit['target']
            if '__' not in target_composed:
                continue
            target_prefix, target_original = target_composed.rsplit('__', 1)
            if target_prefix != target_island:
                continue
            candidate_hits.append((target_original, hit))

        # No significant hit between this protein and the target island:
        # explicitly contribute 0% to the query average.
        if not candidate_hits:
            identities.append(0.0)
            continue

        # Select exactly one best hit per query protein.
        target_original, best_hit = max(candidate_hits, key=lambda item: item[1]['bitscore'])
        identity = float(best_hit['pident'])

        identities.append(identity)
        best_matches.append((record.id, target_original, identity / 100.0))

        # The 90% threshold is used ONLY to define high-similarity
        # homologs / shared proteins and does not affect the average.
        if identity >= PROTEIN_IDENTITY_THRESHOLD:
            high_query.add(record.id)
            high_target.add(target_original)

    avg_identity = sum(identities) / len(seqs_query) if seqs_query else 0.0
    return best_matches, high_query, high_target, avg_identity


def process_pair_protein(args):
    """Process a pair of islands in both directions in protein mode."""
    key1, key2, faa_files, diamond_hits = args

    seqs_by_island = {
        key1: read_fasta_sequences(faa_files[key1]),
        key2: read_fasta_sequences(faa_files[key2])
    }
    total1 = len(seqs_by_island[key1])
    total2 = len(seqs_by_island[key2])

    # A -> B
    matches_ab, matched_a, matched_b_from_ab, identity_ab = _protein_hits_between(
        key1, key2, seqs_by_island, diamond_hits
    )

    # B -> A (required for reciprocal coverage/identity)
    matches_ba, matched_b, matched_a_from_ba, identity_ba = _protein_hits_between(
        key2, key1, seqs_by_island, diamond_hits
    )

    unmatched1 = [r.id for r in seqs_by_island[key1] if r.id not in matched_a]
    unmatched2 = [r.id for r in seqs_by_island[key2] if r.id not in matched_b]

    coverage_1_to_2 = len(matched_a) / total1 if total1 else 0.0
    coverage_2_to_1 = len(matched_b) / total2 if total2 else 0.0

    ## Island similarity: average of the mean identities in both directions.
    ## The matrix is therefore symmetric, while reciprocal coverage is retained.
    avg_identity = (identity_ab + identity_ba) / 2.0

    # --shared uses only best hits exceeding the high-similarity threshold.
    # In contrast, the matrix uses ALL best hits, including identities
    # below 90% and proteins without hits (0%).
    shared_pairs = {}
    for a, b, identity in matches_ab:
        if identity * 100.0 >= PROTEIN_IDENTITY_THRESHOLD:
            shared_pairs[(a, b)] = identity
    for b, a, identity in matches_ba:
        if identity * 100.0 >= PROTEIN_IDENTITY_THRESHOLD:
            shared_pairs.setdefault((a, b), identity)
    high_matches = [(a, b, identity) for (a, b), identity in sorted(shared_pairs.items())]

    summary = (
        f"{key1} vs {key2}: average similarity of {avg_identity:.2f}%.\n"
        f"- Identity A→B: {identity_ab:.2f}% | coverage A→B: {coverage_1_to_2*100:.1f}%\n"
        f"- Identity B→A: {identity_ba:.2f}% | coverage B→A: {coverage_2_to_1*100:.1f}%\n"
        f"- {len(unmatched1)} protein(s) from {key1} without a qualifying homolog in {key2}\n"
        f"- {len(unmatched2)} protein(s) from {key2} without a qualifying homolog in {key1}"
    )

    return (key1, key2, summary, avg_identity,
            coverage_1_to_2, coverage_2_to_1,
            total1, total2, high_matches, unmatched1, unmatched2)

def _normal_interval(start, end):
    """Return an interval with its coordinates in ascending order."""
    return (min(start, end), max(start, end))


def _interval_length(interval):
    """Return the length of an inclusive genomic interval."""
    return interval[1] - interval[0] + 1


def _merge_position_sets(intervals):
    """Return the number of unique positions covered by the intervals."""
    if not intervals:
        return 0
    positions = set()
    for start, end in intervals:
        positions.update(range(start, end + 1))
    return len(positions)


def _select_best_nonredundant_hsps(hsps):
    """
        Sort HSPs by bit score.

        HSPs that add at least one new position to the query and subject are
        retained. Overlapping regions are resolved afterward by query position,
        keeping the HSP with the highest bit score.
        """
    return sorted(hsps, key=lambda h: h['bitscore'], reverse=True)


def _build_nucleotide_position_map(hsps):
    """
        Build exact query and subject position maps from qseq/sseq.

        For each query position, only the best HSP covering that position is
        retained. This prevents double counting of overlapping regions and
        allows calculation of the actual base identity across covered positions.
        """

    query_map = {}
    subject_map = {}

    for hsp in _select_best_nonredundant_hsps(hsps):
        qpos = hsp['qstart']
        spos = hsp['sstart']
        q_step = 1 if hsp['qend'] >= hsp['qstart'] else -1
        s_step = 1 if hsp['send'] >= hsp['sstart'] else -1

        for qbase, sbase in zip(hsp['qseq'], hsp['sseq']):
            q_is_base = qbase != '-'
            s_is_base = sbase != '-'

            current_q = qpos if q_is_base else None
            current_s = spos if s_is_base else None

            if current_q is not None and current_s is not None:
                if current_q not in query_map:
                    identical = qbase.upper() == sbase.upper()
                    query_map[current_q] = (current_s, identical, hsp['bitscore'])
                    subject_map[current_s] = (current_q, identical, hsp['bitscore'])

            if q_is_base:
                qpos += q_step
            if s_is_base:
                spos += s_step

    return query_map, subject_map


def _nucleotide_direction_metrics(query_id, target_id, blast_hits, seq_lengths):
    """
        Calculate identity and coverage for query -> subject without
        double-counting overlapping HSPs.

        Identity is calculated directly from HSP bases (qseq/sseq),
        retaining only the best evidence for each query position.

        Coverage is the number of distinct positions covered in the
        query and subject.
        """
    pair_hits = [h for h in blast_hits.get(query_id, []) if h['target'] == target_id]
    if not pair_hits:
        return 0.0, 0.0, 0.0, []

    query_map, subject_map = _build_nucleotide_position_map(pair_hits)
    if not query_map:
        return 0.0, 0.0, 0.0, []

    total_aligned = len(query_map)
    total_identical = sum(1 for _, (_, identical, _) in query_map.items() if identical) 

    q_len = seq_lengths[query_id]
    s_len = seq_lengths[target_id]
    identity = total_identical / q_len * 100.0 #q_len = total query sequence length
    cov_query = total_aligned / q_len if q_len else 0.0
    cov_subject = len(subject_map) / s_len if s_len else 0.0

    return identity, cov_query, cov_subject, pair_hits


def process_pair_nucleotide(args):
    """
        Process a pair of islands in both directions in nucleotide mode.

        Identity is calculated from the total number of identical bases in
        non-redundant HSPs, normalized by the full query length.

        Coverage is calculated as the union of aligned regions, avoiding
        double counting of overlapping HSPs.
        """
    key1, key2, fna_files, blast_hits, seq_lengths = args
    len1 = seq_lengths[key1]
    len2 = seq_lengths[key2]

    identity_ab, cov1, cov2_from_ab, hsps_ab = _nucleotide_direction_metrics(
        key1, key2, blast_hits, seq_lengths
    )
    identity_ba, cov2, cov1_from_ba, hsps_ba = _nucleotide_direction_metrics(
        key2, key1, blast_hits, seq_lengths
    )

    # For a reciprocal comparison, use both directions.
    avg_identity = (identity_ab + identity_ba) / 2.0

    summary = (
        f"{key1} vs {key2}: nucleotide identity = {avg_identity:.2f}%, "
        f"coverage A→B = {cov1*100:.1f}%, coverage B→A = {cov2*100:.1f}%\n"
        f"- Identity A→B: {identity_ab:.2f}% ({len(hsps_ab)} non-redundant HSPs)\n"
        f"- Identity B→A: {identity_ba:.2f}% ({len(hsps_ba)} non-redundant HSPs)"
    )

    return (key1, key2, summary, avg_identity,
        cov1, cov2,
        len1, len2,
        identity_ab, identity_ba,
        None)


# ==================
#  Similarity network extraction and simplification def's

def extract_similarity_network(summary_file, output_csv, mode):
    """
        Extract network edges from the summary file.

        The pattern is adjusted according to whether the analysis is performed
        at the protein or nucleotide level.
        """

    if mode == 'protein':
        pattern = r"(.+\.faa) vs (.+\.faa): average similarity of ([\d\.]+)%"
    else:
        pattern = r"(.+\.fna) vs (.+\.fna): nucleotide identity = ([\d\.]+)%, coverage"
    with open(summary_file, "r", encoding="utf-8") as infile, \
         open(output_csv, "w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Source", "Target", "Weight", "Group"])
        for line in infile:
            match = re.search(pattern, line)
            if match:
                source = match.group(1).strip()
                target = match.group(2).strip()
                weight = float(match.group(3))
                source_prefix = Path(source).parts[0]
                target_prefix = Path(target).parts[0]
                if source_prefix == target_prefix:
                    continue
                high_threshold = PROTEIN_IDENTITY_THRESHOLD if mode == 'protein' else NUCLEOTIDE_IDENTITY_THRESHOLD
                if weight >= high_threshold:
                    group = f"High (>={high_threshold:.0f}%)"
                elif 50 <= weight < high_threshold:
                    group = f"Medium (50-{high_threshold-1:.0f}%)"
                else:
                    group = "Low (<50%)"
                writer.writerow([source, target, weight, group])

def rename_network_nodes(input_csv):
    """Rename network nodes using simplify_matrix_id."""
    input_path = Path(input_csv)
    output_path = input_path.with_name("renamed_" + input_path.stem + ".csv")
    with open(input_path, newline='', encoding="utf-8") as infile, \
         open(output_path, "w", newline='', encoding="utf-8") as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            row["Source"] = simplify_matrix_id(row["Source"])
            row["Target"] = simplify_matrix_id(row["Target"])
            writer.writerow(row)
    return output_path

def simplify_matrix_id(name):
    """Simplify an island identifier to a short format (e.g., strainA_PI1)."""
    path = Path(name)
    strain = path.parts[0].replace("_Islands_fisher", "")
    filename = path.name.replace(".faa", "").replace(".fna", "")
    island = re.search(r"(Virulence|Resistance|Symbiosis|Metabolic)_Island_(\d+)", filename)
    if not island:
        return name
    island_type = island.group(1)
    number = island.group(2)
    abbreviations = {
        "Virulence": "PI",
        "Resistance": "RI",
        "Symbiosis": "SI",
        "Metabolic": "MI"
    }
    return f"{strain}_{abbreviations[island_type]}{number}"

# ===================
#  Relationship classification (protein and nucleotide)


def classify_relationship_protein(cov_ab, cov_ba):
    """Classify the relationship between two islands based on protein coverage."""
    if cov_ab <= COVERAGE_EPSILON and cov_ba <= COVERAGE_EPSILON:
        return "No shared proteins"
    if cov_ab >= PROTEIN_CONTAINMENT_COVERAGE and cov_ba >= PROTEIN_CONTAINMENT_COVERAGE:
        return "Equivalent islands"
    if abs(cov_ab - 1) < COVERAGE_EPSILON and cov_ba < 1 - COVERAGE_EPSILON:
        return "A is contained in B"
    if abs(cov_ba - 1) < COVERAGE_EPSILON and cov_ab < 1 - COVERAGE_EPSILON:
        return "B is contained in A"
    return "Partial overlap"


def classify_relationship_nucleotide(cov_ab, cov_ba, identity_ab, identity_ba, id_threshold=NUCLEOTIDE_IDENTITY_THRESHOLD, eps=COVERAGE_EPSILON):
    """Classify the relationship between two islands based on nucleotide identity and coverage."""
    if cov_ab <= eps and cov_ba <= eps:
        return "No shared sequence"

    full_cov = NUCLEOTIDE_CONTAINMENT_COVERAGE
    if cov_ab >= full_cov and cov_ba >= full_cov and identity_ab >= id_threshold and identity_ba >= id_threshold:
        return "Equivalent islands"

    if cov_ab >= full_cov and identity_ab >= id_threshold and cov_ba < full_cov:
        return "A is contained in B"

    if cov_ba >= full_cov and identity_ba >= id_threshold and cov_ab < full_cov:
        return "B is contained in A"

    return "Partial overlap"


# ====================
#  Similarity matrix, heatmap, and dendrogram

def generate_similarity_matrix(results, keys, output_prefix): 
    """Generate and save the pairwise similarity matrix."""
    matrix = pd.DataFrame(100.0, index=keys, columns=keys)
    for key1, key2, _, avg_identity, *_ in results:
        matrix.loc[key1, key2] = avg_identity
        matrix.loc[key2, key1] = avg_identity
    display_names = {key: simplify_matrix_id(key) for key in keys}
    matrix.rename(index=display_names, columns=display_names, inplace=True)
    matrix.to_csv(f"{output_prefix}_similarity_matrix.tsv", sep="\t")
    return matrix

def linkage_to_newick(linkage_matrix, labels):
    """Convert a hierarchical clustering linkage matrix to Newick format."""
    tree = to_tree(linkage_matrix, rd=False)
    def build_newick(node):
        if node.is_leaf():
            return labels[node.id]
        left = build_newick(node.left)
        right = build_newick(node.right)
        return f"({left},{right}):{node.dist:.6f}"
    return build_newick(tree) + ";"

def generate_similarity_heatmap(matrix, output_prefix, max_cells_with_text=100):
    """Generate a similarity heatmap with dendrograms. Figure dimensions and font sizes are adjusted dynamically according to the number of islands."""

    distance_matrix = 1 - (matrix / 100)
    condensed = squareform(distance_matrix.values, checks=False)
    linkage_matrix = linkage(condensed, method="average")
    newick_labels = list(matrix.index)
    newick_tree = linkage_to_newick(linkage_matrix, newick_labels)
    newick_file = Path(f"{output_prefix}_dendrogram.newick")
    with open(newick_file, "w", encoding="utf-8") as f:
        f.write(newick_tree)
    dendro = dendrogram(linkage_matrix, no_plot=True)
    order = dendro["leaves"]
    ordered_matrix = matrix.iloc[order, order]

    n = len(ordered_matrix)
    figsize = max(12, min(80, 12 + n * 0.15))
    figsize = (figsize, figsize)
    fig = plt.figure(figsize=figsize)
    if n > 200:
        width_ratio = [0.5, 9.5]
    else:
        width_ratio = [1.5, 8]
    gs = fig.add_gridspec(2, 2, width_ratios=width_ratio, height_ratios=width_ratio,
                          wspace=0.02, hspace=0.02)
    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_heat = fig.add_subplot(gs[1, 1])
    dendrogram(linkage_matrix, orientation="top", no_labels=True,
               color_threshold=0, above_threshold_color="black", ax=ax_top)
    dendrogram(linkage_matrix, orientation="left", no_labels=True,
               color_threshold=0, above_threshold_color="black", ax=ax_left)
    ax_left.invert_yaxis()
    ax_top.axis("off")
    ax_left.axis("off")
    for spine in ax_top.spines.values():
        spine.set_visible(False)
    for spine in ax_left.spines.values():
        spine.set_visible(False)
    color_map = LinearSegmentedColormap.from_list(
        "match_islands",
        ["#FFF700", "#3DCB7D", "#076285"]
    )
    image = ax_heat.imshow(ordered_matrix, aspect="auto", interpolation="nearest",
                           cmap=color_map, vmin=0, vmax=100)

    if n <= max_cells_with_text:
        if n <= 30:
            cell_fontsize = 12
        elif n <= 45:
            cell_fontsize = 9
        elif n <= 80:
            cell_fontsize = 7
        else:
            cell_fontsize = 5
        for i in range(n):
            for j in range(n):
                ax_heat.text(j, i, f"{ordered_matrix.iloc[i, j]:.1f}",
                             ha="center", va="center", color="black",
                             fontsize=cell_fontsize)

    ax_heat.set_xticks(np.arange(n))
    ax_heat.set_yticks(np.arange(n))
    if n <= 30:
        label_font = 15
    elif n <= 45:
        label_font = 10
    elif n <= 80:
        label_font = 8
    elif n <= 150:
        label_font = 6
    elif n <= 300:
        label_font = 4
    else:
        label_font = 2
    ax_heat.set_xticklabels(ordered_matrix.columns, fontsize=label_font, rotation=90)
    ax_heat.set_yticklabels(ordered_matrix.index, fontsize=label_font)
    ax_heat.yaxis.tick_right()
    ax_heat.yaxis.set_label_position("right")
    ax_heat.tick_params(axis="y", left=False, labelleft=False,
                        right=True, labelright=True)
    cax = fig.add_axes([0.04, 0.15, 0.02, 0.60])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label("Similarity (%)")
    plt.subplots_adjust(left=0.10, right=0.95, bottom=0.05, top=0.98,
                        wspace=0.02, hspace=0.02)
    dpi = 300 if n > 100 else 100
    plt.savefig(f"{output_prefix}_heatmap.pdf", dpi=dpi)
    plt.close()

# ====================
# Global summary

def generate_global_summary(base_dir, renamed_network, global_summary,
                            island_pattern, island_type, analysis_start,
                            start_time, cpu_used, mode):
    """Generate a global summary of the genomic island analysis."""
    elapsed = time.time() - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = int(elapsed % 60)
    base_dir = Path(base_dir)
    ext = ".faa" if mode == 'protein' else ".fna"
    islands_per_strain = {}
    all_strains = []
    for folder in base_dir.iterdir():
        if folder.is_dir():
            all_strains.append(folder.name)
            if mode == 'protein':
                seq_dir = folder / "Amino_acids"
            else:
                seq_dir = folder / "Islands_nucleotides"
            if seq_dir.exists():
                islands = list(seq_dir.glob(island_pattern.replace(".faa", ext).replace(".fna", ext)))
                if islands:
                    islands_per_strain[folder.name] = len(islands)
    strains_without_islands = [strain for strain in all_strains if strain not in islands_per_strain]
    similarity_counts = {"High": 0, "Medium": 0, "Low": 0}
    total_edges = 0
    with open(renamed_network, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_edges += 1
            if "High" in row["Group"]:
                similarity_counts["High"] += 1
            elif "Medium" in row["Group"]:
                similarity_counts["Medium"] += 1
            elif "Low" in row["Group"]:
                similarity_counts["Low"] += 1
    unit = "proteins" if mode == 'protein' else "nucleotides (bp)"
    with open(global_summary, "w", encoding="utf-8") as out:
        out.write(f"GLOBAL SUMMARY OF {island_type.upper()} ISLANDS ANALYSIS ({mode.upper()})\n")
        out.write("=" * 60 + "\n\n")
        out.write(f"Analysis date: {analysis_start.strftime('%Y-%m-%d %H:%M:%S')}\n")
        out.write(f"Island category: {island_type}\n")
        out.write(f"Sequence level: {mode}\n")
        out.write(f"Execution time: {hours:02d}:{minutes:02d}:{seconds:02d}\n")
        out.write(f"CPU cores used: {cpu_used}\n\n")
        out.write(f"Total number of strains analyzed: {len(all_strains)}\n")
        out.write(f"Strains with detected {island_type.lower()} islands: {len(islands_per_strain)}\n")
        out.write(f"Strains without detected {island_type.lower()} islands: {len(strains_without_islands)}\n\n")
        out.write(f"{island_type} islands per strain:\n")
        for strain, count in sorted(islands_per_strain.items()):
            out.write(f"- {strain}: {count} islands\n")
        out.write(f"\nStrains without detected {island_type.lower()} islands:\n")
        for strain in sorted(strains_without_islands):
            out.write(f"- {strain}\n")
        high_threshold = PROTEIN_IDENTITY_THRESHOLD if mode == 'protein' else NUCLEOTIDE_IDENTITY_THRESHOLD
        medium_upper = high_threshold - 1
        out.write("\nSimilarity network summary:\n")
        out.write(f"- Total similarity edges: {total_edges}\n")
        out.write(f"- High similarity (>={high_threshold:.0f}%): {similarity_counts['High']}\n")
        out.write(f"- Medium similarity (50-{medium_upper:.0f}%): {similarity_counts['Medium']}\n")
        out.write(f"- Low similarity (<50%): {similarity_counts['Low']}\n")

# ====================
#  Louvain clustering

def run_louvain_clustering(network_file, output_directory):
    """Perform Louvain community detection on the similarity network."""
    df = pd.read_csv(network_file)
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row["Source"], row["Target"], weight=row["Weight"])
    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        print("No network available or empty.")
        return
    print(f"Running Louvain clustering on {G.number_of_nodes()} nodes and {G.number_of_edges()} edges...")
    partition = community_louvain.best_partition(G, weight="weight", resolution=1.0, random_state=42)
    nodes_df = (pd.DataFrame.from_dict(partition, orient="index", columns=["Community"])
                .reset_index().rename(columns={"index": "Node"}))
    nodes_df.to_csv(output_directory / "island_communities.tsv", sep="\t", index=False)
    community_stats = (nodes_df.groupby("Community")
                       .agg(Number_of_islands=("Node", "count"),
                            Islands=("Node", lambda x: "; ".join(sorted(x))))
                       .reset_index())
    community_stats.to_csv(output_directory / "community_statistics.tsv", sep="\t", index=False)
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    n_communities = nodes_df["Community"].nunique()
    largest = community_stats["Number_of_islands"].max()
    smallest = community_stats["Number_of_islands"].min()
    with open(output_directory / "louvain_summary.txt", "w", encoding="utf-8") as out:
        out.write("LOUVAIN COMMUNITY DETECTION SUMMARY\n")
        out.write("=" * 40 + "\n\n")
        out.write(f"Total nodes: {n_nodes}\n")
        out.write(f"Total edges: {n_edges}\n")
        out.write(f"Communities detected: {n_communities}\n")
        out.write(f"Largest community: {largest} islands\n")
        out.write(f"Smallest community: {smallest} islands\n\n")
        out.write("Community sizes:\n")
        for _, row in community_stats.iterrows():
            out.write(f"- Community {row['Community']}: {row['Number_of_islands']} islands\n")
    print("Louvain clustering completed successfully.")
    print(f"- Results directory: {output_directory}")

# ====================
#  Main

def main():
    parser = argparse.ArgumentParser(
        description="Comparative analysis of genomic islands (protein or nucleotide) predicted by GIPSy2."
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="Directory containing the GIPSy2 output folders.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--protein", action="store_true",
                       help="Protein-level analysis (uses .faa files in Amino_acids/).")
    group.add_argument("-n", "--nucleotide", action="store_true",
                       help="Nucleotide-level analysis (uses .fna files in Islands_nucleotides/).")
    parser.add_argument("-vir", "--virulence", action="store_true",
                        help="Analyze virulence/pathogenicity islands.")
    parser.add_argument("-res", "--resistance", action="store_true",
                        help="Analyze resistance islands.")
    parser.add_argument("-sym", "--symbiosis", action="store_true",
                        help="Analyze symbiotic islands.")
    parser.add_argument("-met", "--metabolic", action="store_true",
                        help="Analyze metabolic islands.")
    parser.add_argument("--all", action="store_true",
                        help="Analyze all genomic island types.")
    parser.add_argument("--plot", action="store_true",
                        help="Generate similarity heatmap.")
    parser.add_argument("--louvain", action="store_true",
                        help="Run Louvain community detection on the similarity network.")
    parser.add_argument("--shared", action="store_true",
                        help="Generate shared and exclusive protein .faa files (protein mode only).")
    parser.add_argument("--from-results", action="store_true",
                        help="Reuse previously generated outputs instead of recomputing pairwise comparisons.")
    args = parser.parse_args()

    analysis_start = datetime.now()
    start_time = time.time()

    if args.protein:
        mode = 'protein'
        seq_dir_name = "Amino_acids"
        file_ext = ".faa"
        process_pair_func = process_pair_protein
        classify_rel_func = classify_relationship_protein
        shared_allowed = True
    else:
        mode = 'nucleotide'
        seq_dir_name = "Islands_nucleotides"
        file_ext = ".fna"
        process_pair_func = process_pair_nucleotide
        classify_rel_func = classify_relationship_nucleotide
        shared_allowed = False

    if args.shared and not shared_allowed:
        print("WARNING: --shared option is only available for protein mode (-p). Ignoring.")
        args.shared = False

    selected_islands = []
    if args.all:
        selected_islands = [
            ("Virulence", f"Virulence_Island_*{file_ext}"),
            ("Resistance", f"Resistance_Island_*{file_ext}"),
            ("Symbiosis", f"Symbiosis_Island_*{file_ext}"),
            ("Metabolic", f"Metabolic_Island_*{file_ext}"),
        ]
    else:
        if args.virulence:
            selected_islands.append(("Virulence", f"Virulence_Island_*{file_ext}"))
        if args.resistance:
            selected_islands.append(("Resistance", f"Resistance_Island_*{file_ext}"))
        if args.symbiosis:
            selected_islands.append(("Symbiosis", f"Symbiosis_Island_*{file_ext}"))
        if args.metabolic:
            selected_islands.append(("Metabolic", f"Metabolic_Island_*{file_ext}"))

    if not selected_islands:
        parser.error("Please specify at least one island category (-vir, -res, -sym, -met or --all).")

    base_dir = args.input.resolve()
    if not base_dir.exists() or not base_dir.is_dir():
        parser.error(f"Input directory not found or not a directory: {base_dir}")

    suffix = "_p" if mode == 'protein' else "_n"
    results_dir = base_dir / f"results_matching_islands{suffix}"
    results_dir.mkdir(exist_ok=True)

    for island_type, island_pattern in selected_islands:
        island_output = results_dir / island_type.lower()
        island_output.mkdir(exist_ok=True)
        print(f"\nAnalyzing {island_type} islands ({mode} mode)...")

        seq_files = {}
        for folder in base_dir.iterdir():
            seq_dir = folder / seq_dir_name
            if seq_dir.exists():
                for seq_file in seq_dir.glob(island_pattern):
                    seq_files[f"{folder.name}/{seq_file.name}"] = seq_file

        if len(seq_files) < 2:
            print(f"Skipping {island_type}: fewer than two islands found.")
            continue

        summary_file = island_output / f"{island_type.lower()}_comparison_results_summary.txt"
        network_file = island_output / f"{island_type.lower()}_cytoscape_network_classified.csv"
        renamed_network = island_output / f"renamed_{island_type.lower()}_cytoscape_network_classified.csv"
        global_summary = island_output / f"{island_type.lower()}_global_analysis_summary.txt"
        coverage_file = island_output / f"{island_type.lower()}_island_pairwise_coverage.csv"
        matrix_file = island_output / f"{island_type.lower()}_similarity_matrix.tsv"

        if args.from_results:
            print(f"Using existing results for {island_type}...")
            if args.shared and mode == "protein":
                blast_output = island_output / "diamond_allvsall.tsv"
                if not blast_output.exists():
                    parser.error(f"Existing DIAMOND alignment file not found:\n"f"{blast_output}\n\n""The --from-results --shared combination requires ""diamond_allvsall.tsv generated during a previous analysis.")
                print("Loading existing DIAMOND all-vs-all alignment results...")
                alignment_hits = parse_alignment_output(blast_output, mode)
                keys = sorted(seq_files.keys())
                tasks = [
                    (
                        keys[i],
                        keys[j],
                        seq_files,
                        alignment_hits
                    )
                    for i in range(len(keys))
                    for j in range(i + 1, len(keys))
                ]
                print(f"Reconstructing protein comparisons using "f"{cpu_count()} CPU cores...")
                print(f"Total pairwise comparisons: {len(tasks)}")
                with Pool(cpu_count()) as pool:
                    results = pool.map(process_pair_protein,tasks)
                generate_shared_protein_files(results, seq_files, island_output, island_type)
                print(f"Shared and exclusive protein files generated successfully.")

            if args.plot:
                if not matrix_file.exists():
                    parser.error(f"Similarity matrix not found:\n{matrix_file}")
                similarity_matrix = pd.read_csv(matrix_file, sep="\t", index_col=0)
                generate_similarity_heatmap(similarity_matrix, island_output / island_type.lower())

            if args.louvain:
                if not renamed_network.exists():
                    parser.error(f"Network file not found:\n{renamed_network}")
                louvain_dir = island_output / "louvain_clustering"
                louvain_dir.mkdir(exist_ok=True)
                run_louvain_clustering(renamed_network, louvain_dir)
            continue

        combined_fasta = island_output / "all_sequences.fasta"
        db_name = island_output / ("diamond_db" if mode == "protein" else "blastn_db")
        blast_output = island_output / ("diamond_allvsall.tsv" if mode == "protein" else "blastn_allvsall.tsv")

        prepare_all_island_fasta(seq_files, combined_fasta, mode)
        print("Running all-vs-all alignment...")
        run_alignment_allvsall(
            combined_fasta, db_name, blast_output, mode,
            threads=cpu_count(), n_targets=len(seq_files)
        )
        print("Parsing alignment output...")
        alignment_hits = parse_alignment_output(blast_output, mode)

        keys = sorted(seq_files.keys())
        if mode == 'protein':
            tasks = [(keys[i], keys[j], seq_files, alignment_hits)
                     for i in range(len(keys)) for j in range(i+1, len(keys))]
        else:
            seq_lengths = {}
            for key, fpath in seq_files.items():
                seq = read_fna_sequence(fpath)
                seq_lengths[key] = len(seq) if seq else 0
            tasks = [(keys[i], keys[j], seq_files, alignment_hits, seq_lengths)
                     for i in range(len(keys)) for j in range(i+1, len(keys))]

        print(f"Starting parallel processing using {cpu_count()} CPU cores")
        print(f"Total pairwise comparisons: {len(tasks)}")

        with Pool(cpu_count()) as pool:
            results = pool.map(process_pair_func, tasks)

        if args.shared and mode == 'protein':
            generate_shared_protein_files(results, seq_files, island_output, island_type)

        with open(summary_file, "w", encoding="utf-8") as f:
            for res in results:
                f.write(res[2] + "\n\n")

        with open(coverage_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            if mode == 'protein':
                writer.writerow([
                    "Island_A", "Island_B", "Total_Genes_A", "Total_Genes_B",
                    "Average_Identity", "Coverage_A_to_B", "Coverage_B_to_A", "Relationship"
                ])
            else:
                writer.writerow([
                    "Island_A", "Island_B", "Total_Length_A", "Total_Length_B",
                    "Average_Identity", "Coverage_A_to_B", "Coverage_B_to_A", "Relationship"
                ])
            for res in results:
                key1, key2, _, avg_identity, cov_ab, cov_ba, n1, n2, identity_ab, identity_ba, _ = res

                if mode == 'protein':
                    rel = classify_rel_func(cov_ab, cov_ba)
                else:
                    rel = classify_rel_func(
                        cov_ab,
                        cov_ba,
                        identity_ab,
                        identity_ba
                    )
                writer.writerow([
                    simplify_matrix_id(key1), simplify_matrix_id(key2),
                    n1, n2, f"{avg_identity:.3f}",
                    f"{cov_ab:.6f}", f"{cov_ba:.6f}", rel
                ])

        extract_similarity_network(summary_file, network_file, mode)
        renamed_network = rename_network_nodes(network_file)
        if network_file.exists():
            network_file.unlink()

        if args.louvain:
            louvain_dir = island_output / "louvain_clustering"
            louvain_dir.mkdir(exist_ok=True)
            run_louvain_clustering(renamed_network, louvain_dir)

        similarity_matrix = generate_similarity_matrix(results, keys, island_output / island_type.lower())
        generate_global_summary(base_dir, renamed_network, global_summary,
                                island_pattern, island_type, analysis_start,
                                start_time, cpu_count(), mode)

        if args.plot:
            generate_similarity_heatmap(similarity_matrix, island_output / island_type.lower())


        print(f"\n{island_type} analysis completed.")
        print("Generated files:")
        print(f"- {summary_file}")
        print(f"- {renamed_network}")
        print(f"- {global_summary}")
        print(f"- {coverage_file}")
        print(f"- {island_type.lower()}_similarity_matrix.tsv")
        if args.plot:
            print(f"- {island_type.lower()}_heatmap.pdf")
            print(f"- {island_type.lower()}_dendrogram.newick")
        if args.louvain:
            print(f"- {louvain_dir}/island_communities.tsv")
            print(f"- {louvain_dir}/community_statistics.tsv")
        if args.shared and mode == 'protein':
            print(f"- {island_output / 'protein_comparison'}")

    print("\nPipeline completed successfully.")
    print("Thanks for using matching islands!")

if __name__ == "__main__":
    main()
