#!/usr/bin/env python3
"""
verify_strand_differentiation.py
=================================
Utility script that (1) generates an antisense FASTA from a sense FASTA,
(2) sanity-checks that the two input files are genuinely distinct at the
sequence level, and (3) after PORT-EK find_k has been run in
--strand_specific mode, compares the resulting k-mer sets to confirm
sense-enriched and antisense-enriched k-mers are truly different.

Usage examples
--------------
# Step 1 – Generate antisense FASTA and check for palindromes:
python verify_strand_differentiation.py generate \
    --sense   hiv_sense.fasta \
    --antisense hiv_antisense.fasta

# Step 2 – After running PORT-EK find_k --strand_specific, compare sets:
python verify_strand_differentiation.py compare \
    --project_dir ./my_portek_project \
    --sense_group   sense \
    --antisense_group antisense \
    -k 15
"""

import argparse
import os
import pickle
import pathlib
import sys
from collections import Counter

# ── BioPython is already a PORT-EK dependency ──────────────────────────────
from Bio import SeqIO


# ============================================================================
# Helpers
# ============================================================================

_COMPLEMENT = str.maketrans("ACGTacgtNnRrYyKkMmSsWwBbDdHhVv",
                              "TGCAtgcaNnYyRrMmKkSsWwVvHhDdBb")


def _rc(seq: str) -> str:
    """Return the reverse complement of a nucleotide string."""
    return seq.translate(_COMPLEMENT)[::-1]


def _load_kmer_set(project_dir: str, k: int, group: str) -> set:
    """Load the pickled k-mer integer set for one group."""
    pkl_path = pathlib.Path(project_dir) / "input" / "indices" / f"{k}mer_{group}_set.pkl"
    if not pkl_path.exists():
        raise FileNotFoundError(
            f"Cannot find {pkl_path}.\n"
            f"Did you run:  portek find_k {project_dir} --strand_specific ?"
        )
    with open(pkl_path, "rb") as fh:
        return pickle.load(fh)


def _decode_kmer(kmer_id: int, k: int) -> str:
    """Convert a PORT-EK integer k-mer back to a nucleotide string."""
    decoding = {"00": "A", "01": "C", "10": "G", "11": "T"}
    bits = bin(kmer_id)[2:].rjust(2 * k, "0")
    return "".join(decoding[bits[i: i + 2]] for i in range(0, len(bits), 2))


# ============================================================================
# Sub-command: generate
# ============================================================================

def cmd_generate(args):
    """
    Write a reverse-complement FASTA (antisense) from a sense FASTA and
    report how many sequences are palindromes.
    """
    sense_path = args.sense
    antisense_path = args.antisense

    if not os.path.isfile(sense_path):
        sys.exit(f"[ERROR] Sense FASTA not found: {sense_path}")

    print(f"\n{'='*60}")
    print("STEP 1 — Generate antisense FASTA")
    print(f"{'='*60}")
    print(f"  Input  (sense)    : {sense_path}")
    print(f"  Output (antisense): {antisense_path}\n")

    records = list(SeqIO.parse(sense_path, "fasta"))
    if not records:
        sys.exit(f"[ERROR] No sequences found in {sense_path}")

    palindromes = []
    out_records = []
    seq_diffs = []

    for rec in records:
        sense_str  = str(rec.seq).upper()
        anti_str   = _rc(sense_str)
        is_palindrome = (sense_str == anti_str)
        if is_palindrome:
            palindromes.append(rec.id)

        # Compute Hamming-like difference (length-normalised)
        min_len = min(len(sense_str), len(anti_str))
        mismatches = sum(a != b for a, b in zip(sense_str[:min_len], anti_str[:min_len]))
        pct_diff = 100.0 * mismatches / min_len if min_len > 0 else 0.0
        seq_diffs.append(pct_diff)

        new_rec = rec.__class__(
            id=f"{rec.id}_RC",
            name=f"{rec.name}_RC",
            description=f"{rec.description} [reverse complement]",
            seq=rec.seq.__class__(anti_str),
        )
        out_records.append(new_rec)

    SeqIO.write(out_records, antisense_path, "fasta")

    # ── Report ─────────────────────────────────────────────────────────────
    n = len(records)
    avg_diff = sum(seq_diffs) / n if n else 0.0
    min_diff = min(seq_diffs) if seq_diffs else 0.0

    print(f"  Sequences processed : {n}")
    print(f"  Palindromes (sense == antisense) : {len(palindromes)}")
    if palindromes:
        print(f"    IDs: {', '.join(palindromes[:5])}"
              + (" ..." if len(palindromes) > 5 else ""))
    print(f"  Avg nucleotide difference (sense vs antisense) : {avg_diff:.1f}%")
    print(f"  Min nucleotide difference                       : {min_diff:.1f}%")

    print()
    if len(palindromes) == n:
        print("  [CRITICAL] ALL sequences are palindromes — the antisense FASTA")
        print("  is identical to the sense FASTA.  k-mer comparison will still")
        print("  return identical results even with --strand_specific.")
    elif palindromes:
        print(f"  [WARNING] {len(palindromes)}/{n} sequences are palindromes.")
        print("  Those sequences will contribute identical k-mers to both groups.")
    else:
        print("  [OK] No palindromes detected.")
        print("  Sense and antisense inputs are genuinely distinct.")

    print(f"\n  Antisense FASTA written to: {antisense_path}")
    print("\nNext steps:")
    print("  1. Copy the antisense FASTA into your project's input/ directory.")
    print("  2. Add it as a separate sample_group in config.yaml (e.g. 'antisense').")
    print("  3. Run:  portek find_k <project_dir> --strand_specific")
    print("  4. Run:  portek find_enriched <project_dir> -k <k>")
    print("  5. Run:  python verify_strand_differentiation.py compare ...")


# ============================================================================
# Sub-command: compare
# ============================================================================

def cmd_compare(args):
    """
    After PORT-EK find_k (--strand_specific), load both k-mer sets and
    report how many k-mers are shared vs. strand-exclusive.
    """
    project_dir    = args.project_dir
    sense_group    = args.sense_group
    antisense_group = args.antisense_group
    k              = args.k

    print(f"\n{'='*60}")
    print("STEP 2 — Compare sense vs antisense k-mer sets")
    print(f"{'='*60}")
    print(f"  Project dir      : {project_dir}")
    print(f"  Sense group      : {sense_group}")
    print(f"  Antisense group  : {antisense_group}")
    print(f"  k                : {k}\n")

    sense_set    = _load_kmer_set(project_dir, k, sense_group)
    antisense_set = _load_kmer_set(project_dir, k, antisense_group)

    shared          = sense_set & antisense_set
    sense_only      = sense_set - antisense_set
    antisense_only  = antisense_set - sense_set

    total_sense     = len(sense_set)
    total_antisense = len(antisense_set)
    total_union     = len(sense_set | antisense_set)

    print(f"  Total unique {k}-mers in sense group         : {total_sense:>10,}")
    print(f"  Total unique {k}-mers in antisense group     : {total_antisense:>10,}")
    print(f"  Union (all distinct k-mers)                  : {total_union:>10,}")
    print(f"  Shared (in BOTH groups)                      : {len(shared):>10,}  "
          f"({100*len(shared)/total_union:.1f}% of union)")
    print(f"  Sense-only k-mers                            : {len(sense_only):>10,}  "
          f"({100*len(sense_only)/total_union:.1f}% of union)")
    print(f"  Antisense-only k-mers                        : {len(antisense_only):>10,}  "
          f"({100*len(antisense_only)/total_union:.1f}% of union)")

    # ── Diagnosis ──────────────────────────────────────────────────────────
    print()
    if len(sense_only) == 0 and len(antisense_only) == 0:
        print("  [CRITICAL] The two k-mer sets are IDENTICAL.")
        print("  The most likely cause is that find_k was run WITHOUT --strand_specific.")
        print("  Re-run:  portek find_k <project_dir> --strand_specific")
        print("  Then re-run find_enriched.")
    elif len(shared) == total_union:
        # same as above, redundant guard
        print("  [CRITICAL] Sets are identical — see above.")
    else:
        overlap_pct = 100.0 * len(shared) / total_union
        if overlap_pct > 90:
            print(f"  [WARNING] {overlap_pct:.1f}% overlap — sets are very similar.")
            print("  Consider checking for palindromic sequences in your input.")
        else:
            print(f"  [OK] Sets differ meaningfully ({overlap_pct:.1f}% shared).")
            print("  Strand-specific enrichment analysis should yield distinct results.")

    # ── Optional: show example strand-exclusive k-mers ─────────────────────
    if args.show_examples and (sense_only or antisense_only):
        n_ex = min(args.show_examples, 10)
        print(f"\n  --- Up to {n_ex} sense-only k-mers ---")
        for kid in list(sense_only)[:n_ex]:
            print(f"    {_decode_kmer(kid, k)}")
        print(f"\n  --- Up to {n_ex} antisense-only k-mers ---")
        for kid in list(antisense_only)[:n_ex]:
            print(f"    {_decode_kmer(kid, k)}")

    # ── Check for reverse-complement pairs in the exclusive sets ───────────
    # A correct strand-specific run should show that most sense-only k-mers
    # have their RC counterpart in the antisense-only set.
    if args.rc_check and sense_only and antisense_only:
        print("\n  --- Reverse-complement pairing check ---")
        _COMPLEMENT_BITS = 0b11

        def _rc_kmer_id(kid, k):
            rc = 0
            for _ in range(k):
                rc = (rc << 2) | ((kid & 0b11) ^ _COMPLEMENT_BITS)
                kid >>= 2
            return rc

        matched = sum(
            1 for kid in sense_only if _rc_kmer_id(kid, k) in antisense_only
        )
        pct_matched = 100.0 * matched / len(sense_only) if sense_only else 0
        print(f"  Sense-only k-mers whose RC appears in antisense-only set : "
              f"{matched}/{len(sense_only)} ({pct_matched:.1f}%)")
        if pct_matched > 80:
            print("  [OK] High RC pairing confirms strand-specific indexing is working.")
        else:
            print("  [INFO] Lower RC pairing may indicate genuine sequence-level")
            print("         differences between isolates across strands.")


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Verify that sense and antisense PORT-EK inputs are genuinely "
            "differentiated, and compare the resulting k-mer sets after a "
            "--strand_specific PORT-EK run."
        )
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── generate ──────────────────────────────────────────────────────────
    gen = sub.add_parser(
        "generate",
        help="Generate an antisense FASTA from a sense FASTA and check for palindromes.",
    )
    gen.add_argument("--sense", required=True,
                     help="Path to the sense (original) FASTA file.")
    gen.add_argument("--antisense", required=True,
                     help="Destination path for the antisense FASTA file.")

    # ── compare ───────────────────────────────────────────────────────────
    cmp = sub.add_parser(
        "compare",
        help=(
            "After running PORT-EK find_k --strand_specific, compare the "
            "k-mer sets of the sense and antisense groups."
        ),
    )
    cmp.add_argument("--project_dir", required=True,
                     help="Path to the PORT-EK project directory.")
    cmp.add_argument("--sense_group", required=True,
                     help="Name of the sense sample group (must match config.yaml).")
    cmp.add_argument("--antisense_group", required=True,
                     help="Name of the antisense sample group (must match config.yaml).")
    cmp.add_argument("-k", type=int, required=True,
                     help="k-mer length to compare.")
    cmp.add_argument("--show_examples", type=int, default=5, metavar="N",
                     help="Print N example strand-exclusive k-mer sequences (default 5).")
    cmp.add_argument("--rc_check", action="store_true",
                     help=(
                         "Check what fraction of sense-only k-mers have their "
                         "reverse complement in the antisense-only set. "
                         "A high fraction confirms correct strand-specific indexing."
                     ))

    args = parser.parse_args()

    if args.command == "generate":
        cmd_generate(args)
    elif args.command == "compare":
        cmd_compare(args)


if __name__ == "__main__":
    main()
