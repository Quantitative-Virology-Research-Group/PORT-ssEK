#!/usr/bin/env python3
"""
revcomp_fasta.py
 
Reverse complement all sequences in a multi-FASTA file.
 
Usage:
    python revcomp_fasta.py input.fasta output.fasta
    python revcomp_fasta.py input.fasta            # writes to stdout
    cat input.fasta | python revcomp_fasta.py -     # read from stdin, write to stdout
"""
 
import sys
import argparse
 
# Standard complement table, including common IUPAC ambiguity codes
COMPLEMENT = str.maketrans(
    "ACGTUacgtuRYSWKMBDHVNryswkmbdhvn-",
    "TGCAAtgcaaYRSWMKVHDBNyrswmkvhdbn-"
)
 
 
def revcomp(seq: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""
    return seq.translate(COMPLEMENT)[::-1]
 
 
def parse_fasta(handle):
    """
    Yield (header, sequence) tuples from a FASTA file handle.
    Sequence lines are concatenated and whitespace is stripped.
    """
    header = None
    seq_chunks = []
 
    for line in handle:
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
            header = line[1:]  # drop the leading '>'
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
 
    if header is not None:
        yield header, "".join(seq_chunks)
 
 
def write_fasta_record(out, header, seq, wrap=70):
    """Write a single FASTA record, wrapping sequence lines at `wrap` chars."""
    out.write(f">{header}\n")
    if wrap and wrap > 0:
        for i in range(0, len(seq), wrap):
            out.write(seq[i:i + wrap] + "\n")
    else:
        out.write(seq + "\n")
 
 
def main():
    parser = argparse.ArgumentParser(
        description="Reverse complement all sequences in a multi-FASTA file."
    )
    parser.add_argument(
        "input",
        help="Input FASTA file path, or '-' to read from stdin"
    )
    parser.add_argument(
        "output",
        nargs="?",
        default=None,
        help="Output FASTA file path (default: stdout)"
    )
    parser.add_argument(
        "--wrap",
        type=int,
        default=70,
        help="Line width for wrapping output sequences (0 = no wrapping, default: 70)"
    )
    parser.add_argument(
        "--suffix",
        default="_revcomp",
        help="Suffix appended to each header (default: '_revcomp'; use '' for none)"
    )
    args = parser.parse_args()
 
    # Open input
    if args.input == "-":
        in_handle = sys.stdin
    else:
        in_handle = open(args.input, "r")
 
    # Open output
    if args.output:
        out_handle = open(args.output, "w")
    else:
        out_handle = sys.stdout
 
    try:
        count = 0
        for header, seq in parse_fasta(in_handle):
            rc_seq = revcomp(seq)
            new_header = f"{header}{args.suffix}" if args.suffix else header
            write_fasta_record(out_handle, new_header, rc_seq, wrap=args.wrap)
            count += 1
 
        if out_handle is not sys.stdout:
            print(f"Reverse complemented {count} sequence(s) -> {args.output}", file=sys.stderr)
    finally:
        if in_handle is not sys.stdin:
            in_handle.close()
        if out_handle is not sys.stdout:
            out_handle.close()
 
 
if __name__ == "__main__":
    main()
