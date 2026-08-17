import argparse
import pickle
from collections import defaultdict

from Bio import SeqIO
from datetime import datetime


parser = argparse.ArgumentParser(
    description="Find all kmers in all sequences from a fasta file. Save their counts and positions."
)
parser.add_argument("in_file", help="path to the .fasta file with input sequences")
parser.add_argument(
    "out_dir",
    help="absolute path to the output directory, which must already exist",
)
parser.add_argument("--k", help="lenght of kmers to find", type=int)
parser.add_argument(
    "--group", help="name of the sample group, cannot include '_'", type=str
)
parser.add_argument(
    "--header_format",
    help="(optional) format of the sequence headers in input fasta files. If the format is 'gisaid' or 'ncbi' accession numbers will be extracted, otherwhise the whole header will be used as sample id.",
    type=str,
)


def _unknown_nuc():
    return "X"


def _find_kmers(seq_list: list, k, out_dir, group):

    encoding = defaultdict(_unknown_nuc)
    encoding["A"] = "00"
    encoding["C"] = "01"
    encoding["G"] = "10"
    encoding["T"] = "11"
    kmer_set = set()
    sample_list = [seq.id for seq in seq_list]

    for idx, seq in enumerate(seq_list):
        seqid = seq.id
        seq = [encoding[nuc] for nuc in seq.seq]
        kmers_dict = {}
        kmers_pos_dict = {}
        for i in range(len(seq) - k):
            kmer = seq[i : i + k]
            if "X" not in kmer:
                kmer = int("".join(kmer), base=2)
                kmer_set.add(kmer)
                if kmer in kmers_dict.keys():
                    kmers_dict[kmer] += 1
                    kmers_pos_dict[kmer].append(i + 1)
                else:
                    kmers_dict[kmer] = 1
                    kmers_pos_dict[kmer] = [i + 1]

        with open(
            f"{out_dir}/{k}mer_indices/{group}_{seqid}_count.pkl", mode="wb"
        ) as out_file:
            pickle.dump(kmers_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            f"{out_dir}/{k}mer_indices/{group}_{seqid}_pos.pkl", mode="wb"
        ) as out_file:
            pickle.dump(kmers_pos_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

    with open(f"{out_dir}/{group}_{k}mer_set.pkl", mode="wb") as out_file:
        pickle.dump(kmer_set, out_file, protocol=pickle.HIGHEST_PROTOCOL)
    with open(f"{out_dir}/{group}_sample_list.pkl", mode="wb") as out_file:
        pickle.dump(sample_list, out_file, protocol=pickle.HIGHEST_PROTOCOL)


def main():
    args = parser.parse_args()
    print(f"Start finding {args.k}-mers in {args.in_file} {datetime.now()}")

    seq_list = list(SeqIO.parse(args.in_file, format="fasta"))
    for seq in seq_list:

        if args.header_format == "gisaid":
            seq.id = seq.id.split("|")[1]
        elif args.header_format == "ncbi":
            seq.id = seq.id.split("|")[0][:-1]

        if "/" in seq.id:
            raise ValueError(
                "Sequence ids cannot contain '/'. If using data from GISAID please use '--header_format gisaid' option."
            )
        seq.seq = seq.seq.upper()

    _find_kmers(seq_list, args.k, args.out_dir, args.group)
    print(f"Done finding {args.k}-mers in {args.in_file} {datetime.now()}")


if __name__ == "__main__":
    main()
