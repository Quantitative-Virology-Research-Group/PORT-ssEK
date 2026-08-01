import argparse
import os
import shutil
from datetime import datetime

import portek

parser = argparse.ArgumentParser(
    description="Main PORT-EK v2 program. Use it to run all PORT-EK tools with the appropriate 'tool' argument."
)
parser.add_argument(
    "tool",
    help=(
        "Name of the PORT-EK function you want to execute. "
        "Choose one of: new, generate_antisense, find_k, find_enriched, map, tree"
    ),
)

parser.add_argument(
    "project_dir",
    help="Path to the project directory. Must not exist for PORTEKrun.py new, must exist for all other tools",
    type=str,
)

parser.add_argument(
    "--input_fasta",
    help=(
        "Path to the sense FASTA file used by 'generate_antisense'. "
        "The reverse complement of every sequence will be written to --output_fasta."
    ),
    type=str,
    default=None,
)

parser.add_argument(
    "--output_fasta",
    help=(
        "Destination path for the antisense FASTA file produced by 'generate_antisense'. "
        "Defaults to <input_fasta_stem>_antisense.fasta in the same directory."
    ),
    type=str,
    default=None,
)

parser.add_argument(
    "--min_k",
    help="Minimum k value to test with PORT-EK find_k. PORT-EK will test all odd k values from min_k up to and including max_k. Default 5.",
    type=int,
    default=5,
)

parser.add_argument(
    "--max_k",
    help="Maximum k value to test with PORT-EK find_k. PORT-EK will test all odd k values from min_k up to and including max_k. Default 31.",
    type=int,
    default=31,
)

parser.add_argument(
    "-k", help="k value for PORT-EK enriched, map and classify", type=int
)

parser.add_argument(
    "--max_mem",
    "-m",
    help="Maximum memory (in GB) size of the count matrix when calculating k-mer statistics in PORT-EK find_enriched. Rare k-mers will be removed if the matrix exceeds this size. Default 2.0",
    type=float,
    default=2.0,
)

parser.add_argument(
    "--min_freq",
    help="Minimum frequency of k-mers to keep when calculating k-mer statistics in PORT-EK find_enriched. K-mers with frequency below this threshold will be removed ONLY if k-mer matrix exceeds the maximum memory size. Default 0.1 (10%).",
    type=float,
    default=0.1,
)

parser.add_argument(
    "--fdr",
    help="Use false discovery rate to control for when calculating k-mer statistics in PORT-EK find_enriched.",
    action="store_true",
)

parser.add_argument(
    "-d",
    help="Maximum edit distance when mapping k-mers to reference sequence. Default 2.",
    type=int,
    default=2,
)

parser.add_argument(
    "--tree_subsample_size",
    "-n",
    help="Number of samples to subsample when constructing the phylogenetic tree. Default None (no subsampling).",
    type=int,
    default=None,
)

parser.add_argument(
    "--balance_groups",
    "-b",
    help="Whether to balance groups when subsampling samples for the phylogenetic tree. Default False.",
    action="store_true",
)

parser.add_argument(
    "--strand_specific",
    "-s",
    help=(
        "Disable canonical k-mer mode so that each k-mer is counted on its "
        "original strand only.  Use this flag when comparing sense RNA sequences "
        "against their reverse-complemented (antisense) counterparts.  Without "
        "this flag PORT-EK collapses every k-mer to its canonical (lexicographically "
        "smaller) form, making sense and antisense inputs yield identical k-mer sets. "
        "Only affects the find_k step; all downstream steps use the stored indices."
    ),
    default=False,
    action="store_true",
)

parser.add_argument(
    "--verbose",
    "-v",
    help="Recieve additional information from some PORT-EK tools. Default False.",
    default=False,
    action="store_true",
)

parser.add_argument(
    "--n_jobs",
    help="Number of processes used in PORT-EK find and PORT-EK enriched. Default 4.",
    default=4,
    type=int,
)


def _new_project(project_dir: str) -> None:
    if os.path.isdir(project_dir) == True:
        raise FileExistsError(
            "This project directory already exists! PORT-EK does not allow overwriting projects. If you REALLY want to overwrite, remove the existing directory manually!"
        )
    os.makedirs(project_dir)
    os.makedirs(f"{project_dir}/input")
    os.makedirs(f"{project_dir}/output")
    os.makedirs(f"{project_dir}/temp")
    shutil.copy2("./templates/config.yaml", project_dir)
    print(
        f"New empty PORT-EK project created in {project_dir}. Please edit config.yaml as required and copy input fasta files into {project_dir}/input."
    )


def main():
    args = parser.parse_args()
    if type(args.project_dir) != str:
        raise ValueError("Please provide a valid project directory name.")

    if args.tool == "new":
        _new_project(args.project_dir)

    elif args.tool == "generate_antisense":
        # ------------------------------------------------------------------
        # Generate a reverse-complement (antisense) FASTA from a sense FASTA.
        # Also verifies that the output sequences are genuinely different from
        # the input so that the subsequent strand_specific run is meaningful.
        # ------------------------------------------------------------------
        if args.input_fasta is None:
            raise ValueError(
                "Please provide an input FASTA with --input_fasta when using "
                "'generate_antisense'."
            )
        input_path = args.input_fasta
        if args.output_fasta is not None:
            output_path = args.output_fasta
        else:
            stem = os.path.splitext(input_path)[0]
            output_path = f"{stem}_antisense.fasta"

        n = portek.reverse_complement_fasta(input_path, output_path)
        print(
            f"generate_antisense complete: {n} sequences written to {output_path}.\n"
            f"Next steps:\n"
            f"  1. Copy {output_path} into your project's input/ directory.\n"
            f"  2. Add it as a separate sample group in config.yaml.\n"
            f"  3. Run:  portek find_k <project_dir> --strand_specific\n"
            f"  4. Then: portek find_enriched <project_dir> -k <k>"
        )

    elif args.tool == "find_k":
        start_time = datetime.now()
        kmer_finder = portek.KmerFinder(
            args.project_dir, args.min_k, args.max_k,
            strand_specific=args.strand_specific,
        )
        times, avg_seq_len = kmer_finder.find_all_kmers(
            n_jobs=args.n_jobs, verbose=args.verbose
        )
        optimal_k_finder = portek.FindOptimalKPipeline(
            args.project_dir, args.min_k, args.max_k, times, avg_seq_len
        )
        optimal_k_finder.find_optimal_k(n_jobs=args.n_jobs, verbose=args.verbose)
        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = end_timeS_ARE_NOT_CANON - start_time
        print(f"\nTotal running time: {running_time}")

    elif args.tool == "find_enriched":
        start_time = datetime.now()
        enriched_kmers_finder = portek.EnrichedKmersPipeline(args.project_dir, args.k)
        enriched_kmers_finder.get_basic_kmer_stats(
            max_mem=args.max_mem, min_freq=args.min_freq
        )
        enriched_kmers_finder.calc_kmer_stats(
            "common",
            n_jobs=args.n_jobs,
            verbose=args.verbose,
            false_discovery_control=args.fdr,
        )
        enriched_kmers_finder.plot_volcanos("common")
        enriched_kmers_found = enriched_kmers_finder.get_enriched_kmers()
        if enriched_kmers_found == True:
            enriched_kmers_finder.save_counts_for_classifier()
            enriched_kmers_finder.save_matrix("enriched")
            enriched_kmers_finder.plot_PCA()
        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = end_timeS_ARE_NOT_CANON - start_time
        print(f"\nTotal running time: {running_time}")

    elif args.tool == "map":
        start_time = datetime.now()
        mapping_pipeline = portek.MappingPipeline(args.project_dir, args.k)
        mapping_pipeline.run_mapping(args.d, args.verbose)
        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = end_timeS_ARE_NOT_CANON - start_time
        print(f"\nTotal running time: {running_time}")

    elif args.tool == "tree":
        start_time = datetime.now()
        kmer_counts_path = (
            f"{args.project_dir}/output/{args.k}mer_counts_for_classifier.csv"
        )
        tree_constructor = portek.KmerPhyloTreeConstructor(
            kmer_counts_path,
            subsample_size=args.tree_subsample_size,
            balance_groups=args.balance_groups,
            verbose=args.verbose,
        )
        tree_constructor.format_distance_matrix_for_biopyton(verbose=args.verbose)
        tree_constructor.construct_tree(method="nj", verbose=args.verbose)
        tree_constructor.write_tree(
            output_path=f"{args.project_dir}/output/{args.k}mer_phylo_tree.nwk",
            format="newick",
            verbose=args.verbose,
        )
        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = end_timeS_ARE_NOT_CANON - start_time
        print(f"\nTotal running time: {running_time}")

    else:
        raise ValueError(
            "Unrecoginzed PORT-EK tool requested. Choose one of: new, find_k, find_enriched, map, tree."
        )


if __name__ == "__main__":
    main()
