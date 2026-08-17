from Bio.Phylo._io import write
from Bio.Phylo import TreeConstruction
from Bio.Phylo.BaseTree import Tree
import numpy as np
import pandas as pd
from scipy.spatial import distance


class KmerPhyloTreeConstructor:

    def __init__(
        self,
        kmer_counts_path: str,
        subsample_size: int | None = None,
        balance_groups: bool = False,
        verbose: bool = False,
    ) -> None:
        if verbose:
            print("Loading k-mer counts data...")
        self.kmer_counts_df = pd.read_csv(kmer_counts_path, index_col=0)
        if verbose:
            print(
                f"Constructing phylogenetic tree from k-mer counts with {subsample_size} samples..."
            )
        if subsample_size is not None:
            fraction = subsample_size / len(self.kmer_counts_df)
            self.kmer_counts_df = self.kmer_counts_df.groupby(
                "sample_group", group_keys=False
            )
            per_group_size = subsample_size // len(self.kmer_counts_df)
            if balance_groups:
                self.kmer_counts_df = self.kmer_counts_df.sample(
                    n=per_group_size, random_state=42
                )
            else:
                self.kmer_counts_df = self.kmer_counts_df.sample(
                    frac=fraction, random_state=42
                )
        self.kmer_counts_df = self.kmer_counts_df.drop(columns=["sample_group"])
        self.distance_matrix: list | None = None
        self.tree: Tree | None = None

    def format_distance_matrix_for_biopyton(self, verbose: bool = False) -> None:
        if verbose:
            print("Calculating distance matrix...")
        square_matrix = self._calculate_distance_matrix()
        if verbose:
            print("Formatting distance matrix for Biopython...")
        n = square_matrix.shape[0]
        out_matrix = [square_matrix[i, : i + 1].tolist() for i in range(n)]
        self.distance_matrix = out_matrix

    def _calculate_distance_matrix(self) -> np.ndarray:
        distances = distance.squareform(
            distance.pdist(self.kmer_counts_df.values, metric="euclidean")
        )
        return distances

    def construct_tree(self, method: str, verbose: bool = False) -> Tree | None:
        kmer_dm = TreeConstruction.DistanceMatrix(
            self.kmer_counts_df.index.to_list(), self.distance_matrix
        )
        distance_tree_constructor = TreeConstruction.DistanceTreeConstructor()
        if verbose:
            print(f"Constructing phylogenetic tree using {method} method...")
        if method == "nj":
            kmer_tree = distance_tree_constructor.nj(kmer_dm)
        elif method == "upgma":
            kmer_tree = distance_tree_constructor.upgma(kmer_dm)
        else:
            raise ValueError("Method must be either 'nj' or 'upgma'")

        self.tree = kmer_tree

    def write_tree(self, output_path: str, format: str, verbose: bool = False) -> None:
        if self.tree is None:
            raise ValueError("Tree has not been constructed yet.")
        write(self.tree, output_path, format=format)
        if verbose:
            print(f"Phylogenetic tree saved to {output_path}.")
