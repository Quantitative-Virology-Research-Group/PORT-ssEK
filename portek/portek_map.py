import itertools
import os
import pathlib
import pickle

import pandas as pd

from portek.portek_utils import BasePipeline


class MappingPipeline(BasePipeline):

    def __init__(self, project_dir: str, k: int):
        super().__init__(project_dir, k)
        self.kmer_index = {}
        self.matrices = {}

        try:
            self.matrices["enriched"] = pd.read_csv(
                f"{project_dir}/output/enriched_{self.k}mers_stats.csv", index_col=0
            )
        except:
            raise FileNotFoundError(
                f"No enriched {self.k}-mers table found in {project_dir}output/ ! Please run PORT-EK find_enriched first!"
            )
        self._initialize_mapping_dataframe()
        self._initialize_coverage_dataframe()

    def _initialize_mapping_dataframe(self):
        self.matrices["mapping"] = pd.DataFrame(
            [[[], "", 0, "", ""]],
            columns=[
                "reference_sequence_position",
                "gene",
                "number_of_mismatches",
                "group",
                "exclusivity",
            ],
            index=self.matrices["enriched"].index,
        )

    def _initialize_coverage_dataframe(self):
        self.matrices["coverage"] = pd.DataFrame(
            0,
            columns=[f"{group}_enriched_kmer_coverage" for group in self.sample_groups]
            + ["conserved_kmer_coverage"],
            index=pd.Index(
                range(1, len(self.ref_seq) + 1), name="reference_sequence_position"
            ),
        )

    def run_mapping(self, max_n_mismatch: int, verbose: bool = False) -> None:
        self.index_ref_seq(max_n_mismatch, verbose)

        for kmer in self.matrices["mapping"].index:
            mapping_dict = self.map_kmer_to_index(kmer, max_n_mismatch)
            self.add_kmer_to_mapping_df(max_n_mismatch, kmer, mapping_dict)
            self.update_mapping_df_group(kmer)
            if self.ref_genes:
                self.update_mapping_df_genes(kmer)

        self.fill_coverage_dataframe()
        self.save_mapping_and_coverage(max_n_mismatch)

    def index_ref_seq(self, max_mismatches: int, verbose: bool = False) -> None:
        if self._check_index(max_mismatches) == False:
            self._generate_and_write_index(max_mismatches, verbose)
        else:
            self._load_index(max_mismatches, verbose)

    def _check_index(self, max_mismatches: int) -> bool:
        return pathlib.Path(
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{max_mismatches}.pkl"
        ).exists()

    def _generate_and_write_index(self, max_mismatches, verbose):
        if verbose == True:
            print(
                f"No {self.k}-mer index with maximum distance {max_mismatches} exists for {self.ref_seq_name}, building index."
            )
        kmer_index = self._generate_index(max_mismatches, self.ref_seq)  # type: ignore
        self._write_index(max_mismatches, kmer_index)

        if verbose == True:
            print(
                f"Finished building {self.k}-mer index with maximum distance {max_mismatches} exists for {self.ref_seq_name}"
            )
            for d in range(max_mismatches + 1):
                print(f"Extracted {len(kmer_index[d])} k-mers with distance {d}.")
        self.kmer_index = kmer_index

    def _generate_index(self, max_mismatches: int, ref_seq: str) -> dict:
        kmer_index = {dist: {} for dist in range(max_mismatches + 1)}
        for kmer_start_position in range(0, len(ref_seq) - self.k + 1):
            kmer = ref_seq[kmer_start_position : kmer_start_position + self.k]
            self._save_kmer_to_index(kmer_index, kmer_start_position, 0, kmer)

            for n_mismatches in range(1, max_mismatches + 1):
                sub_kmers = self._generate_sub_kmers(n_mismatches, kmer)
                indel_kmers = self._generate_indel_kmers(
                    n_mismatches, kmer, kmer_start_position, ref_seq
                )
                ambi_kmers = sub_kmers.union(indel_kmers)
                for ambi_kmer in ambi_kmers:
                    self._save_kmer_to_index(
                        kmer_index, kmer_start_position, n_mismatches, ambi_kmer
                    )

        return kmer_index

    def _generate_sub_kmers(
        self,
        n_mismatches: int,
        kmer: str,
    ) -> set:

        sub_kmers = set()
        position_combinations = list(
            itertools.combinations(range(self.k), n_mismatches)
        )
        for positions in position_combinations:
            ambi_kmer = [nuc for nuc in kmer]
            for pos in positions:
                ambi_kmer[pos] = "N"
            sub_kmers.add("".join(ambi_kmer))

        return sub_kmers

    def _generate_indel_kmers(
        self,
        n_mismatches: int,
        kmer: str,
        kmer_start_position: int,
        ref_seq: str,
        generate_deletions: bool = True,
    ) -> set:

        indel_kmers = set()
        for position in range(1, self.k - n_mismatches):
            in_kmer = [nuc for nuc in kmer]
            in_kmer.insert(position, "N" * n_mismatches)
            indel_kmers.add("".join(in_kmer[:-n_mismatches]))
            if generate_deletions == True:
                if position + n_mismatches >= self.k:
                    continue
                if kmer_start_position + self.k + n_mismatches > len(ref_seq):
                    continue
                del_kmer = [nuc for nuc in kmer]
                del_kmer = (
                    del_kmer[:position]
                    + del_kmer[position + n_mismatches :]
                    + [
                        nuc
                        for nuc in ref_seq[
                            kmer_start_position
                            + self.k : kmer_start_position
                            + self.k
                            + n_mismatches
                        ]
                    ]
                )
                indel_kmers.add("".join(del_kmer))

        return indel_kmers

    def _save_kmer_to_index(self, kmer_index, kmer_start_position, n_mismatches, kmer):
        if kmer in kmer_index[n_mismatches].keys():
            kmer_index[n_mismatches][kmer].append(kmer_start_position + 1)
        else:
            kmer_index[n_mismatches][kmer] = [kmer_start_position + 1]

    def _write_index(self, max_mismatches, kmer_index):
        os.makedirs(f"{self.project_dir}/temp/ref_index", exist_ok=True)
        with open(
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{max_mismatches}.pkl",
            mode="wb",
        ) as out_file:
            pickle.dump(kmer_index, out_file)

    def _load_index(self, max_mismatches, verbose):
        if verbose == True:
            print(
                f"{self.k}-mer index with maximum distance {max_mismatches} for {self.ref_seq_name} already exists, loading..."
            )
        with open(
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{max_mismatches}.pkl",
            mode="rb",
        ) as in_file:
            kmer_index = pickle.load(in_file)
        self.kmer_index = kmer_index

    def map_kmer_to_index(
        self,
        kmer: str,
        max_n_mismatch: int,
    ) -> dict:
        mapping_dict = {n: set() for n in range(max_n_mismatch + 1)}
        for n_mismatch in range(max_n_mismatch + 1):
            sub_kmers = self._generate_sub_kmers(n_mismatch, kmer)
            in_kmers = self._generate_indel_kmers(
                n_mismatch, kmer, 0, "", generate_deletions=False
            )
            ambi_kmers = sub_kmers.union(in_kmers)
            for ambi_kmer in ambi_kmers:
                if ambi_kmer in self.kmer_index[n_mismatch].keys():
                    mapping_dict[n_mismatch] = mapping_dict[n_mismatch].union(
                        self.kmer_index[n_mismatch][ambi_kmer]
                    )
        return mapping_dict

    def add_kmer_to_mapping_df(
        self, max_n_mismatch: int, kmer: str, mapping_dict: dict
    ) -> None:

        for n_mismatch in range(max_n_mismatch + 1):
            if len(mapping_dict[n_mismatch]) > 0:
                self.matrices["mapping"].at[kmer, "reference_sequence_position"] = (
                    sorted(list(mapping_dict[n_mismatch]))
                )
                self.matrices["mapping"].at[kmer, "number_of_mismatches"] = n_mismatch
                break

    def update_mapping_df_group(self, kmer: str) -> None:
        self.matrices["mapping"].loc[kmer, ["group", "exclusivity"]] = self.matrices[
            "enriched"
        ].loc[kmer, ["group", "exclusivity"]]

    def update_mapping_df_genes(self, kmer: str) -> None:
        positions: list[int] = self.matrices["mapping"].at[kmer, "reference_sequence_position"]  # type: ignore
        genes = set()
        if positions:
            for pos in positions:
                genes_at_pos = self._find_genes_for_position(pos)
                if genes_at_pos:
                    genes.update(genes_at_pos.split(","))

        self.matrices["mapping"].at[kmer, "gene"] = ", ".join(sorted(genes))

    def _find_genes_for_position(self, position: int) -> str:
        matching_genes = set()

        for gene in self.ref_genes:
            if isinstance(gene["start"], list):
                for start, end in zip(gene["start"], gene["end"]):
                    if start <= position <= end:
                        matching_genes.add(gene["gene"])
            elif gene["start"] <= position <= gene["end"]:
                matching_genes.add(gene["gene"])

        return ",".join(sorted(matching_genes)) if matching_genes else ""

    def fill_coverage_dataframe(self) -> None:
        for kmer in self.matrices["enriched"].index:
            start_positions: list[int] = self.matrices["mapping"].at[
                kmer, "reference_sequence_position"
            ]
            group: str = self.matrices["enriched"].at[kmer, "group"]
            if start_positions:
                for start_position in start_positions:
                    positions = [
                        pos for pos in range(start_position, start_position + self.k)
                    ]
                    self.matrices["coverage"].loc[
                        positions, f"{group}_kmer_coverage"
                    ] += 1

    def save_mapping_and_coverage(self, max_n_mismatch: int) -> None:
        self.matrices["mapping"]["reference_sequence_position"] = self.matrices[
            "mapping"
        ]["reference_sequence_position"].apply(lambda x: ", ".join(map(str, x)))
        self.matrices["mapping"].to_csv(
            f"{self.project_dir}/output/mapping_{self.k}mers_max_{max_n_mismatch}_mismatches.tsv",
            sep="\t",
        )
        self.matrices["coverage"].to_csv(
            f"{self.project_dir}/output/coverage_{self.k}mers_max_{max_n_mismatch}_mismatches.tsv",
            sep="\t",
            index_label="reference_sequence_position",
        )
