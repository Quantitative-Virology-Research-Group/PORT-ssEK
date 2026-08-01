import os
import pathlib
import pickle
import multiprocessing
import numpy as np
import pandas as pd
from time import process_time
from Bio import SeqIO

import portek
from portek.portek_utils import BasePipeline, canonical_kmer_id


class KmerFinder(BasePipeline):
    """
    KmerFinder:
    """

    def _read_seq(self, group, file, header_format):
        in_path = f"{self.project_dir}/input/{file}"
        seq_list = list(SeqIO.parse(in_path, format="fasta"))
        for seq in seq_list:
            if header_format == "gisaid":
                if len(seq.id.split("|")) > 1:
                    seq.id = seq.id.split("|")[1]
            elif header_format == "ncbi":
                seq.id = seq.id.split("|")[0][:-1]
            if "/" in seq.id:
                raise ValueError(
                    "Sequence ids cannot contain '/'. Specify correct header formats in the config file."
                )
            seq.seq = seq.seq.upper()

        sample_list = [seq.id for seq in seq_list]
        self.seq_lists.append(seq_list)
        samplelist_path = f"{self.project_dir}/input/indices/"
        if os.path.exists(samplelist_path) == False:
            os.makedirs(samplelist_path)
        with open(f"{samplelist_path}/{group}_sample_list.pkl", mode="wb") as out_file:
            pickle.dump(sample_list, out_file, protocol=pickle.HIGHEST_PROTOCOL)

    def __init__(
        self, project_dir: str, mink: int, maxk: int, strand_specific: bool = False
    ) -> None:
        super().check_min_max_k(mink, maxk)
        super().__init__(project_dir)

        # When strand_specific=True each k-mer is kept as-is (no canonicalization),
        # so sense and antisense inputs yield distinct k-mer sets.
        # When strand_specific=False (default) each k-mer is collapsed to its
        # canonical (lexicographically smaller) form — the original PORT-EK behaviour
        # for double-stranded data.
        self.strand_specific = strand_specific
        if strand_specific:
            print(
                "Strand-specific mode enabled: k-mers will NOT be canonicalized. "
                "Sense and antisense k-mer sets will be kept separate.",
                flush=True,
            )
        else:
            print(
                "Canonical mode (default): k-mers are collapsed to their canonical "
                "form.  Use --strand_specific to compare sense vs. antisense strands.",
                flush=True,
            )

        self.seq_lists = []
        for i, group in enumerate(self.sample_groups):
            self._read_seq(group, self.input_files[i], self.header_format[i])

    def _unknown_nuc(self):
        return "X"

    def _find_kmers(
        self,
        k: int,
        group: str,
        seq_list: list,
        strand_specific: bool = False,
        verbose: bool = False,
    ):
        start_time = process_time()
        if verbose == True:
            print(f"Finding all {k}-mers in {group} sequences.", flush=True)
        kmer_set = set()

        group_size = len(seq_list)
        avg_dict = {}
        freq_dict = {}
        kmers_pos_dict = {}
        for seq in seq_list:
            seqid = seq.id
            seq = portek.encode_seq(seq.seq)
            kmers_dict = {}

            for i in range(0, len(seq) - k + 1):
                kmer = seq[i : i + k]
                if "X" not in kmer:
                    kmer = int("".join(kmer), base=2)
                    # --- Strand-awareness ---
                    # In canonical (default) mode collapse the k-mer to the
                    # lexicographically smaller of itself and its RC so that
                    # sense/antisense are treated as equivalent.  In
                    # strand_specific mode keep the raw integer so the two
                    # strands produce genuinely different k-mer sets.
                    if not strand_specific:
                        kmer = canonical_kmer_id(kmer, k)
                    kmer_set.add(kmer)

                    if kmer in kmers_dict.keys():
                        kmers_dict[kmer] += 1
                        kmers_pos_dict[kmer].append(i + 1)
                        avg_dict[kmer] += 1 / group_size
                    else:
                        kmers_dict[kmer] = 1
                        if kmer in freq_dict.keys():
                            freq_dict[kmer] += 1 / group_size
                            avg_dict[kmer] += 1 / group_size
                            kmers_pos_dict[kmer].append(i + 1)
                        else:
                            freq_dict[kmer] = 1 / group_size
                            avg_dict[kmer] = 1 / group_size
                            kmers_pos_dict[kmer] = [i + 1]

            with open(
                f"{self.project_dir}/input/indices/{k}mers/{group}_{seqid}_count.pkl",
                mode="wb",
            ) as out_file:
                pickle.dump(kmers_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/indices/{k}mer_{group}_set.pkl", mode="wb"
        ) as out_file:
            pickle.dump(kmer_set, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/indices/{k}mer_{group}_freq_dict.pkl", mode="wb"
        ) as out_file:
            pickle.dump(freq_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/indices/{k}mer_{group}_avg_dict.pkl", mode="wb"
        ) as out_file:
            pickle.dump(avg_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            f"{self.project_dir}/input/indices/{k}mer_{group}_pos_dict.pkl",
            mode="wb",
        ) as out_file:
            pickle.dump(kmers_pos_dict, out_file, protocol=pickle.HIGHEST_PROTOCOL)

        if verbose == True:
            print(
                f"Done finding all {k}-mers in {group} sequences.",
                flush=True,
            )

        return (k, process_time() - start_time)

    def find_all_kmers(
        self, n_jobs: int = 4, verbose: bool = False
    ) -> tuple[dict, np.floating]:
        print(
            f"Finding all k-mers of lengths {self.mink} to {self.maxk} in {len(self.sample_groups)} files.",
            flush=True,
        )
        k_to_test = [k for k in range(self.mink, self.maxk + 1, 2)]
        k_to_test.reverse()

        find_kmers_pool_input = []
        for k in k_to_test:
            indices_path = f"{self.project_dir}/input/indices/{k}mers/"
            if os.path.exists(indices_path) == False:
                os.makedirs(indices_path)
            groups = zip(self.sample_groups, self.seq_lists)
            for group, seq_list in groups:
                find_kmers_pool_input.append(
                    [k, group, seq_list, self.strand_specific, verbose]
                )

        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            times = pool.starmap(self._find_kmers, find_kmers_pool_input, chunksize=1)

        print("Done finding all k-mers!")

        time_dict = {k: 0.0 for k in k_to_test}
        for time in times:
            time_dict[time[0]] += time[1]

        avg_seq_len = np.mean(
            [len(seq) for seq_list in self.seq_lists for seq in seq_list]
        )
        if verbose == True:
            print(f"Average sequence length: {avg_seq_len}")

        return time_dict, avg_seq_len


class FindOptimalKPipeline(BasePipeline):
    """
    FindOptimalKPipeline:
    """

    def __init__(
        self, project_dir: str, mink: int, maxk: int, times: dict, avg_seq_len: float
    ) -> None:
        super().check_min_max_k(mink, maxk)
        super().__init__(project_dir)

        self.times = times
        self.avg_seq_len = avg_seq_len

    def _calc_metrics(self, k: int, verbose: bool = False):

        start_time = process_time()

        if verbose == True:
            print(f"Calculating metrics for {k}-mers.", flush=True)

        total_kmer_num = self.avg_seq_len - k + 1
        all_kmer_num = 4**k
        kmer_uniqueness_prob = np.exp(-1 * total_kmer_num**2 / (2 * all_kmer_num))
        spec = kmer_uniqueness_prob

        sample_list = super().load_sample_list()
        tot_samples = len(sample_list)

        group_len_dict = {
            group: sum(1 for s in sample_list if s.split("_")[0] == group)
            for group in self.sample_groups
        }

        # H(F) >= H(2/N)  iff  2 <= count <= N-2, so skip entropy entirely
        # and just count kmers in that range using the freq_dicts directly.
        total_counts: dict[int, int] = {}
        in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
            f"{k}mer_*_freq_dict.pkl"
        )
        for filename in in_path:
            group = filename.stem.split("_")[1]
            group_size = group_len_dict[group]
            with open(filename, mode="rb") as in_file:
                freq_dict: dict = pickle.load(in_file)
            for kmer, freq in freq_dict.items():
                total_counts[kmer] = total_counts.get(kmer, 0) + round(
                    freq * group_size
                )

        tot_kmer = len(total_counts)
        sig_kmer = sum(1 for c in total_counts.values() if 2 <= c <= tot_samples - 2)

        if verbose == True:
            min_F = 2 / tot_samples
            min_H = -(min_F * np.log2(min_F) + (1 - min_F) * np.log2(1 - min_F))
            print(
                f"{sig_kmer} of {tot_kmer} {k}-mers left after filtering with entropy threshold of {min_H}.",
                flush=True,
            )

        dt = process_time() - start_time

        if verbose == True:
            print(f"Done calculating metrics for {k}-mers.", flush=True)

        return k, spec, dt, sig_kmer

    def find_optimal_k(self, n_jobs: int = 4, verbose: bool = False):
        print("Finding optimal k.")
        k_to_test = [(k, verbose) for k in range(self.mink, self.maxk + 1, 2)]

        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            results = pool.starmap(self._calc_metrics, k_to_test, chunksize=1)

        result_df = pd.DataFrame(
            0,
            index=range(self.mink, self.maxk + 1, 2),
            columns=["spec", "dt"],
            dtype=float,
        )

        for result in results:
            if result != None:
                if result[3] >= 100:
                    result_df.loc[result[0], "spec"] = result[1]
                    result_df.loc[result[0], "dt"] = result[2]
                else:
                    result_df = result_df.drop(labels=[result[0]])
                    print(
                        f"Fewer than 100 {result[0]}-mers passed the entropy filter, they will be ignored."
                    )

        for k, time in self.times.items():
            if k in result_df.index:
                result_df.loc[k, "dt"] += time
        result_df["dt"] = result_df["dt"].round(3)
        result_df.sort_index(inplace=True)

        distances = self._get_distances_to_line(result_df)
        best_k = result_df.index[np.argmax(distances)]

        if len(result_df) == 0:
            print(
                "\nNo value of k resulted in more than 100 k-mers passing the entropy filter!"
            )
            print("The data set is likely too small, too diverse, or too conserved.")
            return None
        with open(
            f"{self.project_dir}/output/k_selection_results.txt", mode="w"
        ) as out_file:
            out_file.write(f"\nHere are the results of optimal k selection:\n")
            for k in result_df.index:
                out_file.write(
                    f"\nk: {k}, Expected probability of unique random {k}-mers: {round(result_df.loc[k,'spec']*100, 2)}%, CPU time {result_df.loc[k,'dt']} s."
                )
            out_file.write(
                f"\n\nPORT-EK thinks k value of {best_k} is the best for your data!"
            )
        with open(
            f"{self.project_dir}/output/k_selection_results.txt", mode="r"
        ) as out_file:
            if verbose == True:
                print(out_file.read())
            else:
                lines = out_file.readlines()
                tail = lines[:3] + lines[-2:]
                print(*tail)

    def _get_distances_to_line(self, result_df: pd.DataFrame) -> np.ndarray:
        x = np.array(result_df.index, dtype=float)
        y = result_df["spec"].values
        x_n = (x - x[0]) / (x[-1] - x[0])
        y_n = (y - y.min()) / (y.max() - y.min()) if y.max() != y.min() else y

        dx, dy = x_n[-1] - x_n[0], y_n[-1] - y_n[0]
        distances = np.abs(
            dy * x_n
            - dx * y_n
            + x_n[-1] * y_n[0]
            - dx * y_n
            + x_n[-1] * y_n[0]
            - x_n[0] * y_n[-1]
        ) / np.sqrt(dx**2 + dy**2)

        return distances
