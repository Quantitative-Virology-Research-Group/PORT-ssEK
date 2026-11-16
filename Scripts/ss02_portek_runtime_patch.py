"""
Runtime (in-memory) corrections for PORT-EK-v2 -- the repository files are NOT edited.
Usage in any driver script, BEFORE running the pipeline:   import portek_runtime_patch  (after `import portek`)
Fixes applied only to the running Python process:
  1. assign_kmer_group_ava: exact group membership instead of substring match ("A" in "CRF01.AE").
  2. scipy.stats.false_discovery_control: invalid p-values (>1 or NaN) are clipped, so `--fdr` no longer crashes.
  3. FindOptimalKPipeline._get_distances_to_line: standard point-to-chord distance (optional; k can instead be fixed by --min_k = --max_k).
Gap-inflated sequence length (len(seq) counts "-") is avoided WITHOUT code change by removing gaps from the input FASTA files
(scripts/degap_fasta.py); HXB2 mapping is unaffected.
"""
import functools, numpy as np, portek, portek.portek_enriched as _enr, portek.portek_findk as _fk

@functools.lru_cache(maxsize=None)
def _pair_map(p_cols, groups):
    out = {}
    for c in p_cols:
        s = c[: -len("_p-value")]
        for g in groups:
            if s.startswith(g + "-") and s[len(g) + 1:] in groups: out[c] = (g, s[len(g) + 1:]); break
        else: raise ValueError(f"cannot parse group pair from {c}")
    return out

def assign_kmer_group_ava(row, p_cols, avg_cols, freq_cols, err_cols):
    groups = tuple(c[: -len("_avg")] for c in avg_cols)
    max_group = row[avg_cols].idxmax()[: -len("_avg")]
    pm = _pair_map(tuple(p_cols), groups)
    rel_p_cols = [c for c in p_cols if max_group in pm[c]]
    if all(row[freq_cols] > 0.9) and all(abs(row[err_cols]) < 0.1): return "conserved"
    elif all(row[rel_p_cols] < 0.01): return f"{max_group}_enriched"
    return "not_significant"

def _get_distances_to_line(self, df):
    x = np.array(df.index, dtype=float); y = df["spec"].to_numpy(float)
    x_n = (x - x.min()) / (x.max() - x.min()); y_n = (y - y.min()) / (y.max() - y.min())
    dx, dy = x_n[-1] - x_n[0], y_n[-1] - y_n[0]
    return np.abs(dy * x_n - dx * y_n + x_n[-1] * y_n[0] - x_n[0] * y_n[-1]) / np.sqrt(dx ** 2 + dy ** 2)

# ---- (3) --fdr crash: Mann-Whitney can return p slightly > 1 (or NaN) and scipy.false_discovery_control then raises
import scipy.stats as _st
_orig_fdc = _st.false_discovery_control
def _safe_fdc(ps, *args, **kwargs):
    p = np.clip(np.nan_to_num(np.asarray(ps, dtype=float), nan=1.0), 0.0, 1.0)   # NaN -> 1 (conservative), clip to [0, 1]
    return np.clip(_orig_fdc(p, *args, **kwargs), 0.0, 1.0)   # BY-adjusted p-values can exceed 1; cap at 1
_st.false_discovery_control = _safe_fdc
portek.assign_kmer_group_ava = assign_kmer_group_ava      # pipeline calls portek.assign_kmer_group_ava at run time
_orig_dist = _fk.FindOptimalKPipeline._get_distances_to_line
def apply_findk_patch(): _fk.FindOptimalKPipeline._get_distances_to_line = _get_distances_to_line
