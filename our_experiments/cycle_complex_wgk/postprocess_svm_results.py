import os
import re
import argparse
from typing import Dict, List
import numpy as np
import pandas as pd

# Directory containing the raw CSV outputs
TOP_OUTPUT_DIR = os.path.join('.', 'our_experiments', 'cycle_complex_wgk', 'topowgk_outputs')


dname_baseline = {
    'BZR_MD': [0.619355, 0.109972, 0.618951, 0.167287, 0.670000, 0.108232],
    'MSRC_21': [0.865487, 0.043243, 0.862474, 0.041771, 0.987434, 0.032442],
    'MSRC_9': [0.904444, 0.067077, 0.902731, 0.067633, 0.985556, 0.028982],
    'ENZYMES': [0.563333, 0.065505, 0.557807, 0.067417, 0.834833, 0.046983],
    'COX2_MD': [0.575410, 0.091193, 0.573452, 0.119436, 0.600591, 0.100467],
    'DHFR_MD': [0.640506, 0.145119, 0.633295, 0.179665, 0.626222, 0.122307],
    'ER_MD': [0.701111, 0.162795, 0.701813, 0.217219, 0.758720, 0.164358],
    'IMDB-BINARY': [0.717000, 0.121821, 0.716502, 0.121155, 0.784195, 0.194736],
    'IMDB-MULTI': [0.402000, 0.054564, 0.396991, 0.070308, 0.593108, 0.051730],
    'COX2': [0.792553, 0.065194, 0.785728, 0.085700, 0.794260, 0.165638],
    'BZR': [0.841975, 0.164383, 0.841396, 0.155878, 0.870588, 0.118176],
    'DHFR': [0.806579, 0.249315, 0.804528, 0.248347, 0.889338, 0.253118],
    'PROTEINS': [0.772646, 0.032888, 0.770225, 0.034706, 0.835079, 0.032712],
    'Mutagenicity': [0.823733, 0.007661, 0.823752, 0.007742, 0.896813, 0.008614],
    'FIRSTMM_DB': [0.761538, 0.076495, 0.721538, 0.090519, 0.810227, 0.031787],
    'NCI1': [0.850608, 0.010618, 0.850312, 0.010691, 0.910675, 0.011267],
    'NCI109': [0.837772, 0.010161, 0.837347, 0.010254, 0.900845, 0.008164]
}
    

def _parse_filename(fn: str):
    """Parse filenames like: {dataset}_svm_result_niter_tn_hcc_(n_m).csv

    Returns (dataset, niter_tn, niter_hcc) or None if not matching.
    """
    m = re.match(r'^(?P<dataset>.+)_svm_result_niter_tn_hcc_\((?P<n>\d+)_(?P<m>\d+)\)\.csv$', fn)
    if not m:
        return None
    return m.group('dataset'), int(m.group('n')), int(m.group('m'))

def aggregate_topowgk_outputs(input_dir: str = TOP_OUTPUT_DIR) -> List[str]:
    """
    Aggregate per-dataset CSVs into total_{dataset}_svm_results.csv files.

    Scans `input_dir` for files matching the pattern, groups by dataset name,
    inserts two columns `niter_tn` and `niter_hcc` between `gk_gamma` and `acc`,
    and writes a combined CSV `total_{dataset}_svm_results.csv` into `input_dir`.

    Returns list of generated file paths.
    """
    files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    groups: Dict[str, List[str]] = {}
    for f in files:
        parsed = _parse_filename(f)
        if parsed is None:
            continue
        dataset, n, m = parsed
        groups.setdefault(dataset, []).append(f)

    generated = []
    for dataset, flist in groups.items():
        dfs = []
        for f in sorted(flist):
            path = os.path.join(input_dir, f)
            try:
                df = pd.read_csv(path)
            except Exception as e:
                print(f"Warning: failed to read {path}: {e}")
                continue

            parsed = _parse_filename(f)
            if parsed is None:
                continue
            _, niter_tn, niter_hcc = parsed

            # Default: append if cannot find gk_gamma/acc
            if 'gk_gamma' in df.columns and 'acc' in df.columns:
                gk_idx = list(df.columns).index('gk_gamma')
                # insert niter_tn right after gk_gamma
                df.insert(gk_idx + 1, 'niter_tn', niter_tn)
                # insert niter_hcc after niter_tn
                df.insert(gk_idx + 2, 'niter_hcc', niter_hcc)
            elif 'gk_gamma' in df.columns:
                gk_idx = list(df.columns).index('gk_gamma')
                df.insert(gk_idx + 1, 'niter_tn', niter_tn)
                df.insert(gk_idx + 2, 'niter_hcc', niter_hcc)
            elif 'acc' in df.columns:
                acc_idx = list(df.columns).index('acc')
                df.insert(acc_idx, 'niter_tn', niter_tn)
                df.insert(acc_idx + 1, 'niter_hcc', niter_hcc)
            else:
                # neither column present: append to end
                df['niter_tn'] = niter_tn
                df['niter_hcc'] = niter_hcc

            dfs.append(df)

        if not dfs:
            continue

        # Concatenate, align columns
        combined = pd.concat(dfs, ignore_index=True, sort=False)

        # Ensure column order places niter_tn and niter_hcc between gk_gamma and acc if possible
        cols = list(combined.columns)
        if 'gk_gamma' in cols and 'acc' in cols:
            gk_idx = cols.index('gk_gamma')
            # remove niter cols if they exist to reinsert in correct position
            for c in ['niter_tn', 'niter_hcc']:
                if c in cols:
                    cols.remove(c)
            cols[gk_idx+1:gk_idx+1] = ['niter_tn', 'niter_hcc']
            combined = combined[cols]

        out_name = f"total_{dataset}_svm_results.csv"
        out_path = os.path.join(input_dir, out_name)
        combined.to_csv(out_path, index=False)
        generated.append(out_path)
        print(f"Wrote {out_path} ({len(combined)} rows)")

    return generated

def aggregate_svm_results_cli():
    parser = argparse.ArgumentParser(description='Aggregate topowgk SVM CSV outputs per dataset.')
    parser.add_argument('--input-dir', '-i', default=TOP_OUTPUT_DIR, help='Directory with raw CSV outputs')
    parser.add_argument('--dataset', '-d', default=None, help='If set, only aggregate this dataset name')
    args = parser.parse_args()

    if args.dataset:
        # Filter to only files for this dataset by temporarily renaming groups
        all_generated = []
        # build just for specified dataset by scanning filenames
        files = [f for f in os.listdir(args.input_dir) if f.endswith('.csv')]
        matched = [f for f in files if f.startswith(f"{args.dataset}_svm_result_niter_tn_hcc_")]
        if not matched:
            print(f"No files found for dataset: {args.dataset}")
            return
        # write matched files to a temp subdir? Simpler: process by reading matched list
        dfs = []
        for f in sorted(matched):
            path = os.path.join(args.input_dir, f)
            try:
                df = pd.read_csv(path)
            except Exception as e:
                print(f"Warning: failed to read {path}: {e}")
                continue
            parsed = _parse_filename(f)
            if parsed is None:
                continue
            _, niter_tn, niter_hcc = parsed
            if 'gk_gamma' in df.columns and 'acc' in df.columns:
                gk_idx = list(df.columns).index('gk_gamma')
                df.insert(gk_idx + 1, 'niter_tn', niter_tn)
                df.insert(gk_idx + 2, 'niter_hcc', niter_hcc)
            elif 'acc' in df.columns:
                acc_idx = list(df.columns).index('acc')
                df.insert(acc_idx, 'niter_tn', niter_tn)
                df.insert(acc_idx + 1, 'niter_hcc', niter_hcc)
            else:
                df['niter_tn'] = niter_tn
                df['niter_hcc'] = niter_hcc
            dfs.append(df)

        if not dfs:
            print(f"No readable files for dataset: {args.dataset}")
            return
        combined = pd.concat(dfs, ignore_index=True, sort=False)
        cols = list(combined.columns)
        if 'gk_gamma' in cols and 'acc' in cols:
            gk_idx = cols.index('gk_gamma')
            for c in ['niter_tn', 'niter_hcc']:
                if c in cols:
                    cols.remove(c)
            cols[gk_idx+1:gk_idx+1] = ['niter_tn', 'niter_hcc']
            combined = combined[cols]
        out_path = os.path.join(args.input_dir, f"total_{args.dataset}_svm_results.csv")
        combined.to_csv(out_path, index=False)
        print(f"Wrote {out_path} ({len(combined)} rows)")
        return

    # default: aggregate all datasets found
    aggregate_topowgk_outputs(args.input_dir)


def extract_best_from_total_files(input_dir: str = TOP_OUTPUT_DIR,
                                  output_dir: str = None) -> List[str]:
    """For each total_{dataset}_svm_results.csv in `input_dir`, select rows
    with maximum `acc`, `f1_score`, and `roc_auc_score` and write them to
    `best_{dataset}_svm_result.csv` under `output_dir`.

    If a single row is the best for multiple metrics it will be written once
    and the `selected_metric` column will list the metrics (semicolon-separated).

    Returns list of written file paths.
    """
    if output_dir is None:
        output_dir = os.path.join(TOP_OUTPUT_DIR, 'processed_svm_results')
    os.makedirs(output_dir, exist_ok=True)

    files = [f for f in os.listdir(input_dir) if f.startswith('total_') and f.endswith('.csv')]
    written = []
    for f in sorted(files):
        path = os.path.join(input_dir, f)
        try:
            df = pd.read_csv(path)
        except Exception as e:
            print(f"Warning: failed to read {path}: {e}")
            continue

        # Determine which metric columns are present and numeric
        metric_cols = [c for c in ['acc', 'f1_score', 'roc_auc_score'] if c in df.columns and pd.api.types.is_numeric_dtype(df[c])]
        if not metric_cols:
            print(f"No numeric metric columns (acc/f1_score/roc_auc_score) found in {f}; skipping")
            continue

        # Compute maxima for present metrics
        maxima = {c: df[c].max(skipna=True) for c in metric_cols}

        # Selection priority: all metrics (if all present) -> any pair -> any single
        selected_mask = None

        # Try all metrics (only possible if all three present)
        if set(['acc', 'f1_score', 'roc_auc_score']).issubset(set(metric_cols)):
            mask = np.ones(len(df), dtype=bool)
            for c in ['acc', 'f1_score', 'roc_auc_score']:
                mask &= np.isclose(df[c].to_numpy(), maxima[c], equal_nan=False)
            if mask.any():
                selected_mask = mask

        # If no triple, try pairs
        if selected_mask is None:
            from itertools import combinations
            pair_masks = []
            for a, b in combinations(metric_cols, 2):
                mask = np.isclose(df[a].to_numpy(), maxima[a], equal_nan=False) & np.isclose(df[b].to_numpy(), maxima[b], equal_nan=False)
                if mask.any():
                    pair_masks.append(mask)
            if pair_masks:
                # combine all pair masks (rows satisfying any pair)
                combined_pair = np.zeros(len(df), dtype=bool)
                for pm in pair_masks:
                    combined_pair |= pm
                selected_mask = combined_pair

        # If still none, pick any rows that maximize at least one metric
        if selected_mask is None:
            single_masks = [np.isclose(df[c].to_numpy(), maxima[c], equal_nan=False) for c in metric_cols]
            combined_single = np.zeros(len(df), dtype=bool)
            for sm in single_masks:
                combined_single |= sm
            if combined_single.any():
                selected_mask = combined_single

        if selected_mask is None or not selected_mask.any():
            print(f"No rows found for selection in {f}; skipping")
            continue

        selected = df.loc[selected_mask].copy()

        # For each selected row, compute which metrics reach their maxima
        def _row_selected_metrics(row):
            mets = []
            for c in metric_cols:
                try:
                    if np.isclose(row[c], maxima[c], equal_nan=False):
                        mets.append(c)
                except Exception:
                    continue
            return ';'.join(sorted(mets))

        selected['selected_metric'] = selected.apply(_row_selected_metrics, axis=1)

        dataset = f[len('total_'):-len('_svm_results.csv')] if f.endswith('_svm_results.csv') else f
        out_name = f"best_{dataset}_svm_result.csv"
        out_path = os.path.join(output_dir, out_name)
        selected.to_csv(out_path, index=False)
        written.append(out_path)
        print(f"Wrote best file {out_path} ({len(selected)} rows)")

    return written

def select_best_priority_and_generate_profiles(processed_dir: str = None, top_output_dir: str = TOP_OUTPUT_DIR) -> str:
    """Select best entries per dataset from `best_{dataset}_svm_result.csv` files
    in `processed_dir` according to priority (acc, then roc_auc_score, then f1_score),
    write a combined `best_all_dataset_svm_results.csv`, and generate three
    profile CSVs per selected entry (vary alpha, vary gk_gamma, vary niter_hcc).

    Returns the path to the combined CSV.
    """
    # locate processed_dir if not provided
    # processed_svm_results is located under TOP_OUTPUT_DIR
    proc_dir = processed_dir or os.path.join(TOP_OUTPUT_DIR, 'processed_svm_results')
    if not os.path.isdir(proc_dir):
        raise FileNotFoundError(f'processed_svm_results directory not found at {proc_dir}')

    best_files = [f for f in os.listdir(proc_dir) if f.startswith('best_') and f.endswith('.csv')]
    selected_rows = []
    for bf in sorted(best_files):
        path = os.path.join(proc_dir, bf)
        try:
            df = pd.read_csv(path)
        except Exception as e:
            print(f"Warning: cannot read {path}: {e}")
            continue

        # ensure numeric columns handled
        if 'acc' not in df.columns:
            print(f"Skipping {bf}: no 'acc' column")
            continue

        # lexicographic selection: acc -> roc_auc_score -> f1_score
        # Step 1: max acc
        max_acc = df['acc'].max(skipna=True)
        mask = np.isclose(df['acc'].to_numpy(), max_acc, equal_nan=False)
        df_masked = df[mask]

        # Step 2: among those, max roc_auc_score if present
        if 'roc_auc_score' in df_masked.columns and pd.api.types.is_numeric_dtype(df_masked['roc_auc_score']):
            max_roc = df_masked['roc_auc_score'].max(skipna=True)
            mask2 = np.isclose(df_masked['roc_auc_score'].to_numpy(), max_roc, equal_nan=False)
            df_masked = df_masked[mask2]

        # Step 3: among those, max f1_score if present
        if 'f1_score' in df_masked.columns and pd.api.types.is_numeric_dtype(df_masked['f1_score']):
            max_f1 = df_masked['f1_score'].max(skipna=True)
            mask3 = np.isclose(df_masked['f1_score'].to_numpy(), max_f1, equal_nan=False)
            df_masked = df_masked[mask3]

        # add dataset column from filename 'best_{dataset}_svm_result.csv'
        dataset = bf[len('best_'):-len('_svm_result.csv')] if bf.endswith('_svm_result.csv') else bf
        df_masked = df_masked.copy()
        df_masked.insert(0, 'dataset', dataset)
        selected_rows.append(df_masked)

    if not selected_rows:
        raise RuntimeError('No selected rows found in processed best files')

    combined = pd.concat(selected_rows, ignore_index=True, sort=False)

    # select and order required columns for best_all file
    out_cols = ['dataset', 'alpha', 'gk_gamma', 'niter_tn', 'niter_hcc',
                'acc', 'acc_std', 'f1_score', 'f1_score_std', 'roc_auc_score', 'roc_auc_score_std']
    for c in out_cols:
        if c not in combined.columns:
            combined[c] = pd.NA
    best_all = combined[out_cols]

    out_path = os.path.join(proc_dir, 'best_all_dataset_svm_results.csv')
    best_all.to_csv(out_path, index=False)
    print(f"Wrote combined best file {out_path} ({len(best_all)} rows)")

    # Helper: safe string for filenames
    def _safe(x):
        s = str(x)
        return s.replace('.', 'p').replace('/', '_').replace(' ', '_')

    # Generate profile CSVs per selected row
    for _, row in best_all.iterrows():
        ds = row['dataset']
        alpha = row['alpha']
        gk = row['gk_gamma']
        niter_tn = row['niter_tn']
        niter_hcc = row['niter_hcc']

        # find the source file for (niter_tn, niter_hcc)
        src_fname = f"{ds}_svm_result_niter_tn_hcc_({int(niter_tn)}_{int(niter_hcc)}).csv" if (pd.notna(niter_tn) and pd.notna(niter_hcc)) else None
        if src_fname:
            src_path = os.path.join(top_output_dir, src_fname)
        else:
            src_path = None

        if not src_path or not os.path.isfile(src_path):
            print(f"Source file not found for {ds} tn={niter_tn} hcc={niter_hcc}: {src_path}; skipping profile generation for this entry")
            continue

        try:
            src_df = pd.read_csv(src_path)
        except Exception as e:
            print(f"Failed to read source {src_path}: {e}")
            continue

        metrics = ['acc', 'acc_std', 'f1_score', 'f1_score_std', 'roc_auc_score', 'roc_auc_score_std']

        # (1) fixed gk_gamma, niter_tn, niter_hcc: varying alpha
        if 'gk_gamma' in src_df.columns:
            mask_gk = np.isclose(src_df['gk_gamma'].to_numpy(), float(gk), equal_nan=False) if pd.notna(gk) else src_df['gk_gamma'].isna()
            df_alpha = src_df[mask_gk].copy()
            # ensure niter columns present
            try:
                tn_val = int(niter_tn)
            except Exception:
                tn_val = niter_tn
            try:
                hcc_val = int(niter_hcc)
            except Exception:
                hcc_val = niter_hcc
            df_alpha['niter_tn'] = tn_val
            df_alpha['niter_hcc'] = hcc_val
            cols_alpha = ['alpha', 'gk_gamma', 'niter_tn', 'niter_hcc'] + [c for c in metrics if c in df_alpha.columns]
            if not df_alpha.empty:
                out1 = df_alpha[cols_alpha]
                out1_path = os.path.join(proc_dir, f"profile_alpha_{ds}_svm_results.csv")
                write_header = not os.path.exists(out1_path)
                out1.to_csv(out1_path, index=False, mode='a', header=write_header)
                print(f"Appended profile (vary alpha) {out1_path} ({len(out1)} rows)")

        # (2) fixed alpha, niter_tn, niter_hcc: varying gk_gamma
        if 'alpha' in src_df.columns:
            mask_alpha = np.isclose(src_df['alpha'].to_numpy(), float(alpha), equal_nan=False) if pd.notna(alpha) else src_df['alpha'].isna()
            df_gk = src_df[mask_alpha].copy()
            try:
                tn_val = int(niter_tn)
            except Exception:
                tn_val = niter_tn
            try:
                hcc_val = int(niter_hcc)
            except Exception:
                hcc_val = niter_hcc
            df_gk['niter_tn'] = tn_val
            df_gk['niter_hcc'] = hcc_val
            cols_gk = ['alpha', 'gk_gamma', 'niter_tn', 'niter_hcc'] + [c for c in metrics if c in df_gk.columns]
            if not df_gk.empty:
                out2 = df_gk[cols_gk]
                out2_path = os.path.join(proc_dir, f"profile_gk_gamma_{ds}_svm_results.csv")
                write_header = not os.path.exists(out2_path)
                out2.to_csv(out2_path, index=False, mode='a', header=write_header)
                print(f"Appended profile (vary gk) {out2_path} ({len(out2)} rows)")

        # (3) fixed alpha, gk_gamma, niter_tn: varying niter_hcc across files
        # find all files for this dataset with same niter_tn
        pattern_prefix = f"{ds}_svm_result_niter_tn_hcc_({int(niter_tn)}_"
        matched_files = [ff for ff in os.listdir(top_output_dir) if ff.startswith(pattern_prefix) and ff.endswith('.csv')]
        rows_nh = []
        for mf in sorted(matched_files):
            # extract niter_hcc from filename
            parsed = _parse_filename(mf)
            if parsed is None:
                continue
            _, tn, hcc = parsed
            mfpath = os.path.join(top_output_dir, mf)
            try:
                mdf = pd.read_csv(mfpath)
            except Exception:
                continue
            # find row matching alpha and gk_gamma
            if 'alpha' in mdf.columns and 'gk_gamma' in mdf.columns:
                mask_row = np.isclose(mdf['alpha'].to_numpy(), float(alpha), equal_nan=False) & np.isclose(mdf['gk_gamma'].to_numpy(), float(gk), equal_nan=False)
                sel = mdf.loc[mask_row]
                if not sel.empty:
                    r = sel.iloc[0][[c for c in metrics if c in sel.columns]].to_dict()
                    r['niter_hcc'] = hcc
                    rows_nh.append(r)

        if rows_nh:
            df_nh = pd.DataFrame(rows_nh)
            # ensure alpha, gk_gamma, niter_tn present
            try:
                tn_val = int(niter_tn)
            except Exception:
                tn_val = niter_tn
            df_nh['alpha'] = alpha
            df_nh['gk_gamma'] = gk
            df_nh['niter_tn'] = tn_val
            cols_nh = ['alpha', 'gk_gamma', 'niter_tn', 'niter_hcc'] + [c for c in metrics if c in df_nh.columns]
            out3_path = os.path.join(proc_dir, f"profile_niter_hcc_{ds}_svm_results.csv")
            write_header = not os.path.exists(out3_path)
            df_nh.to_csv(out3_path, index=False, mode='a', header=write_header, columns=cols_nh)
            print(f"Appended profile (vary niter_hcc) {out3_path} ({len(df_nh)} rows)")

    return out_path


if __name__ == '__main__':
    aggregate_svm_results_cli()