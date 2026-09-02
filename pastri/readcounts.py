import dataclasses
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pastri.smahtid import SMaHTid

def get_gt(gt):
    """
    Simple renamer
    """
    gt = tuple(gt)
    if None in gt:
        return "NON"
    # Broken for hemi
    if len(gt) == 1:
        return "HEMI" if gt[0] else "REF"
    if len(gt) != 2:
        return "UNK"
    if gt == (0, 0):
        return "REF"
    if gt[0] != gt[1]:
        return "HET"
    return "HOM"


def load_parquet(pq_fn, filtering=True):
    """
    Load variants from a parquet file
    With filtering, automatically subset to passing SVs with a donor-level genotype of 0/0
    return DataFrame of ref/alt read counts. Column names are samples and index are SV IDs
    """
    schema = pq.read_schema(pq_fn)
    all_columns = schema.names
    sample_gt = [_ for _ in all_columns if _.endswith(
        '_GT') and _.count('-') == 0][0]
    columns = [_.split('_')[0] for _ in all_columns
               if _.endswith('_GT') and _.count('-') > 1]
    needed_columns = (
        ['id', 'is_pass', sample_gt]
        + [f'{_}_AD_ref' for _ in columns]
        + [f'{_}_AD_alt' for _ in columns]
    )

    data = pd.read_parquet(pq_fn, columns=needed_columns)

    if filtering:
        data['GT'] = data[sample_gt].apply(get_gt)
        keep = (data['GT'] == 'REF') & (data['is_pass'])
        data = data[keep]
    data.set_index('id', inplace=True)
    
    ref_reads = data[[f'{_}_AD_ref' for _ in columns]]
    alt_reads = data[[f'{_}_AD_alt' for _ in columns]]
    ref_reads.columns = columns
    alt_reads.columns = columns
    return ref_reads, alt_reads

def load_tsv(ref_fn, alt_fn):
    ref_idx = pd.read_csv(ref_fn, sep='\t', nrows=0).columns[0]
    ref_reads = pd.read_csv(ref_fn, sep='\t', dtype={ref_idx: str})
    print(f"Setting --ref row idx to {ref_idx}")
    ref_reads.set_index(ref_idx, inplace=True)
    
    alt_idx = pd.read_csv(alt_fn, sep='\t', nrows=0).columns[0]
    alt_reads = pd.read_csv(alt_fn, sep='\t', dtype={alt_idx: str})
    print(f"Setting --alt row idx to {alt_idx}")
    alt_reads.set_index(alt_idx, inplace=True)

    assert (ref_reads.columns == alt_reads.columns).all(), "--ref/--alt columns do not match"
    assert (ref_reads.index == alt_reads.index).all(), "--ref/--alt indexes do not match"
    return ref_reads, alt_reads
    



class ReadCounts():
    """
    Class for holding reference and alternate read counts, tracking the samples (columns) and index (rows)
    """

    def __init__(self, ref_reads, alt_reads, columns=None, index=None):
        self.ref = np.array(ref_reads)
        self.alt = np.array(alt_reads)

        n_row, n_col = self.ref.shape
        self.col = np.array(
            columns) if columns is not None else np.arange(n_col)
        self.idx = np.array(index) if index is not None else np.arange(n_row)
        assert self.ref.shape == self.alt.shape, "Incompatable ref/alt shape"
        assert n_col == len(self.col), "Invalid column names"
        assert n_row == len(self.idx), "Invalid index names"
    
    @staticmethod
    def from_pq(pq_fn, filtering=True):
        """
        Build ReadCount from a parquet file with SMaHTid columns
        """
        ref_reads, alt_reads = load_parquet(pq_fn, filtering)
        names = [SMaHTid(_) for _ in ref_reads.columns]
        return ReadCounts(ref_reads, alt_reads, names, ref_reads.index)

    @staticmethod
    def from_tsv(ref_fn, alt_fn):
        """
        Build ReadCount from a tsv with SMaHTid columns
        """
        ref_reads, alt_reads = load_tsv(ref_fn, alt_fn)
        names = [SMaHTid(_) for _ in ref_reads.columns]
        return ReadCounts(ref_reads, alt_reads, names, ref_reads.index)

    @staticmethod
    def __matrix_collapse(names, matrix):
        """
        Collapses a 2d matrix by matching names
        """
        new_names, inverse = np.unique(names, return_inverse=True)
        new_matrix = np.zeros((matrix.shape[0], len(new_names)),
                              dtype=matrix.dtype)
        np.add.at(new_matrix, (slice(None), inverse),
                  np.nan_to_num(matrix, nan=0.0))
        return new_names, new_matrix

    def to_tsv(self, out_prefix, compress=True):
        """
        Save to `.ref.tsv[.gz]` and `.alt.tsv[.gz]` that can be loaded with ReadCounts.from_tsv
        """
        suffix = '.gz' if compress else ''
        pd.DataFrame(self.ref, index=self.idx, columns=self.col).to_csv(f"{out_prefix}.ref.tsv{suffix}", sep='\t')
        pd.DataFrame(self.alt, index=self.idx, columns=self.col).to_csv(f"{out_prefix}.alt.tsv{suffix}", sep='\t')

    def min_mask(self, min_coverage=10, min_col=2, **kwargs):
        """
        Return a copy of self with ref/alt counts converted to NA if cell has < min_coverage or
        if rows have fewer than min_col non-NA
        """
        total = self.ref + self.alt

        # Mask cells below min_coverage
        low_cov = total < min_coverage
        ref_masked = self.ref.astype(float).copy()
        alt_masked = self.alt.astype(float).copy()
        ref_masked[low_cov] = np.nan
        alt_masked[low_cov] = np.nan

        # Mask entire rows where fewer than min_col cells survive
        valid_cols = ~low_cov  # shape (n_row, n_col)
        insufficient_rows = valid_cols.sum(axis=1) < min_col  # shape (n_row,)
        ref_masked[insufficient_rows, :] = np.nan
        alt_masked[insufficient_rows, :] = np.nan

        return ReadCounts(ref_masked, alt_masked, self.col, self.idx)

    def collapse_by(self, names):
        """
        Sum reference/alternate read counts that have the same names
        """
        new_names, ref_reads = self.__matrix_collapse(names, self.ref)
        _, alt_reads = self.__matrix_collapse(names, self.alt)
        return ReadCounts(ref_reads, alt_reads, new_names, self.idx)

    def col_present_sum(self):
        """
        Return a pd.Series of the total present (alt != 0) rows across columns
        """
        return pd.Series(((self.alt != 0) & (~np.isnan(self.alt))).sum(axis=0), name='Total', index=self.col)

    def exclude_finished(self, results):
        """
        Helper method for excluding idx that have finished their testing
        """
        return self[~np.isin(self.idx, results[results['finished']]['idx'].unique())]

    def presence_lookup(self):
        """
        Make a presence lookup Series with index of idx and a presence state' off True present, False absent, or np.nan
        To place the presence annotation into a df, use `df['idx'].map(lookup)`
        """
        m_idx = np.array(self.idx)
        unchecked = np.isnan(self.ref).all(
            axis=1) | np.isnan(self.alt).all(axis=1)
        has_alt = np.nansum(self.alt, axis=1) > 0
        conds = [
            unchecked,
            has_alt,
        ]
        choice = [
            None,
            True
        ]
        state = np.select(conds, choice, False).astype(object)
        state[unchecked] = None
        return pd.Series(state, name='state', index=m_idx)

    def to_frame(self):
        """
        Turn these read counts into a flat pandas DataFrame
        """
        n_rows = self.ref.shape[0]

        row_index = np.repeat(self.idx, len(self.col))
        col_index = np.tile(self.col, n_rows)

        return pd.DataFrame({
            'idx':  row_index,
            'col': col_index,
            'ref':    self.ref.ravel(),
            'alt':    self.alt.ravel(),
        })

    def dbg_view(self, idx=None, full=False, present=False):
        """
        Quick viewer of a single row when the col is a SMaHTid
        """
        if idx:
            df = self[self.idx == idx].to_frame()
        else:
            df = self.to_frame()
        expanded = df['col'].apply(lambda x: pd.Series(dataclasses.asdict(x)))
        df = pd.concat([df, expanded], axis=1)
        if present:
            df = df[df['alt'] != 0]
        if not full:
            df = df[['idx', 'col', 'ref', 'alt', 'protocol',
                     'tissue_abv', 'layer', 'dedup_tissue_abv', 'platform_name']]
        return df

    def __getitem__(self, key):
        """
        Slicing
        """
        if not isinstance(key, tuple):
            row_key, col_key = key, slice(None)
        else:
            row_key, col_key = key

        # Normalise a bare integer to a length-1 slice so all axes stay 2-D
        if isinstance(row_key, (int, np.integer)):
            row_key = slice(row_key, row_key + 1)

        return ReadCounts(
            self.ref[row_key][:, col_key],
            self.alt[row_key][:, col_key],
            self.col[col_key],
            self.idx[row_key]
        )

