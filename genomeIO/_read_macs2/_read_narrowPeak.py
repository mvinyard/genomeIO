
import glob
import os
import pandas as pd

def _return_narrowPeak_columns():

    """"""

    return [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "peak",
    ]


def _read_narrowPeak(path):

    """"""

    narrowPeak_cols = _return_narrowPeak_columns()
    narrowPeaks_df = pd.read_csv(path, sep="\t", header=None, names=narrowPeak_cols)

    return narrowPeaks_df


class narrowPeaks:
    def __init__(self, path):

        """"""
        
        self.path = path
        self.DataFrames = {}

    def read_group(
        self,
        columns=["chrom", "chromStart", "chromEnd"],
        colnames=["chr", "start", "end"],
    ):

        """"""

        for n, path in enumerate(glob.glob(self.path)):
            df = _read_narrowPeak(path)
            df["sample"] = name = os.path.basename(path).split(".")[0]
            self.DataFrames[name] = df

        self.df = (
            pd.concat(list(self.DataFrames.values()))[columns]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        self.df.columns = colnames


def _read_multi_narrowPeaks(
    narrowPeak_path,
    columns=["chrom", "chromStart", "chromEnd"],
    colnames=["chr", "start", "end"],
):

    """"""

    nPeaks = narrowPeaks(narrowPeak_path)
    nPeaks.read_group()

    return nPeaks.df