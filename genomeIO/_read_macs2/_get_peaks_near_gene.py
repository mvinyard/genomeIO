import seq_toolkit as seq
import anngtf

from ._read_narrowPeak import _read_multi_narrowPeaks

def _fetch_gene_TSS(gene_strand, gene_bounds):

    """"""

    chrom = gene_bounds["seqname"].values[0]
    start = gene_bounds["start"].values[0]
    end = gene_bounds["end"].values[0]

    if gene_strand == "+":
        return [chrom, start]
    else:
        return [chrom, end]


def _get_window_df(gtf, Gene, chrom, Window):

    Window["df"] = (
        gtf.loc[gtf["seqname"] == Gene["TSS"][0]]
        .loc[gtf["start"] > Window["start"]]
        .loc[gtf["end"] > Window["end"]]
    )

    Window["genes_in_window"] = (
        Window["df"]
        .loc[Window["df"]["feature"] == "gene"]["gene_name"]
        .unique()
        .tolist()
    )
    Window["n_genes"] = len(Window["genes_in_window"])

    return Window


class GeneRegulatoryRegion:
    def __init__(self):

        self._gtf = anngtf.load()
        self.Gene = {}

    def get_gene(self, gene):

        self.Gene["gene"] = self._gtf.loc[self._gtf["feature"] == "gene"].loc[
            self._gtf["gene_name"] == gene
        ]
        self.Gene["strand"] = self.Gene["gene"]["strand"].unique()[0]
        self.Gene["bounds"] = self.Gene["gene"][["seqname", "start", "end"]]
        self.Gene["TSS"] = _fetch_gene_TSS(self.Gene["strand"], self.Gene["bounds"])

    def set_window(self, window_size=[200_000, 200_000]):

        """"""

        self.Window = {}
        self.Window["chrom"] = self.Gene["TSS"][0]
        self.Window["start"] = (window_size[0] * -1) + self.Gene["TSS"][1].astype(int)
        self.Window["end"] = window_size[1] + self.Gene["TSS"][1].astype(int)
        self.Window = _get_window_df(
            self._gtf, self.Gene, self.Window["chrom"], self.Window
        )


def _fetch_gene_regulatory_region_surrounding_gene(
    gene, window_start=200_000, window_end=200_000
):

    """"""

    reg = GeneRegulatoryRegion()
    reg.get_gene(gene)
    reg.set_window(window_size=[window_start, window_end])

    return reg.Gene, reg.Window


def _get_peaks_in_windows(peak_df, WindowDict):

    """"""

    chrom = WindowDict["chrom"]
    window_start = WindowDict["start"]
    window_end = WindowDict["end"]

    window_peak_df = (
        peak_df.loc[peak_df["chr"] == chrom]
        .loc[peak_df["start"] > window_start]
        .loc[peak_df["end"] < window_end]
    )

    return window_peak_df

def _get_peaks_near_gene(narrowPeak_path, gene, window_start=200_000, window_end=200_000):
    
    PeakGene = {}
    
    PeakGene['peak_df'] = peak_df = _read_multi_narrowPeaks(narrowPeak_path)
    PeakGene['gene'], PeakGene['window'] = gene, window = _fetch_gene_regulatory_region_surrounding_gene(gene, 
                                                                                                         window_start, 
                                                                                                         window_end)
    PeakGene['filtered_peak_df'] = _get_peaks_in_windows(peak_df, window)
    PeakGene['filtered_peak_df'].columns = ["Chromosome", "Start", "End"]
    PeakGene['filtered_peak_df'] = PeakGene['filtered_peak_df'].reset_index(drop=True)
    peak_GF = seq.GenomicFeatures(PeakGene['filtered_peak_df'])
    peak_GF.merge()
    
    PeakGene['merged_peak_df'] = peak_GF.merged_df
    
    return PeakGene