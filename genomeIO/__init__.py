# __init__.py

from ._xlsx._df_to_excel import _write_multi_df_to_excel as to_excel
from ._read_macs2._read_narrowPeak import _read_narrowPeak as read_narrowPeak
from ._read_macs2._read_narrowPeak import _read_multi_narrowPeaks as read_multi_narrowPeaks
from ._read_macs2._get_peaks_near_gene import _get_peaks_near_gene as peaks_near_gene