from collections import OrderedDict

from mudata import MuData
from pandas import DataFrame, concat
from singlecellexperiment import from_anndata

from ..MultiAssayExperiment import MultiAssayExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def from_mudata(mudata: MuData) -> MultiAssayExperiment:
    """Read :py:class:`~mudata.MuData` as
    :py:class:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment`.

    The import naively creates sample mapping, each ``experiment`` is considered to be a `sample`.
    We add a sample with the following pattern - ``"unknown_sample_{experiment_name}"`` to
    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`
    All cells from the same experiment are considered to be extracted from the same sample and is
    reflected in
    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_map`.

    Args:
        mudata (MuData): MuData object.

    Raises:
        Exception: If ``mudata`` object is read in backed mode :py:attr:`~mudata.MuData.isbacked`.

    Returns:
        MultiAssayExperiment: MAE representation.
    """

    if mudata.isbacked is True:
        raise Exception("backed mode is currently not supported.")

    experiments = OrderedDict()

    sample_map = DataFrame()
    samples = []

    for asy, adata in mudata.mod.items():
        experiments[asy] = from_anndata(adata)

        colnames = None
        if adata.obs.index.tolist() is not None:
            colnames = adata.obs.index.tolist()
        else:
            colnames = range(len(adata.shape[0]))

        asy_sample = f"unknown_sample_{asy}"

        asy_df = DataFrame(
            {
                "assay": [asy] * len(colnames),
                "primary": [asy_sample] * len(colnames),
                "colname": colnames,
            }
        )

        sample_map = concat([sample_map, asy_df])
        samples.append(asy_sample)

    col_data = DataFrame({"samples": samples}, index=samples)

    return MultiAssayExperiment(
        experiments=experiments,
        col_data=col_data,
        sample_map=sample_map,
        metadata=mudata.uns,
    )
