import anndata
import singlecellexperiment
from anndata import AnnData
from pandas import DataFrame

from ..MultiAssayExperiment import MultiAssayExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def from_anndata(adata: AnnData, name: str = "unknown") -> MultiAssayExperiment:
    """Read :py:class:`~anndata.AnnData` objects as a
    :py:class:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment`.

    Since :py:class:`~anndata.AnnData` does not contain sample information,
    sample named ``"unknown_sample"`` will be added to
    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`.
    All cells are considered to be extracted from this sample and is reflected in
    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_map`.

    Args:
        data (AnnData): An `AnnData` object.
        name (str, optional): Name for the experiment. Defaults to "unknown".

    Returns:
        MultiAssayExperiment: An MAE of ``data``.
    """

    if not isinstance(adata, AnnData):
        raise TypeError("data is not an `AnnData` object.")

    scexpt = singlecellexperiment.from_anndata(adata=adata)

    experiments = {name: scexpt}

    col_data = DataFrame({"samples": ["unknown_sample"]}, index=["unknown_sample"])

    sample_map = DataFrame()
    colnames = None
    if adata.obs.index.tolist() is not None:
        colnames = adata.obs.index.tolist()
    else:
        colnames = range(len(adata.shape[0]))
    sample_map["colname"] = colnames
    sample_map["assay"] = "unknown"
    sample_map["primary"] = "unknown_sample"

    return MultiAssayExperiment(
        experiments=experiments,
        col_data=col_data,
        sample_map=sample_map,
        metadata=adata.uns,
    )


def read_h5ad(path: str) -> MultiAssayExperiment:
    """Read a H5ad file as a
    :py:class:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment`.

    This function reads the h5ad at the ``path`` using :py:func:`~anndata.read_h5ad` and converts
    it into an MAE using :py:func:`~multiassayexperiment.io.anndata.from_anndata`.

    Args:
        path (str): Path to a H5AD file

    Returns:
        MultiAssayExperiment: An MAE from the H5ad file.
    """

    adata = anndata.read_h5ad(path)
    return from_anndata(adata)
