import anndata
import pandas as pd
import singlecellexperiment as sce

from ..MultiAssayExperiment import MultiAssayExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def fromAnnData(adata: anndata.AnnData) -> MultiAssayExperiment:
    """Transform AnnData object to MAE representation.

    Args:
        data (AnnData): MuData object.

    Returns:
        MultiAssayExperiment: MAE from AnnData.
    """
    scexpt = sce.fromAnnData(adata=adata)

    experiments = {"unknown": scexpt}

    coldata = pd.DataFrame({"samples": ["unknown_sample"]}, index=["unknown_sample"])

    sampleMap = pd.DataFrame()
    colnames = None
    if adata.obs.index.tolist() is not None:
        colnames = adata.obs.index.tolist()
    else:
        colnames = range(len(adata.shape[0]))
    sampleMap["colname"] = colnames
    sampleMap["assay"] = "unknown"
    sampleMap["primary"] = "unknown_sample"

    return MultiAssayExperiment(
        experiments=experiments,
        colData=coldata,
        sampleMap=sampleMap,
        metadata=adata.uns,
    )


def readH5AD(path: str) -> MultiAssayExperiment:
    """Convert H5AD to MAE representation

    Args:
        path (str): path to a H5AD file

    Returns:
        MultiAssayExperiment: mae representation
    """

    adata = anndata.read_h5ad(path)
    return fromAnnData(adata)
