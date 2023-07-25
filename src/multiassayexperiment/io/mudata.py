from collections import OrderedDict

import mudata
import pandas as pd
from singlecellexperiment import fromAnnData

from ..MultiAssayExperiment import MultiAssayExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def fromMuData(mudata: mudata.MuData) -> MultiAssayExperiment:
    """Transform MuData object to MAE representation.

    Args:
        mudata (mudata.MuData): MuData object.

    Returns:
        MultiAssayExperiment: MAE representation.
    """

    if mudata.isbacked:
        raise Exception("backed mode is currently not supported")

    experiments = OrderedDict()

    sampleMap = pd.DataFrame()
    samples = []

    for asy, adata in mudata.mod.items():
        experiments[asy] = fromAnnData(adata)

        colnames = None
        if adata.obs.index.tolist() is not None:
            colnames = adata.obs.index.tolist()
        else:
            colnames = range(len(adata.shape[0]))

        asy_sample = f"unknown_sample_{asy}"

        asy_df = pd.DataFrame(
            {
                "assay": [asy] * len(colnames),
                "primary": [asy_sample] * len(colnames),
                "colname": colnames,
            }
        )

        sampleMap = pd.concat([sampleMap, asy_df])
        samples.append(asy_sample)

    coldata = pd.DataFrame({"samples": samples}, index=samples)

    return MultiAssayExperiment(
        experiments=experiments,
        colData=coldata,
        sampleMap=sampleMap,
        metadata=mudata.uns,
    )
