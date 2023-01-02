from typing import Union, MutableMapping

from collections import OrderedDict

import anndata
import singlecellexperiment as sce
from summarizedexperiment import SummarizedExperiment, RangeSummarizedExperiment
import pandas as pd

from ..MultiAssayExperiment import MultiAssayExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def makeMAE(
    experiments: MutableMapping[
        str,
        Union[
            anndata.AnnData,
            sce.SingleCellExperiment,
            SummarizedExperiment,
            RangeSummarizedExperiment,
        ],
    ]
) -> MultiAssayExperiment:
    """Make MAE from list of experiments. 
        naively creates sample map and coldata objects.
        Also converts `AnnData` objects to `SingleCellExperiment`s.

    Args:
        experiments (MutableMapping[str, Union[anndata.AnnData, SingleCellExperiment, SummarizedExperiment, RangeSummarizedExperiment]]): a dictionary of experiments. 

    Raises:
        TypeError: if any of the provided objects are not an expected types.

    Returns:
        MultiAssayExperiment: an MAE from the experiments
    """
    failedExpts = []
    for expname, expt in experiments.items():
        if not (
            isinstance(expt, sce.SingleCellExperiment)
            or isinstance(expt, anndata.AnnData)
            or isinstance(expt, SummarizedExperiment)
            or isinstance(expt, RangeSummarizedExperiment)
        ):
            failedExpts.append(expname)

    if len(failedExpts) > 0:
        raise TypeError(
            f"Experiments {failedExpts} are not compatible, Must be either a AnnData, SingleCellExperiment, SummarizedExperiment objects"
        )

    newExpts = OrderedDict()

    sampleMap = pd.DataFrame()
    samples = []

    for expname, expt in experiments.items():
        if isinstance(expt, anndata.AnnData):
            newExpts[expname] = sce.fromAnnData(expt)
        else:
            newExpts[expname] = expt

        colnames = newExpts[expname].colnames
        asy_sample = f"unknown_sample_{expname}"
        asy_df = pd.DataFrame(
            {
                "assay": [expname] * len(colnames),
                "primary": [asy_sample] * len(colnames),
                "colname": colnames,
            }
        )

        sampleMap = pd.concat([sampleMap, asy_df])
        samples.append(asy_sample)

    coldata = pd.DataFrame({"samples": samples}, index=samples)

    return MultiAssayExperiment(
        experiments=newExpts, colData=coldata, sampleMap=sampleMap,
    )
