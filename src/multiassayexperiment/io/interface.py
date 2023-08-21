from collections import OrderedDict
from typing import MutableMapping, Union

import pandas as pd
import singlecellexperiment as sce
from anndata import AnnData
from summarizedexperiment import SummarizedExperiment

from ..MultiAssayExperiment import MultiAssayExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def make_mae(
    experiments: MutableMapping[
        str,
        Union[AnnData, SummarizedExperiment],
    ]
) -> MultiAssayExperiment:
    """Read a dictionary of experiments into an
    :py:class:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment`.

    The import naively creates sample mapping, each ``experiment`` is considered to be a `sample`.
    We add a sample with the following pattern - ``unknown_sample_{experiment_name}`` to
    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`
    All cells from the same experiment are extracted from the same sample and should reflect in
    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_mapping`.

    Additionally, converts :py:class`~anndata.AnnData` objects to
    :py:class:`~singlecellexperiment.SingleCellExperiment.SingleCellExperiment` objects.

    Args:
        experiments (MutableMapping[str, Union[AnnData, SummarizedExperiment]]): A dictionary of
            experiments with experiment names as keys and the experiments as values.

            each ``experiment`` can be represented as :py:class`~anndata.AnnData` objects or any
            subclass of :py:class`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`.

    Raises:
        TypeError: If any of the provided objects are not an expected types.
        TypeError: If ``experiments`` is not a dictionary.

    Returns:
        MultiAssayExperiment: An MAE from the experiments.
    """

    if not isinstance(experiments, dict):
        raise TypeError("`experiments` is not a dictionary.")

    failedExpts = []
    for expname, expt in experiments.items():
        if not (
            isinstance(expt, AnnData) or issubclass(type(expt), SummarizedExperiment)
        ):
            failedExpts.append(expname)

    if len(failedExpts) > 0:
        raise TypeError(
            f"Experiments '{', '.join(failedExpts)}' are not compatible, Must be either an "
            "AnnData, or a subclass derived from SummarizedExperiment."
        )

    newExpts = OrderedDict()

    sample_map = pd.DataFrame()
    samples = []

    for expname, expt in experiments.items():
        if isinstance(expt, AnnData):
            newExpts[expname] = sce.from_anndata(expt)
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

        sample_map = pd.concat([sample_map, asy_df])
        samples.append(asy_sample)

    col_data = pd.DataFrame({"samples": samples}, index=samples)

    return MultiAssayExperiment(
        experiments=newExpts,
        col_data=col_data,
        sample_map=sample_map,
    )
