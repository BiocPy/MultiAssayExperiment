from collections import OrderedDict
from typing import Any, Dict

from ..MultiAssayExperiment import MultiAssayExperiment, _create_smap_from_experiments

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def make_mae(experiments: Dict[str, Any]) -> MultiAssayExperiment:
    """Create an :py:class:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment` from a dictionary of
    experiment objects. Each experiment is either an :py:class:`~anndata.AnnData` object or a subclass of
    :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`. :py:class:`~anndata.AnnData` objects
    will be converted to a :py:class:`~singlecellexperiment.SingleCellExperiment.SingleCellExperiment`.

    The import naively creates sample mapping, with each ``experiment`` considered to be a
    independent `sample`. We add a sample to
    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`
    in this pattern - ``unknown_sample_{experiment_name}``. All cells from the same experiment are
    considered to be from the same sample and is reflected in
    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_map`.

    Args:
        experiments:
            A dictionary of experiments with experiment names as keys and the
            experiments as values.

            Each ``experiment`` can be either a :py:class:`~anndata.AnnData` object or a
            subclass of :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`.

    Raises:
        TypeError:
            - If any of the provided objects are not an expected types.
            - If ``experiments`` is not a dictionary.

    Returns:
        An MAE from the experiments.
    """
    from singlecellexperiment import SingleCellExperiment
    from anndata import AnnData
    from summarizedexperiment import SummarizedExperiment

    if not isinstance(experiments, dict):
        raise TypeError("'experiments' is not a dictionary.")

    failedExpts = []
    for expname, expt in experiments.items():
        if not (
            isinstance(expt, AnnData) or issubclass(type(expt), SummarizedExperiment)
        ):
            failedExpts.append(expname)

    if len(failedExpts) > 0:
        raise TypeError(
            f"Experiments '{', '.join(failedExpts)}' are not compatible, Must be either an "
            "AnnData, or a subclass derived from `SummarizedExperiment`."
        )

    newExpts = OrderedDict()
    for expname, expt in experiments.items():
        if isinstance(expt, AnnData):
            newExpts[expname] = SingleCellExperiment.from_anndata(expt)
        else:
            newExpts[expname] = expt

    col_data, sample_map = _create_smap_from_experiments(newExpts)

    return MultiAssayExperiment(
        experiments=newExpts,
        column_data=col_data,
        sample_map=sample_map,
    )
