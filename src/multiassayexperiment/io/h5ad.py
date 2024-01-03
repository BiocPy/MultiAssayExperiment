from ..MultiAssayExperiment import MultiAssayExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def read_h5ad(path: str) -> MultiAssayExperiment:
    """Create a :py:class:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment` from a H5AD file.

    This function reads the h5ad at the ``path`` using :py:func:`~anndata.read_h5ad` and converts
    it into an MAE using :py:func:`~multiassayexperiment.io.anndata.from_anndata`.

    Args:
        path:
            Path to a H5AD file

    Returns:
        An MAE from the H5ad file.
    """
    import anndata

    adata = anndata.read_h5ad(path)
    return MultiAssayExperiment.from_anndata(adata)
