import pytest
import multiassayexperiment as mae

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_MAE_fromH5AD():
    tse = mae.readH5AD("tests/data/adata.h5ad")

    assert tse is not None
    assert isinstance(tse, mae.MultiAssayExperiment)

    assert tse.experiment is not None
    assert tse.sampleMap is not None
    assert tse.colData is not None
