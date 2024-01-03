import anndata
import numpy as np
from anndata import AnnData
from mudata import MuData

import multiassayexperiment as mae

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_MAE_fromH5AD():
    tse = mae.read_h5ad("tests/data/adata.h5ad")

    assert tse is not None
    assert isinstance(tse, mae.MultiAssayExperiment)

    assert tse.experiments is not None
    assert tse.sample_map is not None
    assert tse.column_data is not None


# credit: MuData docs
def test_MAE_from_mudata():
    np.random.seed(1)

    n, d, k = 1000, 100, 10

    z = np.random.normal(loc=np.arange(k), scale=np.arange(k) * 2, size=(n, k))
    w = np.random.normal(size=(d, k))
    y = np.dot(z, w.T)

    adata = AnnData(y)
    adata.obs_names = [f"obs_{i+1}" for i in range(n)]
    adata.var_names = [f"var_{j+1}" for j in range(d)]

    d2 = 50
    w2 = np.random.normal(size=(d2, k))
    y2 = np.dot(z, w2.T)

    adata2 = AnnData(y2)
    adata2.obs_names = [f"obs_{i+1}" for i in range(n)]
    adata2.var_names = [f"var2_{j+1}" for j in range(d2)]

    mdata = MuData({"A": adata, "B": adata2})

    muMAE = mae.MultiAssayExperiment.from_mudata(input=mdata)

    assert muMAE is not None
    assert isinstance(muMAE, mae.MultiAssayExperiment)


def test_MAE_make_mae():
    np.random.seed(1)

    n, d, k = 1000, 100, 10

    z = np.random.normal(loc=np.arange(k), scale=np.arange(k) * 2, size=(n, k))
    w = np.random.normal(size=(d, k))
    y = np.dot(z, w.T)

    adata = AnnData(y)
    adata.obs_names = [f"obs_{i+1}" for i in range(n)]
    adata.var_names = [f"var_{j+1}" for j in range(d)]

    d2 = 50
    w2 = np.random.normal(size=(d2, k))
    y2 = np.dot(z, w2.T)

    adata2 = AnnData(y2)
    adata2.obs_names = [f"obs_{i+1}" for i in range(n)]
    adata2.var_names = [f"var2_{j+1}" for j in range(d2)]

    adata3 = anndata.read_h5ad("tests/data/adata.h5ad")

    muMAE = mae.make_mae(experiments={"rna": adata, "spatial": adata2, "multi": adata3})

    assert muMAE is not None
    assert isinstance(muMAE, mae.MultiAssayExperiment)

    assert muMAE.experiments is not None
    assert muMAE.sample_map is not None
    assert muMAE.column_data is not None

    assert len(set(muMAE.sample_map["assay"])) == 3
    assert len(set(muMAE.sample_map["primary"])) == 3
