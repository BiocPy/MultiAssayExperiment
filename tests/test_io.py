import anndata
import multiassayexperiment as mae
import numpy as np
from anndata import AnnData
from mudata import MuData

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


def test_MAE_fromH5AD():
    tse = mae.readH5AD("tests/data/adata.h5ad")

    assert tse is not None
    assert isinstance(tse, mae.MultiAssayExperiment)

    assert tse.experiments is not None
    assert tse.sampleMap is not None
    assert tse.colData is not None


# credit: MuData docs
def test_MAE_fromMuData():
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

    muMAE = mae.fromMuData(mudata=mdata)

    assert muMAE is not None
    assert isinstance(muMAE, mae.MultiAssayExperiment)


def test_MAE_makeMAE():
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

    muMAE = mae.makeMAE(experiments={"rna": adata, "spatial": adata2, "multi": adata3})

    assert muMAE is not None
    assert isinstance(muMAE, mae.MultiAssayExperiment)

    assert muMAE.experiments is not None
    assert muMAE.sampleMap is not None
    assert muMAE.colData is not None

    assert len(muMAE.sampleMap["assay"].unique()) == 3
    assert len(muMAE.sampleMap["primary"].unique()) == 3
