import anndata
import numpy as np
from anndata import AnnData

import multiassayexperiment as mae

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

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


def test_MAE_slice():
    muMAE = mae.make_mae(experiments={"rna": adata, "spatial": adata2, "multi": adata3})

    assert muMAE is not None
    assert isinstance(muMAE, mae.MultiAssayExperiment)

    assert muMAE.experiments is not None
    assert muMAE.sample_map is not None
    assert muMAE.column_data is not None

    assert len(set(muMAE.sample_map["assay"])) == 3
    assert len(set(muMAE.sample_map["primary"])) == 3

    sliced_MAE = muMAE[1:3, 1:2]
    assert sliced_MAE is not None
    assert isinstance(sliced_MAE, mae.MultiAssayExperiment)

    assert sliced_MAE.experiments is not None
    assert sliced_MAE.sample_map is not None
    assert sliced_MAE.column_data is not None

    assert len(set(sliced_MAE.sample_map["assay"])) == 1
    assert len(set(sliced_MAE.sample_map["primary"])) == 1
    assert sliced_MAE.sample_map.shape[0] != muMAE.sample_map.shape[0]
    assert sliced_MAE.sample_map.shape[0] == 1000

    sliced_MAE_assay = muMAE[None, None, ["rna", "spatial"]]
    assert sliced_MAE_assay is not None
    assert isinstance(sliced_MAE_assay, mae.MultiAssayExperiment)

    assert sliced_MAE_assay.experiments is not None
    assert sliced_MAE_assay.sample_map is not None
    assert sliced_MAE_assay.column_data is not None

    assert len(set(sliced_MAE_assay.sample_map["assay"])) == 2
    assert len(set(sliced_MAE_assay.sample_map["primary"])) == 2
    assert sliced_MAE_assay.sample_map.shape[0] == 2000

    sliced_MAE_assay = muMAE[1:3, 0:5, ["rna"]]
    assert sliced_MAE_assay is not None
    assert isinstance(sliced_MAE_assay, mae.MultiAssayExperiment)

    assert sliced_MAE_assay.experiments is not None
    assert sliced_MAE_assay.sample_map is not None
    assert sliced_MAE_assay.column_data is not None

    assert len(set(sliced_MAE_assay.sample_map["assay"])) == 1
    assert len(set(sliced_MAE_assay.sample_map["primary"])) == 1
    assert sliced_MAE_assay.sample_map.shape[0] == 1000


# def test_MAE_slice_dict():
#     muMAE = mae.make_mae(experiments={"rna": adata, "spatial": adata2, "multi": adata3})

#     sliced_MAE_assay = muMAE[
#         {"rna": slice(0, 5)}, {"spatial": slice(0, 10)}, ["rna", "spatial"]
#     ]
#     assert sliced_MAE_assay is not None
#     assert isinstance(sliced_MAE_assay, mae.MultiAssayExperiment)

#     assert sliced_MAE_assay.experiments is not None
#     assert sliced_MAE_assay.sample_map is not None
#     assert sliced_MAE_assay.column_data is not None

#     assert len(sliced_MAE_assay.sample_map["assay"].unique()) == 2
#     assert len(sliced_MAE_assay.sample_map["primary"].unique()) == 2
#     assert sliced_MAE_assay.sample_map.shape[0] == 1010


def test_MAE_subset_by_row():
    muMAE = mae.make_mae(experiments={"rna": adata, "spatial": adata2, "multi": adata3})

    assert muMAE is not None
    assert isinstance(muMAE, mae.MultiAssayExperiment)

    assert muMAE.experiments is not None
    assert muMAE.sample_map is not None
    assert muMAE.column_data is not None

    assert len(set(muMAE.sample_map["assay"])) == 3
    assert len(set(muMAE.sample_map["primary"])) == 3

    sliced_MAE = muMAE.subset_by_row(rows=[10, 2, 5])
    assert sliced_MAE is not None
    assert isinstance(sliced_MAE, mae.MultiAssayExperiment)

    assert sliced_MAE.experiments is not None
    assert sliced_MAE.sample_map is not None
    assert sliced_MAE.column_data is not None

    assert len(set(sliced_MAE.sample_map["assay"])) == 3
    assert len(set(sliced_MAE.sample_map["primary"])) == 3
    assert sliced_MAE.sample_map.shape == (2030, 3)


def test_MAE_subset_by_column():
    muMAE = mae.make_mae(experiments={"rna": adata, "spatial": adata2, "multi": adata3})

    assert muMAE is not None
    assert isinstance(muMAE, mae.MultiAssayExperiment)

    assert muMAE.experiments is not None
    assert muMAE.sample_map is not None
    assert muMAE.column_data is not None

    assert len(set(muMAE.sample_map["assay"])) == 3
    assert len(set(muMAE.sample_map["primary"])) == 3

    sliced_MAE = muMAE.subset_by_column(columns=[1, 2])
    assert sliced_MAE is not None
    assert isinstance(sliced_MAE, mae.MultiAssayExperiment)

    assert sliced_MAE.experiments is not None
    assert sliced_MAE.sample_map is not None
    assert sliced_MAE.column_data is not None

    assert len(set(sliced_MAE.sample_map["assay"])) == 2
    assert len(set(sliced_MAE.sample_map["primary"])) == 2
    assert sliced_MAE.sample_map.shape == (1030, 3)
    assert len(sliced_MAE.experiment_names) == len(muMAE.experiment_names)


def test_MAE_subsetByExpt():
    muMAE = mae.make_mae(experiments={"rna": adata, "spatial": adata2, "multi": adata3})

    assert muMAE is not None
    assert isinstance(muMAE, mae.MultiAssayExperiment)

    assert muMAE.experiments is not None
    assert muMAE.sample_map is not None
    assert muMAE.column_data is not None

    assert len(set(muMAE.sample_map["assay"])) == 3
    assert len(set(muMAE.sample_map["primary"])) == 3

    sliced_MAE = muMAE.subset_by_experiments(experiments=["rna", "spatial"])
    assert sliced_MAE is not None
    assert isinstance(sliced_MAE, mae.MultiAssayExperiment)

    assert sliced_MAE.experiments is not None
    assert len(sliced_MAE.experiments.keys()) == 2
    assert sliced_MAE.sample_map is not None
    assert sliced_MAE.column_data is not None

    assert len(set(sliced_MAE.sample_map["assay"])) == 2
    assert len(set(sliced_MAE.sample_map["primary"])) == 2
    assert sliced_MAE.sample_map.shape == (2000, 3)
