from random import random

import genomicranges
import multiassayexperiment
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from multiassayexperiment import MultiAssayExperiment
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


nrows = 200
ncols = 6
counts = np.random.rand(nrows, ncols)
df_gr = pd.DataFrame(
    {
        "seqnames": [
            "chr1",
            "chr2",
            "chr2",
            "chr2",
            "chr1",
            "chr1",
            "chr3",
            "chr3",
            "chr3",
            "chr3",
        ]
        * 20,
        "starts": range(100, 300),
        "ends": range(110, 310),
        "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"] * 20,
        "score": range(0, 200),
        "GC": [random() for _ in range(10)] * 20,
    }
)

gr = genomicranges.fromPandas(df_gr)

colData_sce = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=["sce"] * 6,
)
colData_se = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=["se"] * 6,
)

sample_map = pd.DataFrame(
    {
        "assay": ["sce", "se"] * 6,
        "primary": ["sample1", "sample2"] * 6,
        "colname": ["sce", "se"] * 6,
    }
)

sample_data = pd.DataFrame(
    {"samples": ["sample1", "sample2"]}, index=["sample1", "sample2"]
)


def test_MAE_creation():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, rowData=df_gr, colData=colData_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        rowData=df_gr.copy(),
        colData=colData_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        colData=sample_data,
        sampleMap=sample_map,
        metadata={"could be": "anything"},
    )

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)


def test_MAE_creation_with_alts():
    tse = SummarizedExperiment(
        assays={"counts": counts}, rowData=df_gr, colData=colData_se
    )

    tsce = SingleCellExperiment(
        assays={"counts": counts},
        rowData=df_gr,
        colData=colData_sce,
        altExps={"alt": tse},
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        rowData=df_gr.copy(),
        colData=colData_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        colData=sample_data,
        sampleMap=sample_map,
        metadata={"could be": "anything"},
    )

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)

    assert mae.experiments is not None
    assert mae.experiment("sce") is not None
    assert mae.assays is not None
    assert mae.colData is not None
    assert mae.sampleMap is not None

    with pytest.raises(Exception):
        mae.colData = None

    with pytest.raises(Exception):
        mae.sampleMap = None

    assert mae.metadata is not None
    mae.metadata = None
    assert mae.metadata is None


def test_MAE_completedcases():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, rowData=df_gr, colData=colData_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        rowData=df_gr.copy(),
        colData=colData_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        colData=sample_data,
        sampleMap=sample_map,
        metadata={"could be": "anything"},
    )

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)

    completed = mae.completeCases()

    assert completed is not None
    assert len(completed) == len(mae.experiments.keys())
    assert completed == [False, False]


def test_MAE_replicated():
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

    mae = multiassayexperiment.makeMAE(experiments={"rna": adata, "spatial": adata2})

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)

    repls = mae.replicated()

    assert repls is not None
    assert len(repls) == len(mae.experiments.keys())
