from random import random

import genomicranges
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from biocframe import BiocFrame
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

import multiassayexperiment
from multiassayexperiment import MultiAssayExperiment

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

gr = genomicranges.GenomicRanges.from_pandas(df_gr)

column_data_sce = BiocFrame(
    {"treatment": ["ChIP", "Input"] * 3},
    row_names=[f"sce_{i}" for i in range(6)],
)

column_data_se = BiocFrame(
    {"treatment": ["ChIP", "Input"] * 3},
    row_names=[f"se_{i}" for i in range(6)],
)

sample_map = BiocFrame(
    {
        "assay": ["sce", "se"] * 6,
        "primary": ["sample1", "sample2"] * 6,
        "colname": [
            "sce_0",
            "se_0",
            "sce_1",
            "se_1",
            "sce_2",
            "se_2",
            "sce_3",
            "se_3",
            "sce_4",
            "se_4",
            "sce_5",
            "se_5",
        ],
    }
)

sample_data = BiocFrame(
    {"samples": ["sample1", "sample2"]}, row_names=["sample1", "sample2"]
)


def test_MAE_creation():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, row_data=df_gr, column_data=column_data_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        row_data=df_gr.copy(),
        column_data=column_data_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        column_data=sample_data,
        sample_map=sample_map,
        metadata={"could be": "anything"},
    )

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)


def test_MAE_creation_with_alts():
    tse = SummarizedExperiment(
        assays={"counts": counts}, row_data=df_gr, column_data=column_data_se
    )

    tsce = SingleCellExperiment(
        assays={"counts": counts},
        row_data=df_gr,
        column_data=column_data_sce,
        alternative_experiments={"alt": tse},
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        row_data=df_gr.copy(),
        column_data=column_data_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        column_data=sample_data,
        sample_map=sample_map,
        metadata={"could be": "anything"},
    )

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)

    assert mae.experiments is not None
    assert mae.experiment("sce") is not None
    assert mae.assays is not None
    assert mae.column_data is not None
    assert mae.sample_map is not None

    with pytest.raises(Exception):
        mae.column_data = None

    with pytest.raises(Exception):
        mae.sample_map = None

    assert mae.metadata is not None
    mae.metadata = {}
    assert mae.metadata is not None


def test_MAE_completedcases():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, row_data=df_gr, column_data=column_data_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        row_data=df_gr.copy(),
        column_data=column_data_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        column_data=sample_data,
        sample_map=sample_map,
        metadata={"could be": "anything"},
    )

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)

    completed = mae.complete_cases()

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

    mae = multiassayexperiment.make_mae(experiments={"rna": adata, "spatial": adata2})

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)

    repls = mae.replicated()

    assert repls is not None
    assert len(repls) == len(mae.experiments.keys())


def test_with_sample_data():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, row_data=gr.to_pandas(), column_data=column_data_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        row_data=gr.to_pandas().copy(),
        column_data=column_data_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        column_data=sample_data,
        sample_map=sample_map,
        metadata={"could be": "anything"},
    )

    expt_with_sample_data = mae.experiment("se", with_sample_data=True)

    assert expt_with_sample_data is not None
    assert expt_with_sample_data.column_data is not None
    print(expt_with_sample_data.column_data)
    assert expt_with_sample_data.column_data.get_column("samples") is not None


def test_MAE_intersect_methods():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, row_data=df_gr, column_data=column_data_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        row_data=df_gr.copy(),
        column_data=column_data_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        column_data=sample_data,
        sample_map=sample_map,
        metadata={"could be": "anything"},
    )

    row_mae = mae.intersect_rows()
    assert row_mae is not None
    assert isinstance(row_mae, MultiAssayExperiment)
    assert row_mae.experiment_names == mae.experiment_names


def test_MAE_names():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, row_data=df_gr, column_data=column_data_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        row_data=df_gr.copy(),
        column_data=column_data_se.copy(),
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce, "se": tse2},
        column_data=sample_data,
        sample_map=sample_map,
        metadata={"could be": "anything"},
    )

    rownames = mae.rownames
    colnames = mae.column_names

    assert rownames is not None
    assert len(rownames) == len(mae.experiments)

    assert colnames is not None
    assert len(colnames) == len(mae.experiments)
