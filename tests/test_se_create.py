import pytest

from multiassayexperiment import MultiAssayExperiment
from singlecellexperiment import SingleCellExperiment
import numpy as np
from random import random
import pandas as pd
from genomicranges import GenomicRanges
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

gr = GenomicRanges.fromPandas(df_gr)

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

sample_data = pd.DataFrame({"samples": ["sample1", "sample2"]})


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
        alterExps={"alt": tse},
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


def test_MAE_creation_fails():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, rowData=df_gr, colData=colData_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        rowData=df_gr.copy(),
        colData=colData_sce.copy(),
    )

    with pytest.raises(Exception):
        mae = MultiAssayExperiment(
            experiments={"sce": tsce, "se": tse2},
            colData=sample_data,
            sampleMap=sample_map,
            metadata={"could be": "anything"},
        )
