from random import random

import genomicranges
import numpy as np
import pandas as pd
import pytest
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

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

column_data_sce = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=["sce"] * 6,
)
column_data_se = pd.DataFrame(
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


def test_MAE_creation_fails():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, row_data=df_gr, column_data=column_data_sce
    )

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        row_data=df_gr.copy(),
        column_data=column_data_sce.copy(),
    )

    with pytest.raises(Exception):
        MultiAssayExperiment(
            experiments={"sce": tsce, "se": tse2},
            column_data=sample_data,
            sample_map=sample_map,
            metadata={"could be": "anything"},
        )


def test_MAE_save():
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

    mdata = mae.to_mudata()

    assert mdata is not None
    assert len(mdata.mod.keys()) == 2


def test_empty_mae():
    mae = MultiAssayExperiment(experiments={})

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)
