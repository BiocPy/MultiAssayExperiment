from random import random

import genomicranges
import numpy as np
import pandas as pd
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
    index=[f"sce_{i}" for i in range(6)],
)
column_data_se = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=[f"se_{i}" for i in range(6)],
)

sample_map = pd.DataFrame(
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

sample_data = pd.DataFrame(
    {"samples": ["sample1", "sample2"]}, index=["sample1", "sample2"]
)

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


def test_access_expt_with_column_data():
    assert mae is not None

    se = mae.experiment("se")
    assert se.shape == tse2.shape

    sce = mae.experiment("sce", with_sample_data=True)
    assert sce.shape == tsce.shape

    assert len(sce.column_data.columns) >= len(tsce.column_data.columns)


def test_access_expt_with_int_index():
    assert mae is not None

    se = mae.experiment(0)
    assert se.shape == tse2.shape

    sce = mae.experiment(1, with_sample_data=True)
    assert sce.shape == tsce.shape

    assert len(sce.column_data.columns) >= len(tsce.column_data.columns)
