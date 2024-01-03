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
    index=["sce"] * 6,
)

sample_map_sce = pd.DataFrame(
    {
        "assay": ["sce"] * 6,
        "primary": ["sample1"] * 6,
        "colname": ["sce"] * 6,
    }
)

sample_column_data_sce = pd.DataFrame({"samples": ["sample1"]}, index=["sample1"])

column_data_se = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=["se"] * 6,
)

sample_map_se = pd.DataFrame(
    {
        "assay": ["se"] * 6,
        "primary": ["sample2"] * 6,
        "colname": ["se"] * 6,
    }
)

sample_column_data_se = pd.DataFrame({"samples": ["sample2"]}, index=["sample2"])


def test_MAE_addExpt():
    tsce = SingleCellExperiment(
        assays={"counts": counts}, row_data=df_gr, column_data=column_data_sce
    )

    mae = MultiAssayExperiment(
        experiments={"sce": tsce},
        column_data=sample_column_data_sce,
        sample_map=sample_map_sce,
        metadata={"could be": "anything"},
    )

    assert mae is not None
    assert isinstance(mae, MultiAssayExperiment)

    tse2 = SummarizedExperiment(
        assays={"counts": counts.copy()},
        row_data=df_gr.copy(),
        column_data=column_data_se,
    )

    mae.add_experiment(
        name="se",
        experiment=tse2,
        sample_map=sample_map_se,
        column_data=sample_column_data_se,
    )

    assert mae is not None
