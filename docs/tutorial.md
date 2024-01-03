# Tutorial

Container class to represent and manage multi-omics genomic experiments.

For more detailed description checkout the [MultiAssayExperiment Bioc/R package](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html))

# Construct an `MultiAssayExperiment`

An MAE contains three main entities,

- Primary information (`col_data`): Bio-specimen/sample information. The ``col_data`` may provide information about patients, cell lines, or
  other biological units.
- Experiments (`experiments`): Genomic data from each experiment. either a `SingleCellExperiment`, `SummarizedExperiment`, `RangeSummarizedExperiment` or
  any class that extends a `SummarizedExperiment`.
- Sample Map (`sample_map`): Map biological units from ``col_data`` to the list of ``experiments``. Must contain columns,

  - **assay** provides the names of the different experiments performed on the
    biological units. All experiment names from ``experiments`` must be present in this column.
  - **primary** contains the sample name. All names in this column must match with row labels from ``col_data``.
  - **colname** is the mapping of samples/cells within each experiment back to its biosample information in ``col_data``.

Lets create these objects

```python
from biocframe import BiocFrame
from iranges import IRanges
import numpy as np
from genomicranges import GenomicRanges
from random import random

nrows = 200
ncols = 6
counts = np.random.rand(nrows, ncols)
gr = GenomicRanges(
    seqnames=[
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
        ] * 20,
    ranges=IRanges(range(100, 300), range(110, 310)),
    strand = ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"] * 20,
    mcols=BiocFrame({
        "score": range(0, 200),
        "GC": [random() for _ in range(10)] * 20,
    })
)

col_data_sce = BiocFrame({"treatment": ["ChIP", "Input"] * 3},
    row_names=["sce"] * 6,
)

col_data_se = BiocFrame({"treatment": ["ChIP", "Input"] * 3},
    row_names=["se"] * 6,
)

sample_map = BiocFrame({
    "assay": ["sce", "se"] * 6,
    "primary": ["sample1", "sample2"] * 6,
    "colname": ["sce", "se"] * 6
})

sample_data = BiocFrame({"samples": ["sample1", "sample2"]}, row_names=["sample1", "sample2"])
```

Then, create various experiment classes,

```python
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

tsce = SingleCellExperiment(
    assays={"counts": counts}, row_data=gr.to_pandas(), column_data=col_data_sce
)

tse2 = SummarizedExperiment(
    assays={"counts": counts.copy()},
    row_data=gr.to_pandas().copy(),
    column_data=col_data_se.copy(),
)
```

Now that we have all the pieces together, we can now create an MAE,

```python
from multiassayexperiment import MultiAssayExperiment

mae = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    column_data=sample_data,
    sample_map=sample_map,
    metadata={"could be": "anything"},
)
```

To make your life easier, we also provide methods to naively create sample mapping from experiments.

**_This is not a recommended approach, but if you don't have sample mapping, then it doesn't matter._**

```python
import multiassayexperiment
maeObj = multiassayexperiment.make_mae(experiments={"sce": tsce, "se": tse2})
```

## Import `MuData` and `AnnData` as `MultiAssayExperiment`

If you have datasets stored as `MuData`, these can be easily converted to an MAE using the `from_mudata` method.

Lets first construct `AnnData`` objects and then an MAE

```python
import multiassayexperiment as mae
import numpy as np
from anndata import AnnData

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
```

we can now construct a `MuData` object and convert that to an MAE

```python
from mudata import MuData
from multiassayexperiment import MultiAssayExperiment
mdata = MuData({"rna": adata, "spatial": adata2})

maeObj = MultiAssayExperiment.from_mudata(input=mdata)
```

Methods are also available to convert an `AnnData` object to `MAE`.

```python
import multiassayexperiment
maeObj = multiassayexperiment.read_h5ad("tests/data/adata.h5ad")
```

# Accessors

Multiple methods are available to access various slots of a `MultiAssayExperiment` object

```python
mae.assays
mae.column_data
mae.sample_map
mae.experiments
mae.metadata
```

## Access experiments

if you want to access a specific experiment

```python
# access a specific experiment
mae.experiment("se")
```

This does not include the sample data stored in the MAE. If you want to include this information

***Note: This creates a copy of the experiment object.***

```python
expt_with_sampleData = maeObj.experiment(experiment_name, with_sample_data=True)
```

# Slice a `MultiAssayExperiment`

`MultiAssayExperiment` allows subsetting by `rows`, `columns`, and `experiments`. Samples are automatically sliced during this operation.

The structure for slicing,

```
mae[rows, columns, experiments]
```

- rows, columns: accepts either a slice, list of indices or a dictionary to specify slices per experiment.
- experiments: accepts a list of experiment names to subset to.

## Slice by row and column slices

```python
maeObj[1:5, 0:4]
```

## Slice by rows, columns, experiments

```python
maeObj[1:5, 0:4, ["spatial"]]
```

Checkout other methods that perform similar operations - `subset_by_rows`, `subset_by_columns` & `subset_by_experiments`.

# Helper methods

## completedCases

This method returns a boolean vector that specifies which bio specimens have data across all experiments.

```python
maeObj.completed_cases()
```

## replicated

replicated identifies bio specimens that have multiple observations per experiment.

```python
maeObj.replicated()
```
