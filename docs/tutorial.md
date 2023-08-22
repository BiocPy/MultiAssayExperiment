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
import pandas as pd
import numpy as np
from genomicranges import GenomicRanges

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

gr = GenomicRanges.from_pandas(df_gr)

col_data_sce = pd.DataFrame(
    {
        "treatment": ["ChIP", "Input"] * 3,
    },
    index=["sce"] * 6,
)

col_data_se = pd.DataFrame(
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
```

Then, create various experiment classes,

```python
from singlecellexperiment import SingleCellExperiment
from summarizedExperiment import SummarizedExperiment

tsce = SingleCellExperiment(
    assays={"counts": counts}, row_data=df_gr, col_data=col_data_sce
)

tse2 = SummarizedExperiment(
    assays={"counts": counts.copy()},
    row_data=df_gr.copy(),
    col_data=col_data_se.copy(),
)
```

Now that we have all the pieces together, we can now create an MAE,

```python
from multiassayexperiment import MultiAssayExperiment

maeObj = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    col_data=sample_data,
    sample_map=sample_map,
    metadata={"could be": "anything"},
)
```

To make your life easier, we also provide methods to naively create sample mapping from experiments.

**_This is not a recommended approach, but if you don't have sample mapping, then it doesn't matter._**

```python
maeObj = mae.make_mae(experiments={"sce": tsce, "se": tse2})
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
mdata = MuData({"rna": adata, "spatial": adata2})

maeObj = mae.from_mudata(mudata=mdata)
```

Methods are also available to convert an `AnnData` object to `MAE`.

```python
maeObj = mae.read_h5ad("tests/data/adata.h5ad")
```

# Accessors

Multiple methods are available to access various slots of a `MultiAssayExperiment` object

```python
maeObj.assays
maeObj.col_data
maeObj.sample_map
maeObj.experiments
maeObj.metadata
```

## Access experiments

if you want to access a specific experiment

```python
# access a specific experiment
maeObj.experiment(experiment_name)
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
maeObj[rows, columns, experiments]
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

## Specify slices per experiment

You can specify slices by experiment, rest of the experiments are not sliced.

```python
maeObj[{"rna": slice(0,10)}, {"spatial": slice(0,5)}, ["spatial"]]
```

Checkout other methods that perform similar operations - `subset_by_rows`, `subset_by_columns` & `subset_by_experiments`.

# Helper methods

## completedCases

This method returns a boolean vector that specifies which biospecimens have data across all experiments.

```python
maeObj.completed_cases()
```

## replicated

replicated identifies biospecimens that have multiple observations per experiment.

```python
maeObj.replicated()
```
