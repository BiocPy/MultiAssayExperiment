---
file_format: mystnb
kernelspec:
  name: python
---

# Multiple experiments

`MultiAssayExperiment` (MAE) simplifies the management of multiple experimental assays conducted on a shared set of specimens.

:::{note}
These classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.
:::

## Construction

An MAE contains three main entities,

- **Primary information** (`column_data`): Bio-specimen/sample information. The `column_data` may provide information about patients, cell lines, or other biological units. Each row in this table represents an independent biological unit. It must contain an `index` that maps to the 'primary' in `sample_map`.

- **Experiments** (`experiments`): Genomic data from each experiment. either a `SingleCellExperiment`, `SummarizedExperiment`, `RangedSummarizedExperiment` or any class that extends a `SummarizedExperiment`.

- **Sample Map** (`sample_map`): Map biological units from `column_data` to the list of `experiments`. Must contain columns,
    - **assay** provides the names of the different experiments performed on the biological units. All experiment names from experiments must be present in this column.
    - **primary** contains the sample name. All names in this column must match with row labels from col_data.
    - **colname** is the mapping of samples/cells within each experiment back to its biosample information in col_data.

    Each sample in ``column_data`` may map to one or more columns per assay.

Let's start by first creating few experiments:

```{code-cell}

from random import random

import numpy as np
from biocframe import BiocFrame
from genomicranges import GenomicRanges
from iranges import IRanges

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
    row_names=[f"sce_{i}" for i in range(6)],
)

col_data_se = BiocFrame({"treatment": ["ChIP", "Input"] * 3},
    row_names=[f"se_{i}" for i in range(6)],
)
```

More importantly, we need to provide `sample_map` information:

```{code-cell}
sample_map = BiocFrame({
    "assay": ["sce", "se"] * 6,
    "primary": ["sample1", "sample2"] * 6,
    "colname": ["sce_0", "se_0", "sce_1", "se_1", "sce_2", "se_2", "sce_3", "se_3", "sce_4", "se_4", "sce_5", "se_5"]
})

sample_data = BiocFrame({"samples": ["sample1", "sample2"]}, row_names= ["sample1", "sample2"])

print(sample_map)
```


Finally, we can create an `MultiAssayExperiment` object:

```{code-cell}
from multiassayexperiment import MultiAssayExperiment
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

mae = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    column_data=sample_data,
    sample_map=sample_map,
    metadata={"could be": "anything"},
)

print(mae)
```

### No sample mapping?

If both `column_data` and `sample_map` are `None`, the constructor naively creates sample mapping, with each `experiment` considered to be a independent `sample`. We add a sample to `column_data` in this pattern - ``unknown_sample_{experiment_name}``.

All cells from the each experiment are considered to be from the same sample and is reflected in `sample_map`.

:::{important}
***This is not a recommended approach, but if you don’t have sample mapping, then it doesn’t matter***.
:::

```{code-cell}
mae = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    metadata={"could be": "anything"},
)

print(mae)
```

### Interop with `anndata` or `mudata`

We provide convenient methods to easily convert a `MuData` object into an `MultiAssayExperiment`.

Let's create a mudata object:

```{code-cell}

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

from mudata import MuData
mdata = MuData({"rna": adata, "spatial": adata2})

print(mdata)
```

Lets convert this object to an `MAE`:

```{code-cell}
from multiassayexperiment import MultiAssayExperiment

mae_obj = MultiAssayExperiment.from_mudata(input=mdata)
print(mae_obj)
```


## Getters/Setters

Getters are available to access various attributes using either the property notation or functional style.

```{code-cell}
# access assays
print("experiment names (as property): ", mae.experiment_names)
print("experiment names (functional style): ", mae.get_experiment_names())

# access sample data
print(mae.column_data)
```

Check out the [class documentation](https://biocpy.github.io/MultiAssayExperiment/api/multiassayexperiment.html#multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment) for the full list of accessors and setters.


#### Row or column name accessors

A helper method is available to easily access row or column names across all experiments. This method returns a dictionary with experiment names as keys and the corresponding values, which can be either the row or column names depending on the function:

```{code-cell}
from rich import print as pprint
pprint("row names:", mae.get_row_names())
pprint("column names:", mae.get_column_names())
```

#### Access an experiment

One can access an experiment by name:

```{code-cell}
print(mae.experiment("se"))
```

Additionally you may access an experiment with the sample information included in the column data of the experiment:

:::{note}
This creates a copy of the experiment.
:::

```{code-cell}
expt_with_sample_info = mae.experiment("se", with_sample_data=True)
print(expt_with_sample_info)
```

:::{note}
For consistency with the R MAE's interface, we also provide `get_with_column_data` method, that performs the same operation.
:::

### Setters

::: {important}
All property-based setters are `in_place` operations, with further details discussed in [functional paradigm](../philosophy.qmd#functional-discipline) section.
:::

```{code-cell}
modified_column_data = mae.column_data.set_column("score", range(len(mae.column_data)))
modified_mae = mae.set_column_data(modified_column_data)
print(modified_mae)
```

Now, lets check the `column_data` on the original object.

```{code-cell}
print(mae.column_data)
```


## Subsetting

You can subset `MultiAssayExperiment` by using the subset (`[]`) operator. This operation accepts different slice input types, such as a boolean vector, a `slice` object, a list of indices, or names (if available) to subset.

`MultiAssayExperiment` allows subsetting by three dimensions: `rows`, `columns`, and `experiments`. ***`sample_map` is automatically filtered during this operation***.

### Subset by indices

```{code-cell}
subset_mae = mae[1:5, 0:4]
print(subset_mae)
```

### Subset by experiments dimension

The following creates a subset based on the experiments dimension:

```{code-cell}
subset_mae = mae[1:5, 0:1, ["se"]]
print(subset_mae)
```

:::{note}
If you're wondering about why the experiment "se" has 0 columns, it's important to note that our MAE implementation does not remove columns from an experiment solely because none of the columns map to the samples of interest. This approach aims to prevent unexpected outcomes in complex subset operations.
:::

## Helper functions

The `MultiAssayExperiment` class also provides a few methods for sample management.

### Complete cases

The `complete_cases` function is designed to identify samples that contain data across all experiments. It produces a boolean vector with the same length as the number of samples in `column_data`. Each element in the vector is `True` if the sample is present in all experiments, or `False` otherwise.

```{code-cell}
print(mae.complete_cases())
```

You can use this boolean vector to select samples with complete data across all assays or experiments.

```{code-cell}
subset_mae = mae[:, mae.complete_cases(),]
print(subset_mae)
```


### Replicates

This method identifies 'samples' with replicates within each experiment. The result is a dictionary where experiment names serve as keys, and the corresponding values indicate whether the sample is replicated within each experiment.


```{code-cell}
from rich import print as pprint # mainly for pretty printing
pprint(mae.replicated())
```


### Intersect rows

The `intersect_rows` finds common `row_names` across all experiments and returns a `MultiAssayExperiment` with those rows.

```{code-cell}
common_rows_mae = mae.intersect_rows()
print(common_rows_mae)
```

If you are only interested in finding common `row_names` across all experiments:

```{code-cell}
common_rows = mae.find_common_row_names()
print(common_rows)
```

### Empty MAE

While the necessity of an empty `MultiAssayExperiment` might not be apparent, for the sake of consistency with the rest of the tutorials:

```{code-cell}
mae = MultiAssayExperiment(experiments={})
print(mae)
```
