from collections import OrderedDict
from copy import deepcopy
from typing import Dict, MutableMapping, Optional, Sequence
from warnings import warn

from mudata import MuData
from pandas import DataFrame, concat
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment
from summarizedexperiment.type_checks import is_bioc_or_pandas_frame, is_list_of_type

from .types import SlicerArgTypes, SlicerResult, SlicerTypes, StrOrListStr

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"


class MultiAssayExperiment:
    """Container class for representing and managing multi-omics genomic experiments. Checkout the
    `R/MultiAssayExperiment <https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html>`_
    for more information.

    Attributes:
            experiments (MutableMapping[str, SummarizedExperiment]): A dictionary of
                experiments with experiment names as keys and the experiments as values.

                Each ``experiment`` may be either a
                :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`
                and any class that extends `SummarizedExperiment`.

            col_data (DataFrame]): Bio-specimen/sample information.
                The ``col_data`` may provide information about patients, cell lines, or
                other biological units.

                Each row in this table is an independent biological unit. Must contain an `index`
                that maps to primary in ``sample_map``.

            sample_map (DataFrame): Map biological units from
                :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`
                to the list of
                :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.experiments`.

                Must contain columns "assay", "primary" and "colname".

                - **assay** provides the names of the different experiments performed on the
                    biological units. All experiment names from
                    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.experiments`
                    must be present in this column.
                - **primary** contains the sample name. All names in this column must match with
                    row labels from
                    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`.
                - **colname** is the mapping of samples/cells within each experiment back to its
                    biosample information in
                    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`.

                Each sample in ``col_data`` may map to one or more columns per assay.

                This table can be created automatically in simple usecases, Checkout the
                :py:class:`~multiassayexperiment.io.interface.make_mae`, or import functions to
                read data as ``MultiAssayExperiment`` like
                :py:class:`~multiassayexperiment.io.mudata.from_mudata` and
                :py:class:`~multiassayexperiment.io.anndata.from_anndata`.

            metadata (MutableMapping, optional): Additional study level metadata. Defaults to None.
    """

    def __init__(
        self,
        experiments: MutableMapping[str, SummarizedExperiment],
        col_data: DataFrame,
        sample_map: DataFrame,
        metadata: Optional[MutableMapping] = None,
    ) -> None:
        """Construct an MAE."""
        self._validate_experiments(experiments)
        self._validate_sample_map(
            sample_map=sample_map, col_data=col_data, experiments=experiments
        )
        self._sample_map = sample_map
        self._col_data = col_data
        self._experiments = experiments

        self._metadata = metadata

    def _validate_experiments(
        self, experiments: MutableMapping[str, SummarizedExperiment]
    ):
        """Internal method to validate experiments.

        Raises:
            TypeError: If experiments is not a :py:class:`~dict`.
        """
        if not isinstance(experiments, dict):
            raise TypeError("experiments must be a dictionary.")

    def _validate_col_data(self, col_data: DataFrame):
        """Internal method to validate ``col_data``.

        Args:
            col_data (DataFrame): Column data.

                ``col_data`` may be either a :py:class:`~pandas.DataFrame` or
                :py:class:`~biocframe.BiocFrame.BiocFrame` object.

        Raises:
            TypeError: If object is not an expected type.
        """
        if not is_bioc_or_pandas_frame(col_data):
            raise TypeError(
                "`col_data` must be either a pandas DataFrame or a BiocFrame object."
            )

        if isinstance(col_data, DataFrame):
            if col_data.index is None:
                raise ValueError("`col_data` must have an index column.")
        else:
            if col_data.row_names is None:
                raise ValueError("`col_data` must have row names or labels.")

    def _validate_sample_map_with_col_data(
        self, sample_map: DataFrame, col_data: DataFrame
    ):
        """Internal method to validate ``sample_map`` and ``col_data``.

        Args:
            sample_map (DataFrame): Sample mapping.
            col_data (DataFrame): Column data.

        Raises:
            ValueError: If any of the checks fail.
        """
        # check if unique samples is same as in sample data
        _samples = list(sample_map["primary"])
        _sample_set = set(_samples)
        _sample_diff = _sample_set.difference(col_data.index.tolist())
        if len(_sample_diff) > 0:
            raise ValueError(
                "'primary' from `sample_map` has unknown samples not present in `col_data`."
            )

        if len(_sample_set) != col_data.shape[0]:
            warn("'primary' from `sample_map` & `col_data` has missing samples.")

    def _validate_sample_map_with_Expts(
        self,
        sample_map: DataFrame,
        experiments: MutableMapping[str, SummarizedExperiment],
    ):
        """Internal method to validate ``sample_map`` and ``experiments``.

        Args:
            sample_map (DataFrame): Sample mapping.
            experiments (MutableMapping[str, SummarizedExperiment]): Experiments.

        Raises:
            ValueError: If any of the checks fail.
        """
        # check if all assay names are in experiments
        smapUniqueAssaynames = set(sample_map["assay"])
        UniqueExperimentname = set(list(experiments.keys()))

        if (len(UniqueExperimentname) != len(smapUniqueAssaynames)) or (
            UniqueExperimentname != smapUniqueAssaynames
        ):
            raise ValueError(
                "'assays' from sample_map does not match with `experiments`."
            )

        # check if colnames exist
        agroups = sample_map.groupby(["assay"])
        for group, rows in agroups:
            if group not in experiments:
                raise ValueError(
                    f"Experiment '{group}' exists in `sample_map` but not in `experiments`."
                )

            gcol_data = experiments[group].col_data

            if set(rows["colname"].tolist()) != set(gcol_data.index.tolist()):
                raise ValueError(
                    f"Experiment '{group}' does not contain all columns in `sample_map`."
                )

    def _validate_sample_map(
        self,
        sample_map: DataFrame,
        col_data: DataFrame,
        experiments: MutableMapping[str, SummarizedExperiment],
    ):
        """Validate sample map.

        Args:
            sample_map (DataFrame): Sample map.
            col_data (DataFrame): Column data.
            experiments (MutableMapping[str, SummarizedExperiment]): Experiments.

        Raises:
            TypeError, ValueError: If any of the checks fail.
        """
        if not is_bioc_or_pandas_frame(sample_map):
            raise TypeError(
                "`sample_map` must be either a pandas `DataFrame` or a `BiocFrame` object."
            )

        if not set(["assay", "primary", "colname"]).issubset(list(sample_map.columns)):
            raise ValueError(
                "`sample_map` does not contain required columns: 'assay', 'primary/ and `'colname'."
            )

        self._validate_sample_map_with_col_data(sample_map, col_data)
        self._validate_sample_map_with_Expts(sample_map, experiments)

    def _validate(self):
        """Internal method to validate the object

        Raises:
            ValueError: If attributes don't match expectations.
        """
        self._validate_experiments(self._experiments)
        self._validate_col_data(self._col_data)
        self._validate_sample_map(self._sample_map, self._col_data, self._experiments)

    @property
    def experiments(
        self,
    ) -> Dict[str, SummarizedExperiment]:
        """Get experiments.

        Returns:
            Dict[str, SummarizedExperiment]: A dictionary of all experiments, with experiment
            names as keys and experiment data as value.
        """

        return self._experiments

    @experiments.setter
    def experiments(
        self,
        experiments: MutableMapping[str, SummarizedExperiment],
    ):
        """Set new experiments.

        Args:
            experiments (MutableMapping[str, SummarizedExperiment]): New experiments to set.
                A dictionary of experiments with experiment names as keys and the experiments
                as values.

                Each ``experiment`` may be either a
                :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`
                and any class that extends `SummarizedExperiment`.
        """

        self._validate_experiments(experiments)
        self._validate_sample_map_with_Expts(self._sample_map, experiments)
        self._experiments = experiments

    def experiment(
        self, name: str, with_sample_data: bool = False
    ) -> SummarizedExperiment:
        """Get experiment by name.

        If ``with_sample_data`` is True, a copy of the experiment object is returned.

        Args:
            name (str): Experiment name.
            with_sampleData (bool, optional): Whether to merge column data of the experiment with
                sample data from the MAE. Defaults to False.

        Raises:
            ValueError: If experiment name does not exist.

        Returns:
            SummarizedExperiment: A class that extends
            :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`.
        """
        if name not in self._experiments:
            raise ValueError(f"Experiment '{name}' does not exist.")

        expt = self.experiments[name]

        if with_sample_data is True:
            expt = deepcopy(expt)

            subset_map = self.sample_map[self.sample_map["assay"] == name]
            subset_map = subset_map.set_index("colname")

            expt_col_data = expt.col_data
            new_col_data = concat([subset_map, expt_col_data], axis=1)
            expt.col_data = new_col_data

        return expt

    @property
    def sample_map(self) -> DataFrame:
        """Get sample map between experiments and sample metadata.

        Returns:
            DataFrame: a DataFrame with sample mapping information.
        """
        return self._sample_map

    @sample_map.setter
    def sample_map(self, sample_map: DataFrame):
        """Set new sample mapping.

        Args:
            sample_map (DataFrame): New sample map.
        """
        self._validate_sample_map(sample_map, self._col_data, self._experiments)
        self._sample_map = sample_map

    @property
    def col_data(self) -> DataFrame:
        """Get sample metadata.

        Returns:
            DataFrame: Sample metadata.
        """
        return self._col_data

    @col_data.setter
    def col_data(self, col_data: DataFrame):
        """Set sample metadata.

        Args:
            col_data (DataFrame): New sample metadata.
        """
        self._validate_col_data(col_data)
        self._validate_sample_map_with_col_data(self._sample_map, col_data)
        self._col_data = col_data

    @property
    def assays(
        self,
    ) -> Dict[str, SummarizedExperiment]:
        """Get experiments.

        Alias to the
        :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.experiments`.

        Returns:
            Dict[str, SummarizedExperiment]: All experiments.
            A dictionary of experiments with experiment names as keys and the experiments
            as values.

            Each ``experiment`` may be either a
            :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`
            and any class that extends `SummarizedExperiment`.
        """
        return self.experiments

    @property
    def metadata(self) -> Optional[Dict]:
        """Get metadata.

        Returns:
            Optional[Dict]: Metadata if available.
        """
        return self._metadata

    @metadata.setter
    def metadata(self, metadata: MutableMapping):
        """Set metadata.

        Args:
            metadata (MutableMapping): New metadata object.
        """
        self._metadata = metadata

    def _subset_experiments(
        self, subset: StrOrListStr
    ) -> Dict[str, SummarizedExperiment]:
        """Internal method to subset experiments.

        Args:
            subset (StrOrListStr): May be an single experiment name to keep.
                Alternatively, ``subset`` may be a list of experiment names.

        Returns:
            Dict[str, SummarizedExperiment]: A dictionary with experiment names as keys
            and the subsetted experiment data as value.
        """
        if isinstance(subset, str):
            subset = [subset]

        if not is_list_of_type(subset, str):
            raise ValueError("All provided experiment names must be `strings`.")

        newExpt = OrderedDict()

        for texpt in subset:
            if texpt not in self.experiments:
                raise ValueError(
                    f"experiment {texpt} does not exist. should be {list(self.experiments.keys())}"
                )
            newExpt[texpt] = self.experiments[texpt]

        return newExpt

    def _slice(
        self,
        args: SlicerArgTypes,
    ) -> SlicerResult:
        """Internal method to slice by index.

        Args:
            args (SlicerArgTypes): Indices or names to slice. Tuple
                contains slices along dimensions (rows, columns, experiments).

                Each element in the tuple, might be either a integer vector (integer positions),
                boolean vector or :py:class:`~slice` object. Defaults to None.

        Raises:
            ValueError: Too many or too few slices.

        Returns:
            SlicerResult: Sliced row, cols and assays.
        """

        if len(args) == 0:
            raise ValueError("`args` must contain at least one slice.")

        rowIndices = args[0]
        colIndices = None
        exptIndices = None

        if len(args) > 1:
            colIndices = args[1]

        if len(args) > 2:
            exptIndices = args[2]

            if exptIndices is not None:
                if isinstance(exptIndices, str):
                    exptIndices = [exptIndices]

                if not all([isinstance(x, str) for x in exptIndices]):
                    raise ValueError(
                        "All `experiments` in the 3rd slice for `args` must be strings."
                    )

        if len(args) > 3:
            raise ValueError("`args` contains too many slices.")

        subsetExpts = self.experiments.copy()

        if exptIndices is not None:
            subsetExpts = self._subset_experiments(exptIndices)

        if rowIndices is not None:
            if isinstance(rowIndices, dict):
                incorrect = set(list(rowIndices.keys())).difference(
                    list(subsetExpts.keys())
                )
                if len(incorrect) > 0:
                    raise ValueError(
                        f"Incorrect experiment names provided: {', '.join(incorrect)}."
                    )

                for expname, expt in subsetExpts.items():
                    if expname in rowIndices:
                        subsetExpts[expname] = expt[rowIndices[expname], :]
                    else:
                        subsetExpts[expname] = expt

            elif isinstance(rowIndices, slice) or all(
                isinstance(x, int) for x in rowIndices
            ):
                for expname, expt in subsetExpts.items():
                    subsetExpts[expname] = expt[rowIndices, :]
            else:
                raise TypeError("Row indices is not an expected type.")

        if colIndices is not None:
            if isinstance(colIndices, dict):
                incorrect = set(list(colIndices.keys())).difference(
                    list(subsetExpts.keys())
                )
                if len(incorrect) > 0:
                    raise ValueError(
                        f"Incorrect experiment names provided: {', '.join(incorrect)}."
                    )

                for expname, expt in subsetExpts.items():
                    if expname in colIndices:
                        subsetExpts[expname] = expt[:, colIndices[expname]]
                    else:
                        subsetExpts[expname] = expt

            elif isinstance(colIndices, slice) or all(
                isinstance(x, int) for x in colIndices
            ):
                for expname, expt in subsetExpts.items():
                    subsetExpts[expname] = expt[:, colIndices]
            else:
                raise TypeError("Columns indices is not an expected type.")

        # filter sample_map
        subsetColnames = []
        subsetsample_map = DataFrame()
        for expname, expt in subsetExpts.items():
            subsetColnames.extend(expt.colnames)
            subsetsample_map = concat(
                [
                    subsetsample_map,
                    self.sample_map[
                        (self.sample_map["assay"] == expname)
                        & (self.sample_map["colname"].isin(expt.colnames))
                    ],
                ]
            )

        # filter col_data
        subsetcol_data = self.col_data[
            self.col_data.index.isin(subsetsample_map["primary"].unique().tolist())
        ]

        return SlicerResult(subsetExpts, subsetsample_map, subsetcol_data)

    def subset_by_experiments(self, subset: StrOrListStr) -> "MultiAssayExperiment":
        """Subset by experiment(s).

        Args:
            subset (StrOrListStr): May be an single experiment name to keep.
                Alternatively, ``subset`` may be a list of experiment names.

        Returns:
            MultiAssayExperiment: A new `MultiAssayExperiment` with the subset experiments.
        """
        sresult = self._slice(args=(None, None, subset))
        return MultiAssayExperiment(
            sresult.experiments, sresult.col_data, sresult.sample_map, self.metadata
        )

    def subset_by_row(self, subset: SlicerTypes) -> "MultiAssayExperiment":
        """Subset by rows.

        Args:
            subset (SlicerTypes): Row indices or names to slice.

                May be either a integer vector (integer positions),
                boolean vector or :py:class:`~slice` object. Defaults to None.

        Returns:
            MultiAssayExperiment: A new `MultiAssayExperiment` with the subset.
        """
        sresult = self._slice(args=(subset, None, None))
        return MultiAssayExperiment(
            sresult.experiments, sresult.col_data, sresult.sample_map, self.metadata
        )

    def subset_by_column(self, subset: SlicerTypes) -> "MultiAssayExperiment":
        """Subset by column.

        Args:
            subset (SlicerTypes): Column indices or names to slice.

                May be either a integer vector (integer positions),
                boolean vector or :py:class:`~slice` object. Defaults to None.

        Returns:
            MultiAssayExperiment: A new `MultiAssayExperiment` with the subset.
        """
        sresult = self._slice(args=(None, subset, None))
        return MultiAssayExperiment(
            sresult.experiments, sresult.col_data, sresult.sample_map, self.metadata
        )

    def __getitem__(self, args: SlicerArgTypes) -> "MultiAssayExperiment":
        """Subset a `MultiAssayExperiment`.

        Args:
            args (SlicerArgTypes): Indices or names to slice. Tuple
                contains slices along dimensions (rows, columns, experiments).

                Each element in the tuple, might be either a integer vector (integer positions),
                boolean vector or :py:class:`~slice` object. Defaults to None.

        Raises:
            ValueError: Too many or too few slices.

        Returns:
            MultiAssayExperiment: A new sliced `MultiAssayExperiment` object with the subsets.
        """
        sresult = self._slice(args=args)
        return MultiAssayExperiment(
            sresult.experiments, sresult.col_data, sresult.sample_map, self.metadata
        )

    def __str__(self) -> str:
        pattern = (
            f"Class MultiAssayExperiment with {len(self.experiments.keys())} experiments and {len(self.col_data)} samples \n"
            f"  experiments: "
        )

        for expname, expt in self.experiments.items():
            pattern = f"{pattern} \n    {expname}: {str(expt)}"
        return pattern

    def complete_cases(self) -> Sequence[bool]:
        """Identify samples that have data across all experiments.

        Returns:
            Sequence[bool]: A boolean vector, where each element is
            True if sample is present in all experiments or False.
        """
        vec = []
        for x in self.col_data.index.tolist():
            subset = self.sample_map[self.sample_map["primary"] == x]

            vec.append(len(subset["assay"].unique()) == len(self.experiments.keys()))

        return vec

    def replicated(self) -> Dict[str, Dict[str, Sequence[bool]]]:
        """Identify samples with replicates within each experiment.

        Returns:
            Dict[str, Dict[str, Sequence[bool]]]: A dictionary where experiment names
            are keys and values specify if the sample is replicated within each experiment.
        """
        replicates = {}
        allSamples = self.col_data.index.tolist()
        for expname, expt in self.experiments.items():
            if expname not in replicates:
                replicates[expname] = {}

                for s in allSamples:
                    replicates[expname][s] = []

            colnames = expt.colnames
            smap = self.sample_map[self.sample_map["assay"] == expname]

            for x in colnames:
                colmap = smap[smap["colname"] == x]
                for s in allSamples:
                    replicates[expname][s].append(s in colmap["primary"])

        return replicates

    def add_experiment(
        self,
        name: str,
        experiment: SummarizedExperiment,
        sample_map: DataFrame,
        col_data: Optional[DataFrame] = None,
    ):
        """Add a new experiment to `MultiAssayExperiment`.

        ``sample_map`` must be provided to map the cells or sample from this experiment back to
        :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`. This
        will be appended to the existing
        :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_map`.

        Optionally, ``col_data`` may be provided to add new sample information to the
        `MultiAssayExperiment`

        Args:
            name (str): Name of the new experiment.
            experiment (SummarizedExperiment): The experiment to add.
                Must extend
                :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment` class.

            sample_map (DataFrame): Sample map to append to the MAE.

            col_data (DataFrame, optional): Sample data to append to the MAE. Defaults to None.

        Raises:
            ValueError: If ``name`` is an existing experiment name in ``experiments``.
        """

        if name in self.experiments:
            raise ValueError(f"An experiment with {name} already exists.")

        self._validate_col_data(col_data)

        new_experiments = self.experiments.copy()
        new_experiments[name] = experiment

        new_col_data = col_data
        if new_col_data is not None:
            new_col_data = concat([self.col_data, col_data], axis=0)

        new_sample_map = concat([self.sample_map, sample_map], axis=0)

        self._validate_experiments(new_experiments)
        self._validate_sample_map(
            sample_map=new_sample_map,
            col_data=new_col_data,
            experiments=new_experiments,
        )

        self._experiments = new_experiments
        self._sample_map = new_sample_map
        self._col_data = new_col_data

    def to_mudata(self) -> MuData:
        """Transform `SingleCellExperiment` object to :py:class:`~mudata.MuData`.

        Returns:
            MuData: A `MuData` representation.
        """

        exptsList = OrderedDict()

        for expname, expt in self.experiments.items():
            if isinstance(expt, SingleCellExperiment):
                obj, adatas = expt.to_anndata(alts=True)

                exptsList[expname] = obj

                if adatas is not None:
                    for aname, aexpt in adatas.items():
                        exptsList[f"{expname}_{aname}"] = aexpt
            elif isinstance(expt, SummarizedExperiment):
                exptsList[expname] = expt.to_anndata()
            else:
                print(f"Experiment: '{expname}' is not supported!")

        return MuData(exptsList)
