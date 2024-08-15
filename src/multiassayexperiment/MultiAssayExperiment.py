from collections import OrderedDict, namedtuple
from typing import Any, Dict, List, Optional, Sequence, Union
from warnings import warn

import biocframe
import biocutils as ut
import summarizedexperiment as se

__author__ = "jkanche"
__copyright__ = "jkanche"
__license__ = "MIT"

SlicerResult = namedtuple("SlicerResult", ["experiments", "sample_map", "column_data"])


def _sanitize_frame(frame):
    if se._frameutils.is_pandas(frame):
        frame = biocframe.from_pandas(frame)

    return frame


def _validate_experiments(experiments):
    if not isinstance(experiments, dict):
        raise TypeError("experiments must be a dictionary.")

    for k, v in experiments.items():
        if not hasattr(v, "shape"):
            raise ValueError(f"experiment: {k} is not supported.")

        if v.column_names is None:
            warn(
                f"Experiment '{k}' does not contain column (cell/sample) names.",
                UserWarning,
            )


def _validate_column_data(column_data):
    if column_data is None:
        raise ValueError("'column_data' cannot be None.")

    if not isinstance(column_data, biocframe.BiocFrame):
        raise TypeError("'column_data' is not a `BiocFrame` object.")

    if column_data.row_names is None:
        raise ValueError("`column_data` must have row names or labels.")

    if len(set(column_data.row_names)) != len(column_data.row_names):
        warn("'column_data' has duplicate row_names.", UserWarning)


def _validate_sample_map_with_column_data(sample_map, column_data):
    # check if all samples are from primary exist in col data
    _samples = sample_map.get_column("primary")
    _sample_set = set(_samples)
    _sample_diff = _sample_set.difference(column_data.row_names)
    if len(_sample_diff) > 0:
        raise ValueError(
            "`sample_map`'s 'primary' contains samples not represented by 'row_names' from `column_data`."
        )

    if len(_sample_set) != column_data.shape[0]:
        warn("'primary' from `sample_map` & `column_data` mismatch.", UserWarning)


def _validate_sample_map_with_expts(sample_map, experiments):
    # check if all assay names are in experiments
    smap_unique_assays = set(sample_map.get_column("assay"))
    unique_expt_names = set(list(experiments.keys()))

    if (len(unique_expt_names) != len(smap_unique_assays)) or (
        unique_expt_names != smap_unique_assays
    ):
        warn(
            "'experiments' contains names not represented in 'sample_map' or vice-versa.",
            UserWarning,
        )

    # check if colnames exist
    agroups = sample_map.split("assay")
    for grp, rows in agroups.items():
        if grp not in experiments:
            warn(
                f"Experiment '{grp}' exists in `sample_map` but not in `experiments`.",
                UserWarning,
            )

        if set(rows.get_column("colname")) != set(experiments[grp].column_names):
            raise ValueError(
                f"Experiment '{grp}' does not contain all columns mentioned in `sample_map`."
            )


def _validate_sample_map(sample_map, column_data, experiments):
    if sample_map is None:
        raise ValueError("'sample_map' cannot be None.")

    if not isinstance(sample_map, biocframe.BiocFrame):
        raise TypeError("'sample_map' is not a `BiocFrame` object.")

    if not set(["assay", "primary", "colname"]).issubset(sample_map.column_names):
        raise ValueError(
            "'sample_map' does not contain required columns: 'assay', 'primary' and 'colname'."
        )

    _validate_column_data(column_data)
    _validate_sample_map_with_column_data(sample_map, column_data)
    _validate_sample_map_with_expts(sample_map, experiments)


def _create_smap_from_experiments(experiments):
    _all_assays = []
    _all_primary = []
    _all_colnames = []
    samples = []

    for expname, expt in experiments.items():
        colnames = expt.column_names
        asy_sample = f"unknown_sample_{expname}"
        _all_assays.extend([expname] * len(colnames))
        _all_primary.extend([asy_sample] * len(colnames))
        _all_colnames.extend(colnames)

        samples.append(asy_sample)

    sample_map = biocframe.BiocFrame(
        {"assay": _all_assays, "primary": _all_primary, "colname": _all_colnames}
    )
    col_data = biocframe.BiocFrame({"samples": samples}, row_names=samples)

    return col_data, sample_map


class MultiAssayExperiment:
    """Container class for representing and managing multi-omics genomic experiments.

    Checkout the
    `R/MultiAssayExperiment <https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html>`_
    for more information.
    """

    def __init__(
        self,
        experiments: Dict[str, Any],
        column_data: biocframe.BiocFrame = None,
        sample_map: biocframe.BiocFrame = None,
        metadata: Optional[dict] = None,
        validate: bool = True,
    ) -> None:
        """Initialize an instance of ``MultiAssayExperiment``.

        You may also initialize an ``MultiAssayExperiment`` using
        :py:class:`~multiassayexperiment.io.interface.make_mae` or by
        transform from :py:class:`~multiassayexperiment.io.mudata.from_mudata` and
        :py:class:`~multiassayexperiment.io.anndata.from_anndata` objects.

        If both ``column_data`` and ``sample_map`` are None, the constructor naively creates
        sample mapping, with each ``experiment`` considered to be a independent `sample`.
        We add a sample to :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`
        in this pattern - ``unknown_sample_{experiment_name}``. All cells from the same experiment are
        considered to be from the same sample and is reflected in
        :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_map`.

        Args:
            experiments:
                A dictionary containing experiments, with experiment names as keys and
                the experiments as values.

                Each ``experiment`` may be either a
                :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`
                or any class that extends ``SummarizedExperiment``.

            column_data:
                Bio-specimen/sample information.

                ``column_data`` may provide information about patients, cell lines, or other biological units.

                Each row in this table represents an independent biological unit. It must contain an `index`
                that maps to the 'primary' in
                :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_map`.

                Defaults to None.

            sample_map:
                Map biological units from
                :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.column_data`
                to the list of experiments.

                Must contain columns "assay", "primary", and "colname".

                - `assay` provides the names of the different experiments performed on the biological units.
                    All experiment names from
                    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.experiments` must
                    be present in this column.
                - `primary` contains the sample name. All names in this column must match with row labels from
                    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.column_data`.
                - `colname` is the mapping of samples/cells within each experiment back to its biosample information in
                    :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.column_data`.

                Each sample in ``column_data`` may map to one or more columns per assay.

                Defaults to None.

            metadata:
                Additional study-level metadata. Defaults to None.

            validate:
                Internal use only.
        """
        self._experiments = experiments if experiments is not None else {}
        self._metadata = metadata if metadata is not None else {}

        if sample_map is not None and column_data is not None:
            self._sample_map = _sanitize_frame(sample_map)
            self._column_data = _sanitize_frame(column_data)
        elif sample_map is None and column_data is None:
            # make a sample map
            self._column_data, self._sample_map = _create_smap_from_experiments(
                self._experiments
            )
        else:
            raise ValueError(
                "Either 'sample_map' or 'column_data' is `None`. Either both should be provided or set both to `None`."
            )

        if validate:
            _validate_experiments(self._experiments)
            _validate_column_data(self._column_data)
            _validate_sample_map(self._sample_map, self._column_data, self._experiments)

    def _define_output(self, in_place: bool = False) -> "MultiAssayExperiment":
        if in_place is True:
            return self
        else:
            return self.__copy__()

    #########################
    ######>> Copying <<######
    #########################

    def __deepcopy__(self, memo=None, _nil=[]):
        """
        Returns:
            A deep copy of the current ``MultiAssayExperiment``.
        """
        from copy import deepcopy

        _expts_copy = deepcopy(self._experiments)
        _sample_map_copy = deepcopy(self._sample_map)
        _column_data_copy = deepcopy(self._column_data)
        _metadata_copy = deepcopy(self.metadata)

        current_class_const = type(self)
        return current_class_const(
            experiments=_expts_copy,
            column_data=_column_data_copy,
            sample_map=_sample_map_copy,
            metadata=_metadata_copy,
        )

    def __copy__(self):
        """
        Returns:
            A shallow copy of the current ``MultiAssayExperiment``.
        """
        current_class_const = type(self)
        return current_class_const(
            experiments=self._experiments,
            column_data=self._column_data,
            sample_map=self._sample_map,
            metadata=self._metadata,
        )

    def copy(self):
        """Alias for :py:meth:`~__copy__`."""
        return self.__copy__()

    ##########################
    ######>> Printing <<######
    ##########################

    def __repr__(self) -> str:
        """
        Returns:
            A string representation.
        """
        output = f"{type(self).__name__}("
        output += ", experiments=" + ut.print_truncated_list(self._experiments)
        output += ", column_data=" + self._column_data.__repr__()
        output += ", sample_map=" + self._sample_map.__repr__()

        if len(self._metadata) > 0:
            output += ", metadata=" + ut.print_truncated_dict(self._metadata)

        output += ")"
        return output

    def __str__(self) -> str:
        """
        Returns:
            A pretty-printed string containing the contents of this object.
        """
        output = f"class: {type(self).__name__} containing {len(self.experiment_names)} experiments\n"

        for idx in range(len(self.experiment_names)):
            expt_name = self.experiment_names[idx]
            expt = self._experiments[expt_name]
            output += f"[{idx}] {expt_name}: {type(expt).__name__} with {expt.shape[0]} rows and {expt.shape[1]} columns \n"  # noqa

        output += f"column_data columns({len(self._column_data.column_names)}): "
        output += f"{ut.print_truncated_list(self._column_data.column_names)}\n"

        output += f"sample_map columns({len(self._sample_map.column_names)}): "
        output += f"{ut.print_truncated_list(self._sample_map.column_names)}\n"

        output += f"metadata({str(len(self.metadata))}): {ut.print_truncated_list(list(self.metadata.keys()), sep=' ', include_brackets=False, transform=lambda y: y)}\n"  # noqa

        return output

    #############################
    ######>> experiments <<######
    #############################

    def get_experiments(self) -> Dict[str, Any]:
        """Access experiments.

        Returns:
            A dictionary of all experiments, with experiment
            names as keys and experiment data as value.
        """

        return self._experiments

    def set_experiments(
        self, experiments: Dict[str, Any], in_place: bool = False
    ) -> "MultiAssayExperiment":
        """Set new experiments.

        Args:
            experiments:
                New experiments to set. A dictionary of experiments with experiment names as keys and
                the experiments as values.

                Each ``experiment`` may be either a
                :py:class:`~summarizedexperiment.SummarizedExperiment.SummarizedExperiment`
                or any class that extends ``SummarizedExperiment``.

            in_place:
                Whether to modify the ``MultiAssayExperiment`` in place.

        Returns:
            A modified ``MultiAssayExperiment`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        _validate_experiments(experiments)
        _validate_sample_map_with_expts(self._sample_map, experiments)

        output = self._define_output(in_place)
        output._experiments = experiments
        return output

    @property
    def experiments(
        self,
    ) -> Dict[str, Any]:
        """Alias for :py:meth:`~get_experiments`."""
        return self.get_experiments()

    @experiments.setter
    def experiments(
        self,
        experiments: Dict[str, Any],
    ):
        """Alias for :py:meth:`~set_experiments` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'experiments' is an in-place operation, use 'set_experiments' instead",
            UserWarning,
        )
        self.set_experiments(experiments, in_place=True)

    @property
    def assays(self) -> Dict[str, Any]:
        """Alias for :py:meth:`~get_experiments`."""
        return self.get_experiments()

    ##################################
    ######>> experiment names <<######
    ##################################

    def get_experiment_names(self) -> List[str]:
        """Get experiment names.

        Returns:
            List of experiment names.
        """
        return list(self._experiments.keys())

    def set_experiment_names(
        self, names: List[str], in_place: bool = False
    ) -> "MultiAssayExperiment":
        """Replace :py:attr:`~experiments`'s names.

        Args:
            names:
                New names.

            in_place:
                Whether to modify the ``MultiAssayExperiment`` in place.

        Returns:
            A modified ``MultiAssayExperiment`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        current_names = self.get_experiment_names()
        if len(names) != len(current_names):
            raise ValueError(
                "Length of 'names' does not match the number of `experiments`."
            )

        new_experiments = OrderedDict()
        for idx in range(len(names)):
            new_experiments[names[idx]] = self._experiments.pop(current_names[idx])

        output = self._define_output(in_place)
        output._experiments = new_experiments
        return output

    @property
    def experiment_names(self) -> List[str]:
        """Alias for :py:meth:`~get_experiment_names`."""
        return self.get_experiment_names()

    @experiment_names.setter
    def experiment_names(self, names: List[str]):
        """Alias for :py:meth:`~set_experiment_names` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'experiment_names' is an in-place operation, use 'set_experiment_names' instead",
            UserWarning,
        )
        self.set_experiment_names(names, in_place=True)

    #####################################
    ######>> experiment accessor <<######
    #####################################

    def experiment(self, name: Union[int, str], with_sample_data: bool = False) -> Any:
        """Get an experiment by name.

        Args:
            name:
                Name or index position of the experiment.

            with_sample_data:
                Whether to merge column data of the experiment with
                :py:attr:`~sample_data` from the MAE.

                Defaults to False.

        Raises:
            AttributeError:
                If the experiment name does not exist.
            IndexError:
                If index is greater than the number of experiments.

        Returns:
            The experiment object.

            If ``with_sample_data`` is `True`, a copy of the experiment object is returned.
        """
        _name = name
        if isinstance(name, int):
            if name < 0:
                raise IndexError("Index cannot be negative.")

            if name > len(self.experiment_names):
                raise IndexError("Index greater than the number of assays.")

            _name = self.experiment_names[name]
            expt = self._experiments[_name]
        elif isinstance(name, str):
            if name not in self._experiments:
                raise ValueError(f"'{name}' is not a valid experiment name.")

            expt = self.experiments[name]
        else:
            raise TypeError(
                f"'experiment' must be a string or integer, provided '{type(name)}'."
            )

        if with_sample_data is True:
            assay_splits = self.sample_map.split("assay", only_indices=True)
            subset_map = self.sample_map[assay_splits[_name],]
            subset_map = subset_map.set_row_names(subset_map.get_column("colname"))

            expt_column_data = expt.column_data
            new_column_data = biocframe.merge(
                [subset_map, expt_column_data], join="outer"
            )

            new_column_data = biocframe.merge(
                [new_column_data, self._column_data], join="left"
            )

            return expt.set_column_data(new_column_data, in_place=False)

        return expt

    def get_with_column_data(self, name: str) -> Any:
        """Alias to :py:meth:`~experiment`.

        Consistency with Bioconductor's naming of the same function.
        """
        return self.experiment(name, with_sample_data=True)

    ############################
    ######>> sample map <<######
    ############################

    def get_sample_map(self) -> biocframe.BiocFrame:
        """Acess sample map.

        Returns:
            A :py:class:`~biocframe.BiocFrame.BiocFrame` with sample mapping information.
        """
        return self._sample_map

    def set_sample_map(
        self, sample_map: biocframe.BiocFrame, in_place: bool = False
    ) -> "MultiAssayExperiment":
        """Set new sample mapping.

        Args:
            sample_map:
                New sample map.

            in_place:
                Whether to modify the ``MultiAssayExperiment`` in place.

        Returns:
            A modified ``MultiAssayExperiment`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        sample_map = _sanitize_frame(sample_map)
        _validate_sample_map(sample_map, self._column_data, self._experiments)

        output = self._define_output(in_place)
        output._sample_map = sample_map
        return output

    @property
    def sample_map(self) -> biocframe.BiocFrame:
        """Alias for :py:meth:`~get_sample_map`."""
        return self.get_sample_map()

    @sample_map.setter
    def sample_map(self, sample_map: biocframe.BiocFrame):
        """Alias for :py:meth:`~set_sample_map` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'sample_map' is an in-place operation, use 'set_sample_map' instead",
            UserWarning,
        )
        self.set_sample_map(sample_map, in_place=True)

    #############################
    ######>> column_data <<######
    #############################

    def get_column_data(self) -> biocframe.BiocFrame:
        """Get sample metadata.

        Returns:
            A :py:class:`~biocframe.BiocFrame.BiocFrame` containing sample metadata.
        """
        return self._column_data

    def set_column_data(
        self, column_data: biocframe.BiocFrame, in_place: bool = False
    ) -> "MultiAssayExperiment":
        """Set new sample metadata.

        Args:
            column_data:
                New sample metadata.

            in_place:
                Whether to modify the ``MultiAssayExperiment`` in place.

        Returns:
            A modified ``MultiAssayExperiment`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        column_data = _sanitize_frame(column_data)

        _validate_column_data(column_data)
        _validate_sample_map_with_column_data(self._sample_map, column_data)

        output = self._define_output(in_place)
        output._column_data = column_data
        return output

    @property
    def column_data(self) -> biocframe.BiocFrame:
        """Alias for :py:meth:`~get_column_data`."""
        return self.get_column_data()

    @column_data.setter
    def column_data(self, column_data: biocframe.BiocFrame):
        """Alias for :py:meth:`~set_column_data` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'column_data' is an in-place operation, use 'set_column_data' instead",
            UserWarning,
        )
        self.set_column_data(column_data, in_place=True)

    ###########################
    ######>> metadata <<#######
    ###########################

    def get_metadata(self) -> dict:
        """
        Returns:
            Dictionary of metadata for this object.
        """
        return self._metadata

    def set_metadata(
        self, metadata: dict, in_place: bool = False
    ) -> "MultiAssayExperiment":
        """Set additional metadata.

        Args:
            metadata:
                New metadata for this object.

            in_place:
                Whether to modify the ``MultiAssayExperiment`` in place.

        Returns:
            A modified ``MultiAssayExperiment`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """
        if not isinstance(metadata, dict):
            raise TypeError(
                f"`metadata` must be a dictionary, provided {type(metadata)}."
            )
        output = self._define_output(in_place)
        output._metadata = metadata
        return output

    @property
    def metadata(self) -> dict:
        """Alias for :py:attr:`~get_metadata`."""
        return self.get_metadata()

    @metadata.setter
    def metadata(self, metadata: dict):
        """Alias for :py:attr:`~set_metadata` with ``in_place = True``.

        As this mutates the original object, a warning is raised.
        """
        warn(
            "Setting property 'metadata' is an in-place operation, use 'set_metadata' instead",
            UserWarning,
        )
        self.set_metadata(metadata, in_place=True)

    #########################
    ######>> subset <<#######
    #########################

    def _normalize_column_slice(self, columns: Union[str, int, bool, Sequence, slice]):
        _scalar = None
        if columns != slice(None):
            columns, _scalar = ut.normalize_subscript(
                columns, len(self._column_data), self._column_data.row_names
            )

        return columns, _scalar

    def _filter_sample_map(self, columns: Union[str, int, bool, Sequence, slice]):
        _samples_to_filter = self._column_data[columns,].row_names

        column_names_to_keep = {}
        for i in self.experiment_names:
            column_names_to_keep[i] = []

        for _, row in self._sample_map:
            if row["primary"] in _samples_to_filter:
                column_names_to_keep[row["assay"]].append(row["colname"])

        return column_names_to_keep

    def subset_experiments(
        self,
        rows: Optional[Union[str, int, bool, Sequence]],
        columns: Optional[Union[str, int, bool, Sequence]],
        experiment_names: Union[str, int, bool, Sequence],
    ) -> Dict[str, Any]:
        """Subset experiments.

        Args:
            rows:
                Row indices to subset.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

            columns:
                Column indices (from :py:attr:`~column_data`) to subset.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

            experiment_names:
                Experiment names to keep.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

                Check :py:attr:`~experiment_names` for a list of valid experiment names.

        Returns:
            A dictionary with experiment names as keys
            and the subsetted experiment data as value.
        """
        _expts_copy = self._experiments.copy()
        if experiment_names is None:
            experiment_names = slice(None)

        if experiment_names != slice(None):
            expts, _ = ut.normalize_subscript(
                experiment_names, len(self.experiment_names), self.experiment_names
            )

            to_keep = [self.experiment_names[idx] for idx in expts]

            new_expt = OrderedDict()
            for texpt in to_keep:
                new_expt[texpt] = _expts_copy[texpt]

            _expts_copy = new_expt

        if rows != slice(None):
            for k, v in _expts_copy.items():
                _expts_copy[k] = v[rows,]

        columns, _ = self._normalize_column_slice(columns)
        if columns != slice(None):
            _col_dict = self._filter_sample_map(columns)

            for k, v in _expts_copy.items():
                if k in _col_dict:
                    if len(_col_dict[k]) != 0:
                        _matched_indices = ut.match(_col_dict[k], v.column_names)
                    else:
                        _matched_indices = []
                    _expts_copy[k] = v[:, list(_matched_indices)]

        return _expts_copy

    def _generic_slice(
        self,
        rows: Optional[Union[str, int, bool, Sequence]] = None,
        columns: Optional[Union[str, int, bool, Sequence]] = None,
        experiments: Optional[Union[str, int, bool, Sequence]] = None,
    ) -> SlicerResult:
        """Slice ``MultiAssayExperiment`` along the rows and/or columns, based on their indices or names.

        Args:
            rows:
                Rows to be extracted.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

            columns:
                Columns to be extracted.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

            experiment:
                Experiments to extract.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

                Check :py:attr:`~experiment_names` for a list of valid experiment names.

        Returns:
            The sliced tuple containing the new sample_map, column_data and experiments
            for use in downstream methods.
        """

        if rows is None:
            rows = slice(None)

        if columns is None:
            columns = slice(None)
        columns, _ = self._normalize_column_slice(columns)

        # filter column_data
        _new_column_data = self._column_data[columns,]

        if experiments is None:
            experiments = slice(None)

        _new_experiments = self.subset_experiments(
            experiment_names=experiments, rows=rows, columns=columns
        )

        # filter sample_map
        smap_indices_to_keep = []
        for expname, expt in _new_experiments.items():
            counter = 0
            for _, row in self._sample_map:
                if (
                    row["assay"] == expname
                    and row["primary"] in _new_column_data.row_names
                    and row["colname"] in expt.column_names
                ):
                    smap_indices_to_keep.append(counter)
                counter += 1
        _new_sample_map = self._sample_map[list(set(smap_indices_to_keep)),]

        return SlicerResult(_new_experiments, _new_sample_map, _new_column_data)

    def subset_by_experiments(
        self, experiments: Union[str, int, bool, Sequence]
    ) -> "MultiAssayExperiment":
        """Subset by experiment(s).

        Args:
            experiments:
                Experiments to extract.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

                Check :py:attr:`~experiment_names` for a list of valid experiment names.

        Returns:
            A new `MultiAssayExperiment` with the subset experiments.
        """
        sresult = self._generic_slice(experiments=experiments)
        return MultiAssayExperiment(
            sresult.experiments, sresult.column_data, sresult.sample_map, self.metadata
        )

    def subset_by_row(
        self, rows: Union[str, int, bool, Sequence]
    ) -> "MultiAssayExperiment":
        """Subset by rows.

        Args:
            rows:
                Rows to be extracted.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

        Returns:
            A new `MultiAssayExperiment` with the subsetted rows.
        """
        sresult = self._generic_slice(rows=rows)
        return MultiAssayExperiment(
            sresult.experiments, sresult.column_data, sresult.sample_map, self.metadata
        )

    def subset_by_column(
        self, columns: Union[str, int, bool, Sequence]
    ) -> "MultiAssayExperiment":
        """Subset by column.

        Args:
            columns:
                Columns to be extracted.

                Integer indices, a boolean filter, or (if the current object is
                named) names specifying the ranges to be extracted, see
                :py:meth:`~biocutils.normalize_subscript.normalize_subscript`.

        Returns:
            A new `MultiAssayExperiment` with the subsetted columns.
        """
        sresult = self._generic_slice(columns=columns)
        return MultiAssayExperiment(
            sresult.experiments, sresult.column_data, sresult.sample_map, self.metadata
        )

    def __getitem__(self, args: tuple) -> "MultiAssayExperiment":
        """Subset a `MultiAssayExperiment`.

        Args:
            args:
                Tuple containing slices along dimensions (rows, columns, experiments).

                Each element in the tuple, might be either a integer vector (integer positions),
                boolean vector or :py:class:`~slice` object. Defaults to None.

        Raises:
            ValueError:
                Too many or too few slices.

        Returns:
            A new sliced `MultiAssayExperiment` object with the subsets.
        """

        if isinstance(args, tuple):
            if len(args) == 0:
                raise ValueError("At least one slice argument must be provided.")

            if len(args) == 1:
                sresult = self._generic_slice(rows=args[0])
                return MultiAssayExperiment(
                    sresult.experiments,
                    sresult.column_data,
                    sresult.sample_map,
                    self.metadata,
                )
            elif len(args) == 2:
                sresult = self._generic_slice(rows=args[0], columns=args[1])
                return MultiAssayExperiment(
                    sresult.experiments,
                    sresult.column_data,
                    sresult.sample_map,
                    self.metadata,
                )
            elif len(args) == 3:
                sresult = self._generic_slice(
                    rows=args[0], columns=args[1], experiments=args[2]
                )
                return MultiAssayExperiment(
                    sresult.experiments,
                    sresult.column_data,
                    sresult.sample_map,
                    self.metadata,
                )
            else:
                raise ValueError(
                    f"`{type(self).__name__}` only supports 3-dimensional slicing along rows, columns and/or experiments."  # noqa
                )

        raise TypeError("'args' must be a tuple")

    ################################
    ######>> miscellaneous <<#######
    ################################

    def complete_cases(self) -> Sequence[bool]:
        """Identify samples that have data across all experiments.

        Returns:
            A boolean vector same as the number of samples in 'column_data',
            where each element is True if sample is present in all experiments or False.
        """
        vec = []
        for x in self._column_data.row_names:
            _primary = self._sample_map.get_column("primary")

            smap_indices_to_keep = []
            for rdx in range(len(_primary)):
                if _primary[rdx] == x:
                    smap_indices_to_keep.append(rdx)

            subset = self.sample_map[list(set(smap_indices_to_keep)),]
            vec.append(set(subset.get_column("assay")) == set(self.experiment_names))

        return vec

    def replicated(self) -> Dict[str, Dict[str, Sequence[bool]]]:
        """Identify samples with replicates within each experiment.

        Returns:
            A dictionary where experiment names
            are keys and values specify if the sample is replicated within each experiment.
        """
        replicates = {}
        all_samples = self._column_data.row_names
        for expname, expt in self._experiments.items():
            if expname not in replicates:
                replicates[expname] = {}

                for s in all_samples:
                    replicates[expname][s] = [False] * expt.shape[1]

            smap_indices_to_keep = []

            _assay = self._sample_map.get_column("assay")
            for adx in range(len(_assay)):
                if _assay[adx] == expname:
                    smap_indices_to_keep.append(adx)

            subset_smap = self.sample_map[list(set(smap_indices_to_keep)),]

            counter = 0
            for _, row in subset_smap:
                if row["assay"] == expname:
                    replicates[expname][row["primary"]][counter] = True
                counter += 1

        return replicates

    def find_common_row_names(self) -> List[str]:
        """Finds common row names across all experiments."""

        _common = None
        for _, expt in self._experiments.items():
            if _common is None:
                _common = set(expt.row_names)
            else:
                _common = set(expt.row_names).intersection(_common)

        return _common

    def intersect_rows(self) -> "MultiAssayExperiment":
        """Finds common row names across all experiments and filters the MAE to these rows.

        Returns:
            A new sliced `MultiAssayExperiment` object with the filtered rows.
        """
        _common_row_names = self.find_common_row_names()
        sresult = self._generic_slice(rows=_common_row_names)
        return MultiAssayExperiment(
            sresult.experiments,
            sresult.column_data,
            sresult.sample_map,
            self.metadata,
        )

    ####################################
    ######>> row or column names <<#####
    ####################################

    def get_row_names(self) -> Dict[str, Optional[ut.Names]]:
        """
        Returns:
            Dictionary, with experiment names as keys, and row names as values.
        """
        _all_row_names = {}
        for expname, expt in self._experiments.items():
            _all_row_names[expname] = expt.row_names

        return _all_row_names

    @property
    def rownames(self) -> Dict[str, Optional[ut.Names]]:
        """Alias for :py:attr:`~get_row_names`, provided for back-compatibility."""
        return self.get_row_names()

    def get_column_names(self) -> Dict[str, Optional[ut.Names]]:
        """
        Returns:
            Dictionary, with experiment names as keys, and the column names as values.
        """
        _all_row_names = {}
        for expname, expt in self._experiments.items():
            _all_row_names[expname] = expt.column_names

        return _all_row_names

    @property
    def columnnames(self) -> Dict[str, Optional[ut.Names]]:
        """Alias for :py:attr:`~get_column_names`, provided for back-compatibility."""
        return self.get_column_names()

    @property
    def colnames(self) -> Dict[str, Optional[ut.Names]]:
        """Alias for :py:attr:`~get_column_names`, provided for back-compatibility."""
        return self.get_column_names()

    @property
    def column_names(self) -> Dict[str, Optional[ut.Names]]:
        """Alias for :py:attr:`~get_column_names`, provided for back-compatibility."""
        return self.get_column_names()

    #################################
    ######>> add experiment <<#######
    #################################

    def add_experiment(
        self,
        name: str,
        experiment: Any,
        sample_map: biocframe.BiocFrame,
        column_data: Optional[biocframe.BiocFrame] = None,
        in_place: bool = False,
    ) -> "MultiAssayExperiment":
        """Add a new experiment to `MultiAssayExperiment`.

        ``sample_map`` must be provided to map the columns from this experiment to
        :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.column_data`.
        This will be appended to the existing
        :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_map`.

        Optionally, ``column_data`` may be provided to add new sample information.

        Args:
            name:
                Name of the new experiment.

            experiment:
                The experiment to add.

            sample_map:
                Sample map to append to the MAE.

            column_data:
                Sample data to append to the MAE.

                Defaults to None.

            in_place:
                Whether to modify the ``MultiAssayExperiment`` in place.
                Defaults to False.

        Returns:
            A modified ``MultiAssayExperiment`` object, either as a copy of the original
            or as a reference to the (in-place-modified) original.
        """

        if name in self.experiments:
            raise ValueError(f"An experiment with {name} already exists.")

        _new_column_data = self._column_data
        if column_data is not None:
            column_data = _sanitize_frame(column_data)
            _new_column_data = ut.combine_rows(self._column_data, column_data)

        sample_map = _sanitize_frame(sample_map)
        _new_sample_map = ut.combine_rows(self._sample_map, sample_map)

        _new_experiments = self._experiments.copy()
        _new_experiments[name] = experiment

        _validate_column_data(_new_column_data)
        _validate_experiments(_new_experiments)
        _validate_sample_map(
            sample_map=_new_sample_map,
            column_data=_new_column_data,
            experiments=_new_experiments,
        )

        output = self._define_output(in_place)
        output._experiments = _new_experiments
        output._sample_map = _new_sample_map
        output._column_data = _new_column_data

        return output

    #################################
    ######>> mudata interop <<#######
    #################################

    def to_mudata(self):
        """Transform ``MultiAssayExperiment`` object to :py:class:`~mudata.MuData`.

        Returns:
            A `MuData` representation.
        """
        from mudata import MuData
        from singlecellexperiment import SingleCellExperiment

        exptsList = OrderedDict()

        for expname, expt in self._experiments.items():
            if isinstance(expt, SingleCellExperiment):
                obj, adatas = expt.to_anndata(include_alternative_experiments=True)

                exptsList[expname] = obj

                if adatas is not None:
                    for aname, aexpt in adatas.items():
                        exptsList[f"{expname}_{aname}"] = aexpt
            elif isinstance(expt, se.SummarizedExperiment):
                exptsList[expname] = expt.to_anndata()
            else:
                print(f"Experiment: '{expname}' is not supported!")

        return MuData(exptsList)

    @classmethod
    def from_mudata(cls, input: "mudata.MuData") -> "MultiAssayExperiment":
        """Create a ``MultiAssayExperiment`` object from :py:class:`~mudata.MuData`.

        The import naively creates sample mapping, each ``experiment`` is considered to be a `sample`.
        We add a sample with the following pattern - ``"unknown_sample_{experiment_name}"`` to
        :py:attr:`~col_data`. All cells from the same experiment are considered to be extracted from
        the same sample and is reflected in :py:attr:`~sample_map`.

        Args:
            input:
                MuData object.

        Raises:
            Exception:
                If ``mudata`` object is read in backed mode :py:attr:`~mudata.MuData.isbacked`.

        Returns:
            ``MultiAssayExperiment`` object.
        """

        from singlecellexperiment import SingleCellExperiment

        if input.isbacked is True:
            raise Exception("backed mode is currently not supported.")

        experiments = OrderedDict()

        _all_assays = []
        _all_primary = []
        _all_colnames = []
        samples = []

        for asy, adata in input.mod.items():
            experiments[asy] = SingleCellExperiment.from_anndata(adata)

            colnames = None
            if adata.obs.index.tolist() is not None:
                colnames = adata.obs.index.tolist()
            else:
                colnames = range(len(adata.shape[0]))

            asy_sample = f"unknown_sample_{asy}"

            _all_assays.extend([asy] * len(colnames))
            _all_primary.extend([asy_sample] * len(colnames))
            _all_colnames.extend(colnames)

            samples.append(asy_sample)

        sample_map = biocframe.BiocFrame(
            {"assay": _all_assays, "primary": _all_primary, "colname": _all_colnames}
        )
        col_data = biocframe.BiocFrame({"samples": samples}, row_names=samples)

        return cls(
            experiments=experiments,
            column_data=col_data,
            sample_map=sample_map,
            metadata=input.uns,
        )

    @classmethod
    def from_anndata(
        cls, input: "anndata.AnnData", name: str = "unknown"
    ) -> "MultiAssayExperiment":
        """Create a ``MultiAssayExperiment`` from :py:class:`~anndata.AnnData`.

        Since :py:class:`~anndata.AnnData` does not contain sample information,
        sample named ``"unknown_sample"`` will be added to
        :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.col_data`.
        All cells are considered to be extracted from this sample and is reflected in
        :py:attr:`~multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment.sample_map`.

        Args:
            input:
                An ``AnnData`` object.

            name:
                Name for the experiment.

                Defaults to "unknown".

        Returns:
            An ``MultiAssayExperiment``.
        """
        from singlecellexperiment import SingleCellExperiment

        scexpt = SingleCellExperiment.from_anndata(input=input)

        experiments = {name: scexpt}

        col_data = biocframe.BiocFrame(
            {"samples": ["unknown_sample"]}, row_names=["unknown_sample"]
        )

        colnames = None

        if input.obs.index.tolist() is not None:
            colnames = input.obs.index.tolist()
        else:
            colnames = range(len(input.shape[0]))

        sample_map = biocframe.BiocFrame(
            {
                "colname": colnames,
                "assay": ["unknown"] * len(colnames),
                "primary": ["unknown_sample"] * len(colnames),
            }
        )

        return cls(
            experiments=experiments,
            column_data=col_data,
            sample_map=sample_map,
            metadata=input.uns,
        )


@ut.extract_row_names.register(MultiAssayExperiment)
def _rownames_mae(x: MultiAssayExperiment):
    return x.get_row_names()


@ut.extract_column_names.register(MultiAssayExperiment)
def _colnames_mae(x: MultiAssayExperiment):
    return x.get_column_names()
