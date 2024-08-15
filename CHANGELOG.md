# Changelog

## Version 0.4.1 - 0.4.3

- Access an experiment by index.
- Helper methods to create sample mapping if not provided.
- Subset operations on samples.
- Update sphinx configuration to run snippets in the documentation.

## Version 0.4.0

This is a complete rewrite of the package, following the functional paradigm from our [developer notes](https://github.com/BiocPy/developer_guide#use-functional-discipline).

- `column_data` and `sample_map` are expected to be `BiocFrame` objects and will be converted if a pandas `DataFrame` is provided. This will allows us to reduce complexity and provide consistent downstream operations.
- A `sample_map` will be created by default if both `column_data` and `sample_map` are None.
- A warning is raised if `column_names` are empty for an experiment.
- A warning is raised if `column_data` contains duplicate row names.
- Streamlines subset operation.
- Reduce dependency on a number of external packages.
- Update docstrings, tests and docs.

## Version 0.3.0

This release migrates the package to a more palatable Google's Python style guide. A major modification to the package is with casing, all `camelCase` methods, functions and parameters are now `snake_case`.

In addition, docstrings and documentation has been updated to use sphinx's features of linking objects to their types. Sphinx now also documents private and special dunder methods (e.g. `__getitem__`, `__copy__` etc). Intersphinx has been updated to link to references from dependent packages.

In addition, pyscaffold has been updated to use "myst-parser" as the markdown compiler instead of rst.
As part of the pyscaffold setup, one may use pre-commits to run some of the routine tasks of linting and formatting before every commit. While this is sometimes annoying and can be ignored with `--no-verify`, it brings some consistency to the code base.

## Version 0.2.0

- Rewriting MAE class
  - add more methods especially slicing
- more robust tests
- updated documentation, tutorial

## Version 0.1

- Initialize MAE class
- Tests
- Documentation
