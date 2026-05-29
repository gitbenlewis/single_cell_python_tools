# Development

## Repository Layout

- [`pyproject.toml`](../pyproject.toml): package metadata and install requirements.
- [`sctl.yaml`](../sctl.yaml): conda environment for the Scanpy stack.
- [`src/single_cell_python_tools`](../src/single_cell_python_tools): package source.
- [`Example_notebooks`](../Example_notebooks): notebook examples and generated example outputs.
- [`docs`](.): canonical Markdown docs for the GitHub Pages docs hub.

## Documentation Maintenance

The docs are static Markdown. When public functions, signatures, or exported
names change, update the relevant docs page and [API Reference](api-reference.md)
in the same change.

The docs hub stages these files into the Sphinx site at build time. Do not copy
large notebook outputs or example data into the docs directory.

## Public API Notes

Current public names include several legacy spellings. Preserve them in docs and
examples unless the package code changes intentionally:

- `basic_filitering`
- `filter_cells_by_anotated_QC_gene`
- `_plots_depreciated`

