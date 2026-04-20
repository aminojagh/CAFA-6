# CAFA

Reproducible CAFA-style protein function prediction project, built independently from Kaggle-provided competition files.

## Goal

- recreate CAFA-style public artifacts from pinned public sources
- validate recreated artifacts against `comp_data/` as reference only
- train a sequence-based GO prediction baseline
- evaluate models on a recreated time-split benchmark

## Structure

- `src/cafa/`: project source code
- `tests/`: automated tests for implemented modules
- `main.ipynb`: the project run interface

## Status

The project is being implemented incrementally, one coherent commit at a time.

Current implemented scope:

- shared domain types and `ProjectConfig`
- GO ontology parsing and traversal
- ontology tests
- notebook cells that load and exercise the GO graph
