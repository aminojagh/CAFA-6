# CAFA-6

CAFA-6 is a Python project for building a complete protein function prediction pipeline around Gene Ontology annotation. The intended scope spans data creation from pinned public sources, preprocessing and validation, benchmark construction, model training, model evaluation, and deployment-oriented inference workflows.

## Project Goal

The project aims to provide an end-to-end system for protein function prediction that covers:

- reproducible dataset creation and preprocessing
- ontology-aware artifact validation and benchmark preparation
- training both baseline and new prediction models
- evaluating models with CAFA-style metrics and benchmark logic
- serving inference workflows for downstream use

## Current Scope

The repository currently includes the core foundations for the larger pipeline:

- shared data and configuration types
- source snapshot resolution utilities
- Swiss-Prot train taxonomy extraction from the pinned 2025_03 release archive
- Gene Ontology parsing, canonicalization, traversal, and propagation helpers
- deterministic readers and writers for core artifact formats
- validation utilities and automated unit tests
- a practical `main.ipynb` pipeline for the currently implemented stages

## Repository Layout

```text
src/cafa/                 Project source code
tests/                    Automated tests
main.ipynb                Primary project notebook
notebooks/                Experimental and exploratory notebooks
requirements.txt          Python dependencies
LICENSE                   Apache 2.0 license
```

## Setup

Create and activate the Conda environment, then install dependencies:

```bash
conda create -n CAFA python=3.11
conda activate CAFA
python3 -m pip install -r requirements.txt
```

For terminal-based usage, set the source path explicitly:

```bash
export PYTHONPATH=src
```

## Development Workflow

The codebase is being built incrementally. `main.ipynb` is the operational notebook for the stages implemented so far: it uses pinned real sources, recreates actual artifacts under `recreated_comp_data/`, and validates them against the reference data. The broader pipeline layers for train terms, train sequences, test extraction, benchmark assembly, training, evaluation, and inference will be added to that same flow as they are implemented.

## Testing

Run the test suite with:

```bash
conda activate CAFA
PYTHONPATH=src python3 -m unittest discover -s tests -p "test_*.py"
```

## License

This project is licensed under the Apache License 2.0. See `LICENSE`.
