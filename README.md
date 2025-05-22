# VQE algorithm performance analysis on various molecules

Repository for code implementation and exact results data for the bachelor thesis "VQE algorithm performance analysis on various molecules" at Vilnius University, Software Engineering program, 2025.

The paper is (will be) accessible on eLABa platform: [here (will be)](./)

Done: MAPPINGS/

Done: ANSATZE/hardware/

Done: examples.ipynb

# WIP

## Table of contents

- [Examples in the paper, not the experiment](./examples.ipynb)
- [Experiment. Mappings](./mappings.ipynb)
- [Experiment. Ansatze](./ansatze.ipynb)
- [Experiment. Optimizers](./optimizers.ipynb)
- [Tryout. On qhardware](./qhardware.ipynb)

## Manage python environments (notes)

1. Using Venv (without conda)

```bash
conda deactivate
python3 -m venv qiskit-env
source qiskit-env/bin/activate
pip install -r requirements.txt
deactivate
```

```bash
conda create --prefix ./tequila-env python=3.10
conda activate ./tequila-env
pip install tequila-basic # version 1.9.9, latest
conda install madtequila -c kottmann
```
