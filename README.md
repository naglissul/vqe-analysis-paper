# ORIGINAL REPOSITORY - VQE algorithm performance analysis on various molecules

----

This is a very messy repository, since the thesis was mostly written in the last week before deadline. It's not changed because of academic fairness. 

Check out this repository for more order and cleaner code:

https://github.com/naglissul/vqe-analysis-paper-clean-code

----

Repository for code implementation and exact results data for the bachelor thesis "VQE algorithm performance analysis on various molecules" at Vilnius University, Software Engineering program, 2025.

The paper is accessible on Vilnius University online library: [here](https://virtualibiblioteka.vu.lt/discovery/fulldisplay?vid=370LABT_VU:VU&search_scope=VU_ETD_PDB&tab=VU_ETD_PDB&docid=alma9913200291908452&lang=en&context=L)

## Table of contents

- [Examples in the paper, not the experiment](./examples.ipynb)
- [Experiment. Mappings](./MAPPINGS/)
- [Experiment. Optimizers](./optimizers.ipynb)
- Experiment. Ansatze and full runs: [Data displayed in graphs](./ANSATZE/FINAL_ALL_FULL_RUNS.ipynb), Data collected for graphs:
  - [kUpCCGSD](./ANSATZE/kUpCCGSD-tequila-env.ipynb) - in the end of file
  - [ADAPT-VQE](./ANSATZE/ADAPT-VQE-simulator.ipynb)
  - [UCCSD](./ANSATZE/UCCSD-simulator.ipynb)
  - [Graph with quantum computer data](./ANSATZE/HEA.ipynb), [Data collected from quantum computer](./ANSATZE/hardware/run.ipynb)

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
