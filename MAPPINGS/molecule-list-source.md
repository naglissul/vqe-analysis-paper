# Molecule list from ChatGPT

> ChatGPT model 4o, 2025-05-22

    **PROMPT**: GIve me the list of molecules with increasing number of fermionic modes, make a dict with fermionic mode value next

Here’s a dict of small molecules and their **fermionic modes** (i.e., 2 × number of spatial orbitals in minimal STO-3G basis):

```python
molecule_modes = {
    "H2": 4,          # 2 spatial orbitals × 2 = 4 spin orbitals
    "HeH+": 4,
    "LiH": 12,
    "BeH2": 14,
    "CH4": 18,
    "NH3": 20,
    "H2O": 14,
    "HF": 12,
    "CO": 20,
    "N2": 20,
    "F2": 20,
    "CH3OH": 24,
    "C2H4": 28,
    "C2H6": 28,
    "C2H2": 28,
}
```

These values depend on the basis set used — this assumes **STO-3G** with **frozen core**, unless stated otherwise.

    **PROMPT**: molecules = {
        "H2": "H .0 .0 .0; H .0 .0 0.735",
        "LiH": "Li .0 .0 .0; H .0 .0 1.6",
        "BeH2": "Be .0 .0 .0; H .0 .0 1.3; H .0 .0 -1.3",
        # Add more...
    }


    Write  them in this format

```python
molecules = {
    "H2": ("H .0 .0 .0; H .0 .0 0.735", 4),
    "HeH+": ("He .0 .0 .0; H .0 .0 0.75", 4),
    "LiH": ("Li .0 .0 .0; H .0 .0 1.6", 12),
    "BeH2": ("Be .0 .0 .0; H .0 .0 1.3; H .0 .0 -1.3", 14),
    "HF": ("H .0 .0 .0; F .0 .0 0.917", 12),
    "H2O": ("O .0 .0 .0; H .0 .0 0.96; H .0 .76 0.48", 14),
    "NH3": ("N .0 .0 .0; H .0 .94 .0; H .0 -.47 .81; H .0 -.47 -.81", 20),
    "CH4": ("C .0 .0 .0; H .63 .63 .63; H -.63 -.63 .63; H -.63 .63 -.63; H .63 -.63 -.63", 18),
    "CO": ("C .0 .0 .0; O .0 .0 1.128", 20),
    "N2": ("N .0 .0 .0; N .0 .0 1.0977", 20),
    "F2": ("F .0 .0 .0; F .0 .0 1.4119", 20),
    "CH3OH": ("C .0 .0 .0; O .0 .0 1.43; H .0 0.94 -.63; H .0 -.94 -.63; H .94 .0 .63; H -.94 .0 .63; H .0 .0 2.46", 24),
    "C2H4": ("C .0 .0 .0; C .0 .0 1.34; H .0 .94 -0.6; H .0 -.94 -0.6; H .0 .94 1.94; H .0 -.94 1.94", 28),
    "C2H6": ("C .0 .0 .0; C .0 .0 1.54; H .94 .0 .0; H -.94 .0 .0; H .0 .94 .0; H .0 -.94 .0; H .0 .0 .94; H .0 .0 -0.94", 28),
    "C2H2": ("C .0 .0 .0; C .0 .0 1.20; H .0 .0 -1.1; H .0 .0 2.3", 28),
}
```
