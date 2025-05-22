# Running on IBM Quantum hardware

Prerequisites:

- Have IBM Quantum platform API key saved (check online tutorials, there is an official Qiskit tutorial for that)

## Run in background

```bash
: > output.txt && nohup python3 hardware_run.py >> output.txt 2>&1 &
```

Output is dynamically written in `output.txt` file.
