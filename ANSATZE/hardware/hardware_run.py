# Running this as background process because of long running time

from lib import h2_hamiltonian_parity, make_hydrogen_hea
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_ibm_runtime import EstimatorV2 as Estimator
from qiskit import transpile
import time
from datetime import datetime
from qiskit_algorithms.optimizers import L_BFGS_B
import numpy as np

print()
print("=========================================")
print("H_2 case of 0.735 Angstrom distance with HEA ansatz from literature, Bravyi-Kitaev mapper, unlimiter iterations")
print("=========================================")
print()

backend_name = "ibm_brisbane"
service = QiskitRuntimeService()
backend = service.backend(backend_name)
print("Done getting backend")
print()
print("==============ANSATZ==================")
print()

ansatz = make_hydrogen_hea()
print(ansatz)

print()
print("============HAMILTONIAN==============")
print()

start = time.time()
hamiltonian = h2_hamiltonian_parity(distance=0.735)
end = time.time()
print(f"{0.735} case took time to make:", end - start)
print(hamiltonian)

print()
print("=============TRANSPILE====================")
print()

ansatz = transpile(
    [ansatz],
    optimization_level=1,
    backend=backend,
)[0]

hamiltonian = hamiltonian.apply_layout(ansatz.layout)

estimator = Estimator(mode=backend)
estimator.options.resilience_level = 1
estimator.options.dynamical_decoupling.enable = True
estimator.options.dynamical_decoupling.sequence_type = "XY4"

print("Transpiled ansatz:")
print(ansatz)
print("Transpiled hamiltonian:")
print(hamiltonian)

print()
print("==============RUN=====================")
print()

optimizer = L_BFGS_B()
initial_point = np.zeros(ansatz.num_parameters)
iter_energies=[]
q_runtimes=[]

def estimate_energy(parameters):
    start = time.time()
    job = estimator.run(
        [(ansatz, hamiltonian, parameters)]
    )
    print(f"{datetime.now()} Job submitted: {job.job_id()}")
    result = job.result()[0].data.evs
    end = time.time()
    print(f"Job time: {end - start}")
    q_runtimes.append(end - start)
    print(job.metrics())
    print(result)
    print(parameters)
    iter_energies.append(result)
    return result

print(datetime.now())
start = time.time()
result = optimizer.minimize(fun=estimate_energy, x0=initial_point)
end = time.time()
print(f"Total time: {end - start}")
print(result)
print("Final energy:", result.fun)
print(iter_energies)
print(q_runtimes)