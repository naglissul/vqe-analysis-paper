from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.mappers import ParityMapper, JordanWignerMapper, BravyiKitaevMapper
from qiskit import QuantumCircuit

###============HAMILTONIANS==================
def hamiltonian_for_molecule_and_mapper(molecule, mapper) -> QuantumCircuit:
   
    driver = PySCFDriver(
        atom=molecule,
        unit=DistanceUnit.ANGSTROM,
        basis='sto3g',
    )
    problem = driver.run()

    if mapper == "JW":
        mapper = JordanWignerMapper()
    elif mapper == "BK":
        mapper = BravyiKitaevMapper()
    elif mapper == "Parity":
        mapper = ParityMapper(num_particles=problem.num_particles)
    else:
        raise ValueError(f"Unknown mapper: {mapper}")
    return mapper.map(problem.hamiltonian.second_q_op())

