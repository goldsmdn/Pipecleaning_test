import pytest
from circuits import SteaneCodeLogicalQubit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, Aer

TEST_X_QUBIT = 4
SHOTS = 100     #Number of shots to run    
SIMULATOR = Aer.get_backend('qasm_simulator')

def test_parity_validation():
    """test that a random errors get fixed"""
    parity_check_matrix =  [[0,0,0,1,1,1,1],
                            [0,1,1,0,0,1,1],
                            [1,0,1,0,1,0,1]]
    for text_X_qubit in range(6):
        qubit = SteaneCodeLogicalQubit(0, parity_check_matrix)
        qubit.set_up_logical_zero()
        qubit.force_X_error(text_X_qubit)   #force X error for testing
        qubit.set_up_ancilla()
        qubit.decode()
        result = execute(qubit, SIMULATOR, shots=SHOTS).result()
        length = len(result.get_counts(qubit))
        assert length == 1
    

       