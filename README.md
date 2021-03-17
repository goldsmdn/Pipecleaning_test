# Pipecleaning_test
Pipecleaning test for MSc project.

The initial test was to check the basic connectivity.

The "Steane_code_old.ipynb" file is a more extensive test with the qubit codes keyed in.  
It sets up two logical qubits in the Steane code and then entangles them with a logical Hadamard and a logical CX.

The Steane_code.ipynb" file is the most recent version.  It uses the circuits.py code.  The Jupyter workbook calls methods in the SteaneCodeLogicalQubit python class:
 - validate the parity check matrix
 - initialise the qubits needed to build one or two logical qubits
 - set up a logical zero
 - force in correct X or Y errors
 - set up ancillas
 - correct errors found in the ancilla.  NB at present errors can only be corrected if one qubit is simulated due to memory constaints. 
 - apply logical gates:
    - X or it flip
    - Hadamard
    - CNOT
 - decode the qubit

An appropriate noise model is applied in both cases.

An incomplete testing framework is in place and will be enhanced.  Documentation is not yet complete.
