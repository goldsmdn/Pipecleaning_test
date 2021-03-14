import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, Aer

class SteaneCodeLogicalQubit(QuantumCircuit):
    """A single logical Qubit for the Steane code.  This makes a circuit to encode the logical zero"""

    def __init__(self, d, parity_check_matrix, *args, **kwargs):
        """
        Initializes a new QuantumCircuit for this logical qubit

        Args:
            d (int): Number of logical "data" qubit. This is a label only.

                def ___validate_parity_matrix___(self.__parity_check_matrix, self.__codewords):

        """
        self.__parity_check_matrix = parity_check_matrix
        # define valid codewords to check they are orthogonal to the parity check matrix.
        self.__codewords = [[0,0,0,0,0,0,0],   
                            [1,0,1,0,1,0,1],
                            [0,1,1,0,0,1,1],
                            [1,1,0,0,1,1,0],
                            [0,0,0,1,1,1,1],
                            [1,0,1,1,0,1,0],
                            [0,1,1,1,1,0,0],
                            [1,1,0,1,0,0,1]]
        
        self.__num_data = 7             #seven data qubits for the Steane code
        self.__num_ancilla = 3          #six ancilla qubits, three x and three z
        self.__num_extra_ancilla = 4    #four extra ancilla for decoding
        self.__qubit_no = str(d)
    
        self.__data = QuantumRegister(self.__num_data, "data " + self.__qubit_no)
        self.__mx = QuantumRegister(self.__num_ancilla, "ancilla_X " + self.__qubit_no)
        self.__mz = QuantumRegister(self.__num_ancilla, "ancilla_Z " + self.__qubit_no)

        # Spare ancillas (e.g. for readout)
        self.__extra_ancilla = QuantumRegister(self.__num_extra_ancilla, name="extra_ancilla" + self.__qubit_no)

        #classical registers
        self.__data_classical = ClassicalRegister(self.__num_data, "measure_data " + self.__qubit_no)
        self.__mx_classical = ClassicalRegister(self.__num_ancilla, "measure_ancilla_X " + self.__qubit_no)
        self.__mz_classical = ClassicalRegister(self.__num_ancilla, "measure_ancilla_Z " + self.__qubit_no)
        self.__extra_ancilla_classical = ClassicalRegister(self.__num_extra_ancilla, "measure_extra_ancilla " + self.__qubit_no)

        # user super to inherit all methods from Quantum circuits
        super().__init__(self.__data, self.__mx, self.__mz, self.__extra_ancilla,
                        self.__data_classical, self.__mx_classical, self.__mz_classical, self.__extra_ancilla_classical)

        self.validate_parity_matrix()

    def validate_parity_matrix(self):
        """validate the parity matrix and return its size"""
        if self.__parity_check_matrix == []:
            raise ValueError('Parity check matrix must be specified')
        
        for parity_row in self.__parity_check_matrix:
            if len(parity_row) != self.__num_data:
                raise ValueError('Parity check matrix rows incorrect length')

        for codeword_row in self.__codewords:
            if len(codeword_row) != self.__num_data:
                raise ValueError("Code word rows incorrect length")
            for parity_row in self.__parity_check_matrix:
                bit_store = False
                for codeword_bit in codeword_row:
                    if codeword_bit not in [0,1]:
                        raise ValueError("Code word entries must be 0 or 1")
                    for parity_bit in parity_row:
                        if parity_bit not in [0,1]:
                            raise ValueError("Parity matrix entries must be 0 or 1")
                        bit_store = bit_store^(bool(codeword_bit) ^ bool(parity_bit))
                if bit_store:
                    raise ValueError("Code word rows must be orthogonal to the parity matrix")

    def set_up_logical_zero(self):
        """Set up logical zero for data qubit"""
        parity_matrix_totals = [ 0 for x in range(self.__num_data)] # define an empty list ready to work out parity_matrix_totals

        for parity_row in self.__parity_check_matrix:
            for index in range(self.__num_data):
                parity_matrix_totals[index] = parity_matrix_totals[index] + parity_row[index]

        count = 0
        for index in range (self.__num_data):
            if parity_matrix_totals[index] == 1:
                count = count + 1
                self.h(self.__data[index])
                for parity_row in self.__parity_check_matrix:
                    if parity_row[index] == 1:              #correct row to build ancilla from
                        for column_number in range(self.__num_data):
                            if column_number != index:
                                if parity_row[column_number] == 1:
                                        self.cx(self.__data[index],self.__data[column_number])

        if count != self.__num_ancilla:
            raise ValueError('Unable to construct matrix as parity matrix does not match the ancilla needed')    

        self.barrier()

    def force_X_error(self,qubit):
        if qubit > self.__num_data - 1 :
            raise ValueError("Qubit index must be in range of data qubits")
        if qubit < 0:
            raise ValueError("Qubit index must be in range of data qubits")
        self.x(self.__data[qubit])
        self.barrier() 

    def force_Z_error(self,qubit):
        if qubit > self.__num_data - 1 :
            raise ValueError("Qubit index must be in range of data qubits")
        if qubit < 0:
            raise ValueError("Qubit index must be in range of data qubits")
        self.z(self.__data[qubit])
        self.barrier() 

    def set_up_ancilla(self):
        #apply hadamards to all ancillas
        for index in range (self.__num_ancilla):
            self.h(self.__mx[index])
            self.h(self.__mz[index])

        #apply CNOT gates according to the parity index
        for index in range (self.__num_ancilla):
            parity_row = self.__parity_check_matrix[index]
            for column_number in range(self.__num_data):
                if parity_row[column_number] ==1:
                    self.cx(self.__mx[index],self.__data[column_number])
                    self.cz(self.__mz[index],self.__data[column_number])

        self.barrier()  

        #apply final hadamards to all ancillas
        for index in range (self.__num_ancilla):
            self.h(self.__mx[index])
            self.h(self.__mz[index])

    def logical_measure(self):
        for index in range(self.__num_ancilla):
            # need to swap measurement qubits so that the measurements match the normal format of the codewords
            self.measure(self.__mx[index],self.__mx_classical[self.__num_ancilla - index-1])
            self.measure(self.__mz[index],self.__mz_classical[self.__num_ancilla - index-1])

        for index in range(self.__num_extra_ancilla):
            self.measure(self.__extra_ancilla[index],self.__extra_ancilla_classical[self.__num_extra_ancilla - index- 1])

        for index in range(self.__num_data):
            # need to swap measurement qubits so that the measurements match the normal format of the codewords
            self.measure(self.__data[index],self.__data_classical[self.__num_data - index-1])

    def correct_errors(self):
        """ produces circuit to correct  errors.  Note, need to swap ancilla bits to match how printed out"""
        transpose_parity = self._transpose_parity()

        qubit_data = {i: {"count":0} for i in range(self.__num_data)}
        single_CX_updates_list = []
   
        for qubit in qubit_data:
        # read the relevant column of the parity check matrix
            count = 0
            bit_list = transpose_parity[qubit]
            for bits in bit_list:
                if bits == 1:
                    count = count + 1
            qubit_data.update({qubit:{"count":count}})

        for qubit in range(self.__num_data):
            bit_list = transpose_parity[qubit]
            qubit_data_item = qubit_data.get(qubit)
            count = qubit_data_item.get("count")
            if count == 1:
                for bit_index in range(self.__num_ancilla): 
                    if bit_list[bit_index] == 1:
                        self.cx(self.__mz[bit_index], self.__data[qubit])
                        single_CX_updates_list.append([qubit, bit_index])

        extra_ancilla = 0   
        #list_of_impacted_qubits = []
        for qubit in range(self.__num_data):
            bit_list = transpose_parity[qubit]
            qubit_data_item = qubit_data.get(qubit)
            count = qubit_data_item.get("count")           
            if count == 2:      
                for bit_index in range(self.__num_ancilla): 
                    first_bit = 0
                    second_bit = 0   
                    if count ==2:   # need a CCNOT gate
                        for bit_index in range(self.__num_ancilla): 
                            if bit_list[bit_index] == 1:
                                if first_bit == 0:
                                    first_bit = bit_index
                                else:
                                    second_bit = bit_index
            ## need to add a ccx gate
                self.ccx(self.__mz[first_bit], self.__mz[second_bit], self.__extra_ancilla[extra_ancilla])
                self.cx(self.__extra_ancilla[extra_ancilla], self.__data[qubit])
                for items in single_CX_updates_list:
                    other_impacted_qubit = items[0]
                    bit_index = items[1]
                    if first_bit == bit_index or second_bit == bit_index:
                        # need a CX gate to reverse out changes from count 1 gates, or these will show the wrong answer.
                        self.cx(self.__extra_ancilla[extra_ancilla], self.__data[other_impacted_qubit]) 
                extra_ancilla = extra_ancilla + 1

        for qubit in range(self.__num_data):
            qubit_data_item = qubit_data.get(qubit)
            count = qubit_data_item.get("count")  
            if count == 3:
                self.ccx(self.__extra_ancilla[0], self.__extra_ancilla[1], self.__extra_ancilla[3]) 
                self.ccx(self.__extra_ancilla[0], self.__extra_ancilla[2], self.__extra_ancilla[3]) 
                self.ccx(self.__extra_ancilla[1], self.__extra_ancilla[2], self.__extra_ancilla[3]) 
                #need to undo impact of all gates made earlier as odd parity by inspection of codes.
                #maybe later could add code to check the changes made, and show they have odd parity.
                for cx_needed in range(self.__num_data):
                    self.cx(self.__extra_ancilla[3], self.__data[cx_needed])

        self.barrier()  

    def decode(self):
        """Uncomputer setting up logical zero for data qubit"""

        parity_matrix_totals = [ 0 for x in range(self.__num_data)] # define an empty list ready to work out parity_matrix_totals

        for parity_row in self.__parity_check_matrix:
            for index in range(self.__num_data):
                parity_matrix_totals[index] = parity_matrix_totals[index] + parity_row[index]

        for index in range (self.__num_data):
            if parity_matrix_totals[index] == 1:
                for parity_row in self.__parity_check_matrix:
                    if parity_row[index] == 1:              #correct row to build ancilla from
                        for column_number in range(self.__num_data):
                            if column_number != index:
                                if parity_row[column_number] == 1:
                                        self.cx(self.__data[index],self.__data[column_number])

        for index in range (self.__num_data):
            if parity_matrix_totals[index] == 1:
                self.h(self.__data[index])

        self.barrier()  

    def _transpose_parity(self):
        """transposes the parity check matrix"""
        column = []
        parity_check_transpose = []
        for row_count in range(self.__num_data):
            for column_count in range(self.__num_ancilla):
                item = self.__parity_check_matrix[column_count][row_count] 
                column.append(item)
            parity_check_transpose.append(column)
            column = []
        return(parity_check_transpose)
