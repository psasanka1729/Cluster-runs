'''

This code takes a MCX gate and number of qubits as input and returns
the Grover operator.

'''


import time
import numpy as np
from scipy.sparse import identity
from scipy import sparse
from scipy.sparse import lil_matrix
import scipy.sparse.linalg
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info.operators import Operator
import sys
np_load_old = np.load

# Loads the MCX operator.
np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k, encoding='bytes')

N = 10

B = np.load('MCX_operator.npy')

MCX = np.matrix(B)

# The operator Ux.
A = np.ones((2**N, 2**N))
U_x = (2/(2**N))*A - np.identity(2**N, dtype = complex)
U_x_sp = sparse.csr_matrix(U_x)

## The operator U_0. This is neeed for the sign adjustment of Grover_reconstructed operator.
U_0 = - np.identity(2 ** N, dtype=complex) 
Target_state = '0'*N
Target_index = int(Target_state, 2)
U_0.itemset((Target_index, Target_index),1)



def Hadamard_gate(): # Hadamad gate acting on one qubit.
    
    return 1/np.sqrt(2)*np.array([[1,1],[1,-1]])

def RY(theta):
    return np.array([[np.cos(theta/2), -np.sin(theta/2)], [np.sin(theta/2), np.cos(theta/2)]])

def PauliZ():
    return np.array([[1,0],[0,-1]])

# H = RY(pi/2)*Z

def Hadamard(Qubit,Noise): 

    '''

    List below will hold gates acting on one qubit. For example, for L = 3,
    the Hadamard gate acting on the qubit 1 is given by = 1 x H x 1, where 
    x is the Kronecker product. Then, qubits_list = [1,H,1].

    ''' 

    qubits_list = [] 
    
    for i in range(N):
        
        if i == Qubit: # Qubit^th position in the list is H.
            
            qubits_list.append(np.matmul(RY(np.pi/2+Noise),PauliZ()))
            
        else: # Other gates are identity operators.
            
            qubits_list.append(np.identity(2))

    '''
    
    The following loop performs the Kronecker product.

    '''        
    
    M = sparse.csr_matrix(qubits_list[0]) # Initializes the final matrix.
    
    for g in range(1,len(qubits_list)):
        
        M = sparse.kron(M,qubits_list[g]) # kronecker product.
        
    return M

'''

This function takes a gate (matrix) acting on the Qubit-th qubit and returns the matrix.

'''

def Multi_Qubit_Gate(Gate, Qubit):

    if Qubit == 0:
        
        M = sparse.csr_matrix(Gate) # Initializes the final matrix.
        
        for i in range(1,N):
        
            M = sparse.kron(M,identity(2)) # kronecker product.
        
        
    else:
        
        M = identity(2)
        
        for i in range(1,N):
            if i == Qubit:
                M = sparse.kron(M,Gate) # kronecker product.
            else:
                M = sparse.kron(M,identity(2)) # kronecker product.
        
    return M      

XH_gates = []
for i in range(N):
    XH_gates.append(['X',i])
XH_gates.append(['H',N-1])

XHR_gates = [['H',N-1]]
for i in range(N-1,-1,-1):
    XHR_gates.append(['X',i])

def U0_reconstructed(EPSILON):

    j = 0
    
    X = np.matrix([[0,1],[1,0]])
    OrH = identity(2**N)

    for i in XH_gates:
    
        if i[0] == 'H':
        
           
            OrH = OrH * Hadamard(i[1],EPSILON*NOISE[j]) # Noise
            j += 1 
        
        elif i[0] == 'X':
        

        
            OrH = OrH*Multi_Qubit_Gate(Rx(EPSILON*NOISE[j]),i[1]) # Noise
            j += 1
            
            
    '''
    
    The MCX gate
    
    '''        
 
        
    OrHR = identity(2**N)
    
    for i in XHR_gates:
    
        if i[0] == 'H':
        
           
            OrHR = OrHR * Hadamard(i[1],EPSILON*NOISE[j]) # Noise
            j += 1
        
        elif i[0] == 'X':
        
            OrHR = OrHR*Multi_Qubit_Gate(Rx(EPSILON*NOISE[j]),i[1]) # Noise
            j += 1


    return -OrH*MCX*OrHR

NOISE = 2*(np.random.rand(10**3)-0.5)

print(np.nonzero(U0_reconstructed(0.0)-U_0))