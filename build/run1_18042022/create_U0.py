#!/usr/bin/env python
# coding: utf-8

# WORKING.

'''

This programs returns the Grover operator for N qubits. First it creates the MCX gate then uses the
identity X^N H^t MCX H^t X^N = MCZ to create the operator U0, where Grover = U_x * U_0.

'''

import time
import numpy as np
from scipy.sparse import identity
import numpy as np
from scipy import sparse
from scipy.sparse import lil_matrix
import scipy.sparse.linalg
import sys

# In[63]:


N = 10

Noise_count = 0
SEED = int(sys.argv[1])
np.random.seed(SEED)
NOISE = 2*(np.random.rand(10**6)-0.5)

# In[64]:


## The operator U_x.
A = np.ones((2**N, 2**N))
U_x = (2/(2**N))*A - np.identity(2**N, dtype = complex)
U_x_sp = sparse.csr_matrix(U_x)

## The operator U_0. This is neeed for the sign adjustment of Grover_reconstructed operator.
U_0 = - np.identity(2 ** N, dtype=complex) 
Target_state = '0'*N
Target_index = int(Target_state, 2)
U_0.itemset((Target_index, Target_index),1)


## G is the Grover operator.
#G = np.matrix(np.matmul(U_x, U_0)) # U_w = U_x and U_s = U_0.


# In[65]:


'''

The following returns a multicontrolled U gate matrix.

Input  : c (list), t(integer), dagger (True/False).
Output : Matrix of the multicontrolled U gate with control qubits c and target qubit t.

'''
def MCU(c,t,U):
    
    '''
    
    A multicontrolled U gate with controls c (list) and target qubit t is given by 
    
    I x I x ... x I x I - PI1 x PI1 x ... x PI1 x PI1 + PI1 x PI1 x ... x PI1 x U.
    
    
    '''
    
    p0 = ['I']*N
    p1 = ['I']*N
    
    if type(c) == list:
        
        for i in c:
            p0[i] = 'PI_1'
            p1[i] = 'PI_1'
    else:
        p0[c] = 'PI_1'
        p1[c] = 'PI_1'
        
    p0[t] = 'I'
    p1[t] = 'U'
    
    I = np.identity(2)
    Z = np.matrix([[1,0],[0,-1]])
    X = np.matrix([[0,1],[1,0]])
    PI_0 = (I+Z)/2
    PI_1 = (I-Z)/2
    
    
    Matrices = {'I':I,'PI_0':PI_0,'U':U, 'PI_1':PI_1}
    

    PI_0_matrix = Matrices[p0[0]]
    for i in range(1,N):
        PI_0_matrix = sparse.kron(PI_0_matrix, Matrices[p0[i]])
        
    PI_1_matrix = Matrices[p1[0]]
    for i in range(1,N):
        PI_1_matrix = sparse.kron(PI_1_matrix, Matrices[p1[i]])
        
    return np.identity(2**N)-PI_0_matrix+PI_1_matrix


# In[66]:


def Rz_matrix(theta):

    return np.matrix([[np.exp(-1j*theta/2),0],[0,np.exp(1j*theta/2)]])

def Ry_matrix(theta):

    return np.matrix([[np.cos(theta/2), -np.sin(theta/2)],[np.sin(theta/2),np.cos(theta/2)]])

def Phase_matrix(alpha):
    
    return np.matrix([[np.exp(alpha*1j),0],[0,np.exp(alpha*1j)]])


# In[67]:


def Rz_matrix(theta):

    return np.matrix([[np.exp(-1j*theta/2),0],[0,np.exp(1j*theta/2)]])

def Rz(Angle, Qubit,Noise):
    
    if Qubit > N -1 :
        
        print("Qubit number exceeds N")
        
    else:    
    
        qubits_list = []
    
        for i in range(N):
        
            if i == Qubit:
            
                qubits_list.append(Rz_matrix(Angle+Noise))
            
            else:
            
                qubits_list.append(np.matrix(np.identity(2)))
    
        M = sparse.csr_matrix(qubits_list[0])
    
        for g in range(1,len(qubits_list)):
        
            M = sparse.kron(M, qubits_list[g]) # kronecker product.
        
        return M


# In[68]:


#N = 3
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


# In[69]:


'''

This function returns a singly controlled unitary gate. 

'''
X = np.matrix([[0,1],[1,0]])
def CU(c,t,Unitary,Noise):
    
    '''
    Creating the matrix PI0 (|0><0|) and PI1 (|1><1|).
    
    '''
    I = np.identity(2)
    Z = np.matrix([[1,0],[0,-1]])
    #X = np.matrix([[0,1],[1,0]])
    PI_0 = (I+Z)/2
    PI_1 = (I-Z)/2
    
    '''
    The following function returns the X gate for theta = pi. Any deviation from pi will
    result in a slightly different gate, which is used to model the noisy X gate.
    
    '''

    def Rx(Noise):
        A = np.cos((np.pi+Noise)/2)
        B = -1j*np.sin((np.pi+Noise)/2)
        return 1j*np.matrix([[A,B],[B,A]])
    
    Matrices = {'I':I,'PI_0':PI_0,'X':Rx(Noise), 'PI_1':PI_1}
    
    
    '''
    
    We will first create two lists p0 and p1 (for PI0 and PI1) with the matrices
    of the Kronecker product of PI0 and PI1.
    
    '''
    p0 = ['I']*N
    p1 = ['I']*N
    
    
    '''
    The string will be modified according to the position of the target and the control qubits.
    
    '''
    
    p0[c] = 'PI_0'
    p1[c] = 'PI_1'
    p1[t] = 'X'

    

    '''  
    Initialize the PI0 and PI1 matrices as the first elemenst of the list p0 and p1,
    then the following loop will perform the Kronecker product.
    
    '''    
    
    
    
    PI_0_matrix = Matrices[p0[0]]
    for i in range(1,N):
        PI_0_matrix = sparse.kron(PI_0_matrix, Matrices[p0[i]])
        
    PI_1_matrix = Matrices[p1[0]]
    for i in range(1,N):
        PI_1_matrix = sparse.kron(PI_1_matrix, Matrices[p1[i]])

    return PI_0_matrix+PI_1_matrix


# In[70]:


I = np.identity(2)
def Rx(Noise):
    A = np.cos((np.pi+Noise)/2)
    B = -1j*np.sin((np.pi+Noise)/2)
    return 1j*np.matrix([[A,B],[B,A]])
#np.around(CU(0,1,X,0.0).A,2)
#Rx(0.0)


# In[71]:


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


''' Opens the gates_list.txt file '''

l = []

file1 = open(str(N)+'_gates_list.txt', 'r')
Lines = file1.readlines()
 

for line in Lines:
    l.append(line.strip())


# In[74]:


gates_list = []
for i in l:
    j = i.split(",") 
    if j[0] == 'W':
        gates_list.append(['W',float(j[1]),float(j[2]),float(j[3]),float(j[4]),int(j[5]),int(j[6])])
        
        #for kk in Two_Qubit_Decomp(float(j[1]),float(j[2]),float(j[3]),float(j[4]),int(j[5]),int(j[6])):
            #gates_list.append(kk)
    elif j[0] == 'CX':
        gates_list.append(['CX',int(j[1]),int(j[2])])
    elif j[0] == 'H':
        gates_list.append(['H',int(j[1])])
    elif j[0] == 'RZ':
        gates_list.append(['RZ',float(j[1]),int(j[2])])



XH_gates = []
for i in range(N):
    XH_gates.append(['X',i])
XH_gates.append(['H',N-1])

XHR_gates = [['H',N-1]]
for i in range(N-1,-1,-1):
    XHR_gates.append(['X',i])




def Two_Qubit_Decomp(alpha, beta, delta, theta, control, target, Noise):

    
    # adding noise to the angles.
    alpha = alpha + Noise[0]
    beta  = beta  + Noise[1]
    delta = delta + Noise[2]
    theta = theta + Noise[3]


    
    '''
    
    The controlled phase gate will be decomposed using the algorithm described in Lemma 5.2
    page 11 of Elementary Gates for Quantum Computation.
    
    '''
    E = Rz_matrix(delta)*Phase_matrix(delta/2)
    
    A = Rz_matrix(alpha)*Ry_matrix(theta/2)
    
    B = Ry_matrix(-theta/2)*Rz_matrix(-(alpha+beta)/2)

    C = Rz_matrix((beta-alpha)/2)
    
    #CX = sparse.csr_matrix(np.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]))
    
    
    
    CPhase = Multi_Qubit_Gate(E, control).A  
   
    A_gate = Multi_Qubit_Gate(A, target).A  
    
    B_gate = Multi_Qubit_Gate(B, target).A
       
    C_gate = Multi_Qubit_Gate(C, target).A  
     

    return [
            ['P',E,control],
            ['A',A,target],
            ['CX',control, target],
            ['B',B,target],
            ['CX',control,target], 
            ['C',C,target]
                ]


# In[ ]:





# In[82]:


'''

The follwoing function produces X gate.

'''
def Rx(Noise):
    A = np.cos((np.pi+Noise)/2)
    B = -1j*np.sin((np.pi+Noise)/2)
    return 1j*np.matrix([[A,B],[B,A]])




# Reconstructs the MCX operator.
def MCX_reconstructed(EPSILON):


    Noise_count = 0
    X = np.matrix([[0,1],[1,0]])

            
            
    '''
    
    The MCX gate
    
    '''        
    
    OrX = identity(2**N)


    for gate1 in gates_list:

        
        if gate1[0] == 'H':
        
            
            OrX = OrX*Hadamard(gate1[1],EPSILON*NOISE[Noise_count]) # Noise
            Noise_count += 1
        
        elif gate1[0] == 'CX':
        
            
            OrX = OrX*CU(gate1[1], gate1[2], X, EPSILON*NOISE[Noise_count]) # Noise 
            Noise_count += 1
        
        elif gate1[0] == 'RZ':   
        
  
            OrX = OrX*Rz(gate1[1], gate1[2],EPSILON*NOISE[Noise_count]) # Noise  
            Noise_count += 1
        
        elif gate1[0] == 'W': 
            

            TQD = Two_Qubit_Decomp(gate1[1],gate1[2],gate1[3],gate1[4],
                                   gate1[5],gate1[6],EPSILON*NOISE[Noise_count:Noise_count+4]) # Noise
            Noise_count += 4
            for ii in TQD: 
                if ii[0] == 'CX':
                    OrX = OrX*CU(ii[1],ii[2],X,EPSILON*NOISE[Noise_count]) # Noise
                    Noise_count += 1
                else: # P,A,B,C
                    OrX = OrX*Multi_Qubit_Gate(ii[1],ii[2])            
      
    # Adjusting the phase of the matrix.
    #M2 = OrX.A
    #M2 = M2/M2[0,0]
    OrX = OrX/OrX[0,0]

    return -OrX

# Reconstructs the U_0 operator using X^N * H^t * MCX * H^t * X^N = MCZ .
def U0_reconstructed(EPSILON):

    Noise_count = 0    
    X = np.matrix([[0,1],[1,0]])
    OrH = identity(2**N)

    for i in XH_gates:
    
        if i[0] == 'H':
        
           
            OrH = OrH * Hadamard(i[1],EPSILON*NOISE[Noise_count]) # Noise
            Noise_count += 1 
        
        elif i[0] == 'X':
        

        
            OrH = OrH*Multi_Qubit_Gate(Rx(EPSILON*NOISE[Noise_count]),i[1]) # Noise
            Noise_count += 1
            
            

    OrHR = identity(2**N)
    
    for i in XHR_gates:
    
        if i[0] == 'H':
        
           
            OrHR = OrHR * Hadamard(i[1],EPSILON*NOISE[Noise_count]) # Noise
            Noise_count += 1
        
        elif i[0] == 'X':
        
            OrHR = OrHR*Multi_Qubit_Gate(Rx(EPSILON*NOISE[Noise_count]),i[1]) # Noise
            Noise_count += 1


    return -OrH*MCX_reconstructed(EPSILON)*OrHR




# In[86]:
# U0 has -1 along the diagonal except the target state which is 1.

EPSILON = 0.03
M = U0_reconstructed(EPSILON)
print(Noise_count)

M = M/M[0,0]

 
np.save(str(EPSILON)+'_U0_operator',M.A)



