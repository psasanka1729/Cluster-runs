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
#SEED = int(sys.argv[1])
#np.random.seed(SEED)
#NOISE = 2*(np.random.rand(10**6)-0.5)

# In[64]:


## The operator U_x.
A = np.ones((2**N, 2**N))
U_x = (2/(2**N))*A - np.identity(2**N, dtype = complex)

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



 
#np.save(str(EPSILON)+'_U0_operator',M.A)





# In[10]:


import numpy
import numpy.linalg

sigma_x=numpy.zeros((2,2),dtype=complex)
sigma_y=numpy.zeros((2,2),dtype=complex)
sigma_z=numpy.zeros((2,2),dtype=complex)
sigma_0=numpy.identity(2,dtype=complex)
sigma_x[0,1]=1.
sigma_x[1,0]=1.
sigma_y[0,1]=-1.j
sigma_y[1,0]=1.j
sigma_z[0,0]=1.
sigma_z[1,1]=-1.
sigma_plus=(sigma_x+1.j*sigma_y)/2.
sigma_minus=(sigma_x-1.j*sigma_y)/2.

def adjoint(psi):
    return psi.conjugate().transpose()

def psi_to_rho(psi):
    return numpy.outer(psi,psi.conjugate())

def exp_val(psi, op):
    return numpy.real(numpy.dot(adjoint(psi),op.dot(psi)))

def norm_sq(psi):
    return numpy.real(numpy.dot(adjoint(psi),psi))

def normalize(psi,tol=1e-9):
    ns=norm_sq(psi)**0.5
    if ns < tol:
        raise ValueError
    return psi/ns

def comm(a,b):
    return a.dot(b)-b.dot(a)

def anti_comm(a,b):
    return a.dot(b)+b.dot(a)

def is_herm(M,tol=1e-9):
    if M.shape[0]!=M.shape[1]:
        return False
    diff=M-adjoint(M)
    return max(numpy.abs(diff.flatten())) < tol

def is_unitary(M,tol=1e-9):
    if M.shape[0]!=M.shape[1]:
        return False
    diff=M.dot(adjoint(M))-numpy.identity((M.shape[0]))
    return max(numpy.abs(diff.flatten())) < tol

def eigu(U,tol=1e-9):
    (E_1,V_1)=numpy.linalg.eigh(U+adjoint(U))
    U_1=adjoint(V_1).dot(U).dot(V_1)
    H_1=adjoint(V_1).dot(U+adjoint(U)).dot(V_1)
    non_diag_lst=[]
    j=0
    while j < U_1.shape[0]:
        k=0
        while k < U_1.shape[0]:
            if j!=k and abs(U_1[j,k]) > tol:
                if j not in non_diag_lst:
                    non_diag_lst.append(j)
                if k not in non_diag_lst:
                    non_diag_lst.append(k)
            k+=1
        j+=1
    if len(non_diag_lst) > 0:
        non_diag_lst=numpy.sort(numpy.array(non_diag_lst))
        U_1_cut=U_1[non_diag_lst,:][:,non_diag_lst]
        (E_2_cut,V_2_cut)=numpy.linalg.eigh(1.j*(U_1_cut-adjoint(U_1_cut)))
        V_2=numpy.identity((U.shape[0]),dtype=V_2_cut.dtype)
        for j in range(len(non_diag_lst)):
            V_2[non_diag_lst[j],non_diag_lst]=V_2_cut[j,:]
        V_1=V_1.dot(V_2)
        U_1=adjoint(V_2).dot(U_1).dot(V_2)

    # Sort by phase
    U_1=numpy.diag(U_1)
    inds=numpy.argsort(numpy.imag(numpy.log(U_1)))

    return (U_1[inds],V_1[:,inds]) # = (U_d,V) s.t. U=V*U_d*V^\dagger


# In[11]:



# In[ ]:





# In[12]:


def Phi_F(operator): 
    
    return (1j*np.log(eigu(operator)[0])).real  # eigu(Gr)[0] = exp(-i * phi_F)


# In[ ]:


#Phi_F(Grover_reconstructed1(0.0))


# In[14]:


L = N // 2 # Length of half cut number of qubits.



'''
    The following function takes a wavefunction as input and returns its entropy.

'''

def Entropy(Wavefunction):




    # Converting the list to a numpy matrix.
    Psi = np.matrix(Wavefunction).reshape(len(Wavefunction),1) # Psi column matrix.

    # Normalizing Psi.
    Psi = Psi/np.linalg.norm(Psi)


      
    
    def psi(s):
        return Psi[(2**L)*s:(2**L)*s + 2**L]   
    
      
    '''
        psi(s_p) is a row matrix/vector. psi(s) is a column matrix/vector.      
        Dimension of rhoA is N/2 x N/2. 
        The element <s|rhoA|sp> is given by psi_sp^\dagger * psi_s.
        
    ''' 

    def rhoA(s,s_p): # <s|rho_A|s_p>

        # psi(s_p)^\dagger * psi(s) is the element of (s,s_p) of rho_AB.  
        return psi(s_p).getH() * psi(s)
    
    
    
    def rhoA_Matrix(N):
    
        M = np.zeros((N,N), dtype = complex) # 0 to N-1.
    
        '''
            rho is Hermitian, it is sufficient to calculate the elements
            above the diagonal. The the elements below the diagonal
            can be replace by the complex cpnjugate of the
            elements above the diagonal.
        '''
        for i in range(N):
            for j in range(N):
            
                if i <= j : # Above the diagonal (i,j) i<j.
                
                    M[i,j] = rhoA(i,j)[0,0]
                
                else: # Below the diagonal (i,j) i>j.
                
                    M[i,j] = np.conjugate(M[j,i])
        return M    
    
    
    '''
        w is the diagonal of the diagonalized matrix rhoA.

    '''
    
    w, v = np.linalg.eig(rhoA_Matrix(N))
    
    w = w.real

    '''
        The following loop calculates S = - sum \lamba_i * log(\lambda_i).

    '''
    
    DL = np.zeros(N) # Creating an array for log w with zeros.
    
    for i in range(len(w)):
    
        if abs(w[i]) < 1.e-8: # log of zero gives nan.
        
            pass # Leave the log(zero) element as zero.
    
        else:
        
            DL[i] = np.log(w[i])
        
    # Entropy = -Tr(rho * log(rho)).        
    return -sum(w*DL)




def Bin2Dec(BinaryNumber): # Converts binary to decimal numbers.
    return int(str(BinaryNumber),2)


def Dec2Bin(DecimalNumber): # Converts decimal to binary numbers.
    return bin(DecimalNumber).replace("0b", "")



List = [i for i in range(2**N)] 


'''
The following function converts all numbers in decimals in the above list 
 from 0 to 2^N -1 to binary.

''' 
def List_Bin(List):
    
    l = []
    
    for i in List:
        
        i_Bin = Dec2Bin(i)
              
        
        '''
        While converting numbers from decimal to binary, for example, 1
         is mapped to 1, to make sure that
        every numbers have N qubits in them, the following loop adds leading 
        zeros to make the
        length of the binary string equal to N. Now, 1 is mapped to 000.....1
         (string of length N).
        
        '''
        
        while len(i_Bin) < N: 
            
            i_Bin = '0'+i_Bin # This loop adds leading zeros.
            
        l.append(i_Bin)
        
    return l





'''
    The following function takes a binary string as input and rolls the qubits by one and
    returns the rolled string.

'''

def Roll_String(Binary_String):
    
    return Binary_String[-1] + Binary_String[:-1]







'''
    The following function takes a wavefunction as input and performs one roll
     on the qubits and
    returns the resultant wavefunction.

'''

def Psi_Roll(Inital_Psi):
    
    
    
    '''
        The following list contains all possible 2^N qubits after one roll 
        is performed on them.
        For example, the first position 0001 is changed to 1000.
    
    '''
    
    # Rolls every string in the list List by one qubit.
    Rl = [Roll_String(i) for i in List_Bin(List)] 

   

    
    ''' 
        The following list contains the previous list but in decimal numbers. For example,
        for N =4, the first position 1 is changed to 8.
        
    
    '''
    
    Rl_d = [Bin2Dec(i) for i in Rl] # Converts the rolled binary string to decimal number.


    '''
        The following loop rearranges the coefficients of Psi after rolling. 
        For example, for N = 4,
        if the first coefficient 0001 is mapped to the eighth coefficient 1000 after
         one rotation of
        the qubits. The coefficient of the rolled Psi in the i ^ th position is in the
         Rl_d[i] ^ th positon
        of the initial Psi.
    
    '''
    
    
    Psi_Rolled = []

    for i in range(2**N): 
        # Rearranging the coefficients according to the list l_d.
        Psi_Rolled.append(Inital_Psi[Rl_d[i]]) 
        
    return Psi_Rolled






'''
    The following function performs specified number of rolls Num on the qubits.

'''

def N_Rolled(Num, Initial_Psi): # Use loop only for postive N.
    
    if Num == 0:
        
        return Initial_Psi
    
    else:
    
        s = Psi_Roll(Initial_Psi) # One roll.
    
        for i in range(Num-1): # Loop performing remaining N-1 rolls.
        
            s = Psi_Roll(s)
        
        return np.matrix(s).reshape(2**N,1) # Returning the rolled wavefunction as a matrix.

def Average_Entropy(Initial_Psi):
    
    list_of_entropies = []
    
    '''
    The loop calculates all the entropies and returns a list containing them.
    
    '''
    for i in range(N):
        
        S = Entropy(N_Rolled(i, Initial_Psi))
        list_of_entropies.append(S)
        
    # Returns the average entropy    
    return sum(list_of_entropies) / len(list_of_entropies)


def V_Matrix(operator):
    return eigu(operator)[1]
    
def Array2List(Arr):
    Arr_l = Arr.tolist()
    l = []
    for i in Arr_l:
        l.append(i[0])
    return l


# In[15]:


np.random.seed(2022)
NOISE= 2*(np.random.rand(10**6)-0.5)


# In[ ]:





EPSILON = float(sys.argv[1])
M = U0_reconstructed(EPSILON)

M = M/M[0,0]

# Grover operator reconstructed with noise.
Op = np.matmul(U_x, M)


f = open(str(EPSILON)+'plot_data.txt', 'w')    
f = open(str(EPSILON)+'plot_data.txt', 'a')


X = str(EPSILON)
Y = Phi_F(Op)
V = eigu(Op)[1]
            
    # file -> epsilon phi_f entropy    
for j in range(2**N):
    f.write(X +'\t'+ str(Y[j].real)+ '\t' + 
        str(Average_Entropy(Array2List(V[:,j:j+1]))) +'\n')   

