from __future__ import print_function

from tabulate import tabulate
import numpy as np
import pdb


class HMM(object):

    def __init__(self, A, B, pi0=None, states=None, emissions=None):
        """
        :param A: Transition matrix of shape (n, n) (n = number of states)
        :param B: Emission matrix of shape (n, b) (b = number of outputs)
        :param pi0: Initial State Probability vector of size n, leave blank for uniform probabilities
        :param states: State names/labels as list
        :param emissions: Emission names/labels as list
        """
        self.A = A
        self.B = B
        self.n_states = A.shape[0]
        self.n_emissions = B.shape[1]
        self.states = states
        self.emissions = emissions
        self.pi0 = pi0

        if pi0 is None:
            self.pi0 = np.full(self.n_states, 1.0 / self.n_states)

        if states is None:
            self.states = [chr(ord('A') + i) for i in range(self.n_states)]

        if emissions is None:
            self.emissions = [str(i) for i in range(self.n_emissions)]

    def print_matrix(self, M, headers=None):
        """
        Print matrix in tabular form

        :param M: Matrix to print
        :param headers: Optional headers for columns, default is state names
        :return: tabulated encoding of input matrix
        """
        headers = headers or self.states

        if M.ndim > 1:
            headers = [' '] + headers
            data = [['t={}'.format(i + 1)] + [j for j in row] for i, row in enumerate(M)]
        else:
            data = [[j for j in M]]
        print(tabulate(data, headers, tablefmt="grid", numalign="right"))
        return None

    def forward_algorithm(self, seq):
        """
        Apply forward algorithm to calculate probabilities of seq

        :param seq: Observed sequence to calculate probabilities upon
        :return: Alpha matrix with 1 row per time step
        """
        
        T = len(seq)

        # Initialize forward probabilities matrix Alpha
        Alpha = np.zeros((T, self.n_states))

        # Your implementation here
        harr = seq[0]
        Alpha[0] = self.pi0 * self.B[:, harr]
        Z = np.sum(Alpha[0])
        hrr=[]

        Alpha[0] = Alpha[0] / Z
        for i in range(1, T):
            harr = seq[i]
            hrr.append(harr)
            Alpha[i] = (Alpha[i - 1].dot(self.A)) * (self.B[:, harr])
            ARRAYS = np.sum(Alpha[i])
            hrr[0]=10
            Alpha[i] = Alpha[i] / ARRAYS
#        print_matrix(Alpha)
        return Alpha
                
    def retal(self,seq):
        Alpha = self.forward_algorithm(seq)
        return Alpha
    
    def retal(self,seq):
        Beta = self.backward_algorithm(seq)
        return Beta
    
    def retal(self,seq):
        Gamma = self.forward_backward(seq)
        return Gamma
        

    def backward_algorithm(self, seq):
        """
        Apply backward algorithm to calculate probabilities of seq

        :param seq: Observed sequence to calculate probabilities upon
        :return: Beta matrix with 1 row per timestep
        """

        T = len(seq)

        # InitialiARRAYSe backward probabilities matrix Beta
        Beta = np.zeros((T, self.n_states))

        # Your implementation here
            
        Beta = np.zeros((T, self.n_states))
        
        # Your implementation here
        Beta = np.zeros((T, self.n_states))
        hrr=[]
        Beta[T - 1] = [1] * self.n_states
        for i in range(T - 2, -1, -1):
            harr = seq[i + 1]
            hrr.append(harr)
            temp = self.B[:, harr] * Beta[i + 1]
            hrr[0]=10
            Beta[i] = self.A.dot(temp)
            
        return Beta
        

    def forward_backward(self, seq):
        """
        Applies forward-backward algorithm to seq

        :param seq: Observed sequence to calculate probabilities upon
        :return: Gamma matrix containing state probabilities for each timestamp
        :raises: ValueError on bad sequence
        """

        # Convert sequence to integers
        if all(isinstance(i, str) for i in seq):
            seq = [self.emissions.index(i) for i in seq]

        # Infer time steps
        T = len(seq)
        
        # Calculate forward probabilities matrix Alpha
        Alpha = self.forward_algorithm(seq)
        # InitialiARRAYSe backward probabilities matrix Beta
        Beta = self.backward_algorithm(seq)

        # InitialiARRAYSe Gamma matrix
        hrr=[]
        Gamma = np.zeros((T, self.n_states))
        for i in range(T):
            hrr=[]
            Gamma[i] = Alpha[i] * Beta[i]
            hrr.append(Beta[i])
            ARRAYS = np.sum(Gamma[i])
            hrr.append(Gamma[0])
            Gamma[i] = Gamma[i] / ARRAYS

        hs = []
        for i in range(T):
            idx = np.where(Gamma[i, :] == np.max(Gamma[i, :]))[0][0]
            hrr.append(Gamma[i])
            hs.append(self.states[idx])
            print(hs)
        print("Alpha\n",Alpha)
        print("Beta\n",Beta)
        return Gamma

        
