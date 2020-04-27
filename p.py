from __future__ import print_function
import numpy as np

from HMM import HMM

"""
Example from ICA 4
"""

A = np.array([[0.15, 0.25,0.25,0.35],
              [0.60, 0.20,0.1,0.1],
              [0.25, 0.2,0.3,0.25],
              [0.10, 0.40,0.4,0.1]])

B = np.array([[0.60, 0.10,0.1,0.1,0.1],
              [0.10, 0.60,0.1,0.1,0.1],
              [0.1, 0.2,0.2,0.2,0.3],
              [0.0, 0.0,0.0,0.5,0.5],])

pi0 = np.array([0.25, 0.25,0.25,0.25])


seq = ['e4', 'e3', 'e2','e2', 'e0', 'e1']

model = HMM(A, B, pi0,
            states=['A', 'B','C','D'],
            emissions=['e0', 'e1','e2','e3','e4'])

res = model.forward_backward(seq)

model.print_matrix(res)
