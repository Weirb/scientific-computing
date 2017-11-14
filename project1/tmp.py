import numpy as np
from itertools import product

def create_matrix1(n):
    a = np.zeros((n,n))

    for i,j in product(range(n), range(n)):
        if i == j:
            a[i,j] = 2
        elif abs(i-j) == 1:
            a[i,j] = -1
    return a

def create_matrix2(n, m):
    a = np.zeros((n,n))

    for i,j in product(range(n), range(n)):
        if i == j:
            a[i,j] = 2 * (i+1)**2 + m
        elif abs(i-j) == 1:
            a[i,j] = -(i+1)**2
    return a

n = 10
m = 0
a = create_matrix1(n)

b = np.ones((n,1))*2.5

x = np.linalg.solve(a, b)
# print(x)

x_guess = np.zeros((n,1))
for i in range(n):
    x_guess[i] = b[i] / ( 2* (i+1)**2 + m)

# print(x_guess)

# print (np.linalg.norm(x - x_guess))

print(np.linalg.eigvals(a))