#!/usr/bin/python

from math import sqrt

n = 3

# additive two vector
def addition_v_v(a, b):
    res = [0,0,0]
    for i in range(n):
        res[i] = a[i] + b[i]
    return res

# sloupcovy vector * radkovy vector
def multiply_v_v(a, b):
    sum = 0
    for i in range(n):
        sum += a[i] * b[i]
    return sum

# matrix * column vector
def muptiply_m_v(a, b):
    res = [0,0,0]
    for i in range(n):
        for j in range(n):
            res[i] += (a[i][j] * b[j])
    return res
        
# substract two vectors
def subtract_v_v(a, b):
    res = [0,0,0]
    for i in range(n):
        res[i] = a[i] - b[i]
    return res

# vector zvetseny o konstantu
def multiply_v_c(a, c):
    res = [0,0,0]
    for i in range(n):
        res[i] = a[i] * c
    return res

def vector_lenght(a):
    sum = 0
    for i in a:
        sum+=(i**2)
    return sqrt(sum)

def gem_vector_lenght(a,A):
    return sqrt(multiply_v_v(a,muptiply_m_v(A,a)))

k = 10
EPSILON = 1e-15

A = [[3,-1,1],[-1,3,-1],[1,-1,3]]
b = [-1, 7, -7]
xk = [0,0,0]
rk = subtract_v_v(b,muptiply_m_v(A,xk))
sk = rk
for i in range(k):
    ak = multiply_v_v(rk, rk) / multiply_v_v(sk, muptiply_m_v(A,sk)) 
    # print("ak:",ak)
    xk = addition_v_v(xk,multiply_v_c(sk,ak))
    # print(multiply_v_c(sk,ak))
    print(xk)
    rk_1 = subtract_v_v(rk,multiply_v_c(muptiply_m_v(A,sk),ak))
    if gem_vector_lenght(rk_1,A) < EPSILON:
        break
    bk = multiply_v_v(rk_1,rk_1)/multiply_v_v(rk,rk)
    sk = addition_v_v(rk_1,multiply_v_c(sk,bk))
    rk = rk_1
print("Result")
print(xk)