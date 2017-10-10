from math import *
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from sympy.abc import x
from sympy.solvers import solve


class Roots(object):

    def __init__(self, func):
        self.lim = 0.00001
        self.expr = func
        self.expr2 = self.expr.diff(x)

    def bisection(self, x1, x2):
        steps = 0
        f = lambdify(x, self.expr)
        while True:
            if (f(x1)>0 and f(x2)<0) or (f(x1)<0 and f(x2)>0):
                xc = x1 + float(x2-x1)/2.0
                if (f(xc)<0 and f(x1)<0) or (f(xc)>0 and f(x1)>0):
                    x1 = xc
                if (f(xc)<0 and f(x2)<0) or (f(xc)>0 and f(x2)>0):
                    x2 = xc
            else:
                return False, steps
            if sqrt((x2-x1)**2)<= self.lim:
                return x1, steps
            steps +=1


    def NRMethod(self, x0):
        steps = 0
        f1 = lambdify(x, self.expr2)
        f = lambdify(x, self.expr)
        while True:
            d = -1* f(x0)/f1(x0)
            x1 = x0+d
            steps += 1
            if sqrt((x1-x0)**2)<= self.lim:
                return x1, steps
            x0 = x1


    def secant(self, x0, x1):
        steps = 0
        f = lambdify(x, self.expr)
        while True:
            slope = (f(float(x0))-f(float(x1)))/(x0-x1)
            y = slope*(x-x0)+f(float(x0))
            x2 = solve(y, x)
            steps += 1
            if sqrt((x1-x0)**2)<= self.lim:
                return x1, steps
            x0 = x1
            x1 = x2[0]

    def library(self):
        return solve(self.expr, x)






A = Roots(cos(x)*sin(3*x))
b, steps = A.bisection(0.1,1.2)
print(b)
print(steps)

c, steps1 = A.NRMethod(1.0)
print(c)
print(steps1)

d, steps2 = A.secant(0.1,1.2)
print(d)
print(steps2)

e = A.library()
print(e)
