# Linear Recurrence Solver
# 2024-05-12
# Luca Valentini
# <Description>

import sympy as sym
#from numbers import Number

# USER INPUT ---

# get coefficients from user: starting with a_{n-1}, a_{n-2}, etc. and each separated by a comma
in_coeffs = input("Enter the coefficients of the recurrence relation in the order <coefficient of n-1th term>, <coefficient of n-2th term> etc. Separate each value with a comma.\n")

# convert coefficients to floats and store them as a tuple
coeffs = [float(term) for term in tuple(in_coeffs.split(','))]

# get initial values from the user: we need as many as there are coefficients, and the users enters them in order a_0, a_1, etc.
in_values = input("Enter the initial values of the recurrence relation in the order <0th term>, <1st term> etc. Separate each value with a comma. There should be " + str(len(coeffs)) + " terms.\n")

# convert coefficients to floats and store them as a tuple
values = [float(term) for term in tuple(in_values.split(','))]

# check we have as many initial values as coefficients. If not, shorten the longer list.
if len(coeffs) > len(values):
    print("Too many coefficients/not enough initial values entered. Excess coefficients removed.")
    coeffs = coeffs[:len(values)]
elif len(values) > len(coeffs):
    print("Too many initial values/not enough coefficients entered. Excess initial values removed.")
    values = values[:len(coeffs)]

# SOLVING RECURRENCE

# The ordinary generating function (OGF) for the recurrence relation, F(x), is an infinite power series which can be expressed as a quotient of two finite-degree polynomials.

# We first construct the polynomial Q(x) which is the denominator of the aforementioned quotient. 
# This has coefficients given by the coefficients of the recurrence relation.

x, n = sym.symbols(("x", "n"))

Q =  1.0
for i in range(0, len(coeffs)):
    Q -= coeffs[i]*x**(i+1)

# We can use Q(x) to find P(x), the numerator of the quotient. Since F=P/Q, we have P = FQ. Since Q is defined in terms of the recurrence relation, FQ will have only finitely
# many non-zero terms, and these will correspond to the terms obtained by multiplying the known coefficients of F and Q given by the user inputs.

partialF =  0
for i in range(0, len(values)):
    partialF += values[i]*x**(i)

# find P and truncate it to the correct degree
P_coeffs = sym.poly(partialF*Q).coeffs()
P = sum(x**i * P_coeffs[-(i+1)] for i in range(len(values)))

# to find the form of the closed form solution, we need to factorise the denomiator Q into the form (1- r_1 x)(1- r_2 x) etc. The r_i are the roots of a different polynomial:
G = x**len(values)
for i in range(0, len(values)):
    G -= coeffs[i]*x**(len(values)-i-1)
roots = sym.roots(G)

# if these roots are r_1, r_2, ... then the closed form solution is A_1*r_1^n + A_2*r_2n + A_3*r_3^n + ...
# we find the A_i using the partial fraction decomposition of P(x)/Q(x)

# because sympy simplifies the fractions in the PFD, it's a challenge to pair up the partial fractions and the corresponding roots
# we split each term in the partial fraction a/cx-b into a tuple (a, b, c)

closed_form_solution = 0

partials = sym.apart(P/Q, full=True).doit().args
print(partials)
for term in partials:
    (a, b, c) = (sym.fraction(term)[0], -1*sym.fraction(term)[1].args[0], sym.fraction(term)[1].args[1]/x)

    #(a)/(cx-b) = (-a/b)/(1+(-c/b)x), so this corresponds to the root r_i = c/b and we set A_i = -a/b 
    
    closed_form_solution += (-a/b) * (c/b)**n

print(closed_form_solution)

def recurrence(n):
    if n == 0:
        return 3
    if n == 1:
        return 1
    else:
        return 5*recurrence(n-1) - 7*recurrence(n-2)
    
for i in range(0, 20):
    print(recurrence(i))
    print(closed_form_solution.subs([(n, i)]))