using KrylovKit: eigsolve, linsolve

A = rand(5, 5)
b = rand(5)

vec, conv = linsolve(A, b, 0, 1)
