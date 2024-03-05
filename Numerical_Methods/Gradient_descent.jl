function gradient_descent(A, b, n)
    # A -> matrix
    # b -> vector 
    # n -> number of iterations

    x = zero(b) # initialize 
    for _ in 1:n
        r = b - A * x # r is the residual
        # Take the subvector space spaned by r
        # Minimize the goal function on this line 
        # Set x to the new position
        alpha = - (x' * A * r - b' * r)/(r' * A * r)
        x = x + alpha * r
    end 
    return x
end

# Very slow

# make a positive definite symmetric matrix
using LinearAlgebra
A = randn(10, 10)
A = A + A'
A = A*A

b = randn(10)

gradient_descent(A, b, 10000)
inv(A) * b

# It works! requires n ~ 10000 for a 10x10 system 