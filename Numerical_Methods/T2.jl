function Eig_power_solver(A, iters)
    # A -> matrix
    # power -> power of A, implement tolerance
    # Should be a odd

    s = size(A)[1]

    # initialize 
    eig_val = zeros(Float64, s)
    eig_vec = []

    for i in  1:s
        x = randn(s)

        for _ in 1:iters
            x = A * x
            x = x/((x' * x) ^ (1/2))
        end 

        val = x' * A * x
        eig_val[i] = val 
        push!(eig_vec, x)
        A = A - val * x * x'
    end 
    return eig_val, eig_vec
end


using LinearAlgebra
A = randn(2, 2)
A = A + A'
A = A*A

eig_val, eig_vec = Eig_power_solver(A, 1)
eig_val_test, eig_vec_test = eigen(A)

println(eig_val, "\t", eig_vec)
println(eig_val_test, "\t",  eig_vec_test)


# JULIA IS 1 INDEXED