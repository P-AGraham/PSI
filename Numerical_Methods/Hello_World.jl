# Type of variables
typeof(2)
typeof(1.0)


factorial(BigInt(100))
typeof(factorial(BigInt(100)))

# Complex numbers
(1 + 2im) + (1 - 2im)
angle(1 + 2im)


# Act a function on a array
sin.([1, 2, 3, 4])
vcat([1, 2], [4])

# Three random number (normal distribution):
A = randn(3)

# Non commutative multiplication by inverse
2/3 
3\2

# Reduce fractions 
6//9

# Loops 

xs = zeros(10)
for i in 1:10
    xs[i]  = 2*i + 1
end
xs

# All is JIT : first compiled then executed

# Functions
using Pkg

Pkg.add("BenchmarkTools")
Pkg.add("WGLMakie")

#using BenchmarkTools
using BenchmarkTools

function f(n_lim)
    primes = ones(n_lim)
    for i in 2:n_lim
        if Bool(primes[i])
            for j in 2:(n_lim ÷ i)
                primes[j * i] = 0
            end
        end
    end
    return primes
end

@btime f(10)


# Plotting 
using Pkg

Pkg.add("WGLMakie")

using WGLMakie

WGLMakie.lines(randn(100, 2))

# Structure
struct Point{T}
    x::T
    y::T
end 

p = Point(1, 2)
p.x, p.y

# Lisp (mathematica is a variation of Lisp ---> good for symbolic manipulation)
using WGLMakie

N = 60
function xy_data(x, y)
    r = sqrt(x^2 + y^2)
    r == 0.0 ? 1f0 : (sin(r)/r)
end
l = range(-10, stop = 10, length = N)
z = Float32[xy_data(x, y) for x in l, y in l]
surface(
    -1..1, -1..1, z,
    colormap = :Spectral
)



using BenchmarkTools

xs = [1, 2, 3, 4]
@benchmark ys = sin.(xs)

using LinearAlgebra
using SparseArrays

sprand(10, 10, 0.2)
sprand(100, 100, 0.2) # special representation

# Inverting a sparse matris is not 
spzeros(10, 10)
sparse(Diagonal(1:10))
sparse(Tridiagonal(2:10, 1:10, 3:11)) # entries lower/diagonal/upper
