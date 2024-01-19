using Optim
using WGLMakie

# Ex: [x, y] -> [sin(x)/x, sin(y)/y]

function my_func(x)
    return sin.(x) ./ x
end;

x_grid = LinRange(-10, +10, 100);

f = Figure()
Axis(f[1,1], title="Data", xlabel="X", ylabel="Y")
lines!(x_grid, my_func(x_grid))
f

# [x, y] -> sin(x)/x + sin(y)/y

function my_func_adapter(x)
    v = my_func(x)
    sum(v)
end;

opt = optimize(my_func_adapter, [-0.])
# Has to be 0. (float)

# The initial guess is -1 and it is passed as a vector (in general it would have multiple components) to the function which needs to return the value at that point and not consider the vector an object to do pointwise operations on. That is the goal of the wrapper my_func_adapter

x_opt = Optim.minimizer(opt)

f = Figure()
Axis(f[1,1], title="Optimizing", xlabel="X", ylabel="Y")
lines!(x_grid, my_func(x_grid))
scatter!(x_opt, my_func(x_opt), color=:red, markersize=20)
f


[1, 2, 3, 4].^2
typeof(my_func(2))
# optimize will find the minimum of a function (to maximize, we need to call meaximise metho)

# https://juliapackages.com/p/optim