using CSV
using DataFrames
using WGLMakie
using Optim

alldata = CSV.read("C:/Users/pgraham1/Documents/GitHub/PSI/Numerical_Methods/T3/data.csv", DataFrame)

alldata[1:3,:]

data = alldata[5:end, :]


f = Figure()
Axis(f[1,1], title="Data", xlabel="X", ylabel="Y")
errorbars!(data.x, data.y, data.sigma_y)
# You can plot multiple things on the same plot using these commands ending with "!".
# The "!" is a common thing in Julia that means "this function modifies the state".
scatter!(data.x, data.y, markersize=10, color=:maroon)
f


function likelihood(params, x, y, sigma_y)
    # "unpack" the parameters
    (b,m) = params
    # compute the model's predictions for the y values
    y_pred = b .+ m .* x
    chi = (y_pred - y) ./ sigma_y
    like = 1 ./ (sqrt(2 * pi) * sigma_y) .* exp.(-0.5 .* chi.^2)
    return prod(like)
end;


function log_likelihood(params, x, y, sigma_y)
    (b,m) = params
    y_pred = b .+ m .* x
    chi = (y_pred .- y) ./ sigma_y
    # Here, I am omitting the 1/(sqrt(2 pi) sigma_y) term because it is constant every time!
    loglike = -0.5 .* chi .^ 2
    return sum(loglike)
end;

opt = optimize(p -> -log_likelihood(p, data.x, data.y, data.sigma_y), [50., 2.])
b_opt, m_opt = Optim.minimizer(opt)

n_b, n_m = 30,40
b_grid = LinRange(0, 100, n_b)
m_grid = LinRange(1.5, 3.5, n_m)
ll_grid = zeros(n_b, n_m)
for i in 1:n_b
    for j in 1:n_m
        ll_grid[i, j] = log_likelihood([b_grid[i], m_grid[j]], data.x, data.y, data.sigma_y)
    end
end

fig, ax, hm = heatmap(b_grid, m_grid, ll_grid, colorrange=(-25, -9))
Colorbar(fig[:, end+1], hm)
scatter!(b_opt, m_opt, color=:red)
fig


N = length(data.x)
b_jack = zeros(N)
m_jack = zeros(N)
for i in 1:N
    # Make copies of our x,y,sigma_y arrays, and then delete one data point (i) each time through the loop
    x_copy = copy(data.x)
    y_copy = copy(data.y)
    s_copy = copy(data.sigma_y)

    deleteat!(x_copy, i)
    deleteat!(y_copy, i)
    deleteat!(s_copy, i)

    # repeat the optimization (using the leave-one-data arrays)
    opt = optimize(p -> -log_likelihood(p, x_copy, y_copy, s_copy), [50., 2.])
    b,m = Optim.minimizer(opt)

    @show b,m
    b_jack[i] = b
    m_jack[i] = m
end

var_b = (N-1)/N * sum((b_jack .- b_opt).^2)
var_m = (N-1)/N * sum((m_jack .- m_opt).^2)
cov_bm = (N-1)/N * sum((b_jack .- b_opt).*(m_jack .- m_opt))