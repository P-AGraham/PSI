using CSV
using DataFrames
using WGLMakie
using Optim

alldata = CSV.read("C:/Users/pgraham1/Documents/GitHub/PSI/Numerical_Methods/T3/data.csv", DataFrame)
data = alldata[5:end, :];

function log_likelihood(params, x, y, sigma_y)
    (b,m) = params
    y_pred = b .+ m .* x
    chi = (y_pred .- y) ./ sigma_y
    loglike = -0.5 .* chi .^ 2
    return sum(loglike)
end;

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
    b, m = Optim.minimizer(opt)

    @show b, m
    b_jack[i] = b
    m_jack[i] = m
end

opt = optimize(p -> -log_likelihood(p, data.x, data.y, data.sigma_y), [50., 2.])
b_opt, m_opt = Optim.minimizer(opt)

var_b = (N-1)/N * sum((b_jack .- b_opt).^2)
var_m = (N-1)/N * sum((m_jack .- m_opt).^2)
cov_bm = (N-1)/N * sum((b_jack .- b_opt).*(m_jack .- m_opt))
cov = [var_b cov_bm; cov_bm var_m]

# [row row; next_column_row nextcolumn_row]

theta = LinRange(-pi, +pi, 200);

X = cos.(theta)
Y = sin.(theta)

scov = sqrt(cov)

# X = (scov*([X Y]'))[1, :]
# Y = (scov*([X Y]'))[2, :]

XY = scov*[X Y]'

f = Figure()
Axis(f[1,1])
scatter!((+).(XY[1,:], b_opt), (+).(XY[2,:], m_opt), markersize=1)
f

