# fit the SSC model to the data for PKS 0637-752 WK7.8

using SyncSSCModel
using SpectralFitting
using Plots

# define the model we want to fit
# includes default parameter values
# specifies which paramters are free
function PKSSSCModel(;
    K = FitParam(1.0, lower_limit = 0.5, upper_limit = 2.0),
    log_B = FitParam(log10(1.25E-5), lower_limit = -8.0, upper_limit = -3.0),
    n_e0 = FitParam(19.0, lower_limit = 0.0, upper_limit = 1.0E2),
    log_radius = FitParam(log10(1.0E22), lower_limit = 21.0, upper_limit = 23.0),
    Γ = FitParam(2.0),
    log_γ_min = FitParam(log10(2.5E3), lower_limit = 1.0),
    log_γ_max = FitParam(log10(4.0E6), upper_limit = 8.0),
    p = FitParam(2.6, lower_limit = 2.0, upper_limit = 3.5),
    log_dL = FitParam(log10(1.26E28), lower_limit = 23.0, upper_limit = 29.0),
    θ = FitParam(60.0 * pi / 180.0, lower_limit = 30*pi/180, upper_limit = 80*pi/180),
    z = FitParam(0.651),
)
    SSCModel{
        typeof(K),
        SpectralFitting.FreeParameters{(:log_B, :n_e0, :θ)},
    }(
        K, log_B, n_e0, p, log_γ_min, log_γ_max, Γ, log_radius, θ, log_dL, z,
    )
end

# read in the data values
begin
    lines = readlines(@__DIR__() * "/PKS0637-752.txt")
    number_expr = r"(-?\d*\.\d*)"
    search_expr = r"^" * number_expr * r"\s+" * number_expr
    data_stacked = map(filter(!isnothing, match.(search_expr, lines))) do m
        parse.(Float64, m.captures)
    end
    # sort just for coherence
    sort!(data_stacked)
    # flatten array
    data = reduce(hcat, data_stacked)
end

# create a dataset of frequency versus flux
dataset = SimpleDataset(
    "PKS0637data",
    10 .^ data[1, :],
    10 .^ data[2, :],
    x_units = SpectralFitting.SpectralUnits.u"Hz",
    x_err = 0.05 .* (10 .^ data[1, :]),
    y_err = 0.1 .* (10 .^ data[2, :]),
)

# define the fitting problem
model = PKSSSCModel()
prob = FittingProblem(model, dataset)

# fit using Levenberg Marquadt
res = fit(prob, LevenbergMarquadt(), autodiff = :finite)

# fit using Nelder Mead
# using OptimizationOptimJL
# res = fit(prob, ChiSquared(), NelderMead())

# plot the data and the best fit model
plot(dataset.x, dataset.y, seriestype = :scatter, xscale = :log10, yscale = :log10, mc=:red, xrange=(1e7, 1e26), yrange=(1e-17, 5e-13), legend = :topleft, label = "Data")

νrange =  10 .^ collect(range(7, 26, 100))
f = invokemodel(νrange, model, res.u)
nonzero = findall(x -> x > 0.0, f)
plot!(νrange[nonzero], f[nonzero], lc=:blue, label = "Best fit model")

# print the result prettily
display(res)
