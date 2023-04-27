# fit the SSC model to the data for Pictor A

using SyncSSCModel
using SpectralFitting
using Plots

# define the model we want to fit
# includes default parameter values
# specifies which paramters are free
function PicASSCModel(;
    K = FitParam(1.0, lower_limit = 0.5, upper_limit = 2.0),
    log_B = FitParam(log10(3.3E-5), lower_limit = -6.0, upper_limit = -2.0),
    n_e0 = FitParam(5.3, lower_limit = 0.0, upper_limit = 1.0E2),
    log_radius = FitParam(log10(7.7E20), lower_limit = 20.0, upper_limit = 22.0),
    Γ = FitParam(1.0),
    γ_min = FitParam(8.7E1),
    γ_max = FitParam(1.0E6),
    p = FitParam(2.48),
    log_dL = FitParam(log10(4.752E26), lower_limit = 26.0, upper_limit = 28.0),
    θ = FitParam(23.0 * pi / 180.0),
    z = FitParam(0.035),
)
    SSCModel{
        typeof(K),
        # SpectralFitting.FreeParameters{(:K, :B, :p, :radius, :θ, :dL, :z, :n_e0)},
        # SpectralFitting.FreeParameters{(:K, :log_B, :p, :θ, :n_e0,)},
        SpectralFitting.FreeParameters{(:log_B, :n_e0,)},
    }(
        K,
        log_B,
        n_e0,
        p,
        γ_min,
        γ_max,
        Γ,
        log_radius,
        θ,
        log_dL,
        z,
    )
end

νrange =  10 .^ collect(range(7, 26, 100))

model = PicASSCModel()
flux = invokemodel(νrange, model)
# flux = invokemodel!(flux, νrange, model)

# find inices where flux is non zero
# nonzero = findall(x -> x > 0.0, flux)

begin
    p = plot(νrange[1:end-1], flux[1:end-1], xscale=:log10, yscale=:log10, xlabel="ν (Hz)", ylabel="Flux (units)", label = "Model", legend = :topleft)
    display(p)
end

# read in the data values
begin
    lines = readlines(@__DIR__() * "/PicA.txt")
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
    "PicAdata",
    10 .^ data[1, :],
    10 .^ data[2, :],
    x_units = SpectralFitting.SpectralUnits.u"Hz",
    x_err = 0.05 .* (10 .^ data[1, :]),
    y_err = 0.1 .* (10 .^ data[2, :]),
)

# Overplot dataset as large points, colored red, without lines
plot!(dataset.x, dataset.y, seriestype = :scatter, markersize = 3, markerstrokewidth = 0, mc=:red, label = "Data")

# create an instance of the model
# model = SSCModel()

# νrange =  10 .^ collect(range(7, 26, 100))
# f = invokemodel(νrange, model)

# plot(νrange, f)

# begin
    # base case
    # plot!(dataset.x, dataset.y)
    # plot!(νrange, f)
    # display(p)
# end

prob = FittingProblem(model, dataset)
res = fit(prob, LevenbergMarquadt(), autodiff = :finite)
# using OptimizationOptimJL
# res = fit(prob, ChiSquared(), NelderMead())

plot(dataset.x, dataset.y, seriestype = :scatter, xscale = :log10, yscale = :log10, mc=:red, xrange=(1e7, 1e26), yrange=(1e-14, 1e-11), legend = :topleft, label = "Data")
# plot!(res, lc=:blue)

f = invokemodel(νrange, model, res.u)
plot!(νrange, f, lc=:blue, label = "Best fit model")

# print the result prettily
display(res)
