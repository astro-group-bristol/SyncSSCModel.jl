# fit the SSC model to the data for PKS 0637-752 WK7.8

using SyncSSCModel
using SpectralFitting
using Plots

# define the model we want to fit
# includes default parameter values
# specifies which paramters are free
function PKSSSCModel(;
    K = FitParam(1.0),
    B = FitParam(1.25E-6),
    n_e0 = FitParam(19.0),
    radius = FitParam(1.0E22),
    Γ = FitParam(2.0),
    γ_min = FitParam(2.5E3),
    γ_max = FitParam(4.0E6),
    p = FitParam(2.6),
    dL = FitParam(1.26E28),
    θ = FitParam(60.0 * pi / 180.0),
    z = FitParam(0.651),
)
    SSCModel{
        typeof(K),
        # SpectralFitting.FreeParameters{(:K, :B, :p, :radius, :θ, :dL, :z, :n_e0)},
        SpectralFitting.FreeParameters{(:K, :B, :p, :θ, :n_e0,)},
    }(
        K,
        B,
        n_e0,
        p,
        γ_min,
        γ_max,
        Γ,
        radius,
        θ,
        dL,
        z,
    )
end

νrange =  10 .^ collect(range(7, 26, 100))

model = PKSSSCModel()
flux = invokemodel(νrange, model)
# flux = invokemodel!(flux, νrange, model)

# find indices where flux is greater than 0
indices = findall(x -> x > 0, flux)

begin
    p = plot(νrange[indices], flux[indices], xscale=:log10, yscale=:log10, xlabel="ν (Hz)", ylabel="Flux (units)", label = "Model", legend = :topleft)
    display(p)
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
res = fit(prob, LevenbergMarquadt()) # , autodiff = :finite)
# using OptimizationOptimJL
# res = fit(prob, ChiSquared(), NelderMead())
# print the result prettily
display(res)

plot(dataset.x, dataset.y, seriestype = :scatter, xscale = :log10, yscale = :log10, mc=:red, xrange=(1e7, 1e26), yrange=(1e-17, 5e-13), legend = :topleft, label = "Data")
# plot!(res, lc=:blue)

f = invokemodel(νrange, model, res.u)
plot!(νrange, f, lc=:blue, label = "Best fit model")
