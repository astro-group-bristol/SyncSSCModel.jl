using Parameters
using PhysicalConstants
using Unitful
using UnitfulAstro
using SpectralFitting
using Plots

import PhysicalConstants.CODATA2018: c_0, σ_e, m_e, h

qconvert(T::Type, q::Unitful.Quantity) = convert(T, q.val)

const h_cgs_f = qconvert(Float64, uconvert(u"erg*s", h))
const m_e_cgs_f = qconvert(Float64, uconvert(u"g", m_e))
const c_cgs_f = qconvert(Float64, uconvert(u"cm/s", c_0))
const σ_e_cgs_f = qconvert(Float64, uconvert(u"cm^2", σ_e))

# module Callmodel

# Minor functions
"""
Cyclotron energy in units of m_e*c^2
"""
function ϵ_B(mps)
    return(mps.B/4.414E13)
end

"""
Magnetic field energy density
"""
function u_B(mps)
    return((mps.B)^2/(8.0*pi))
end

"""
Lorentz factor
"""
function γ(ϵ, mps)
    return(sqrt(ϵ / ϵ_B(mps)))
end

"""
Spectral Index / power law index
"""
function α(mps)
    return((mps.p - 1.0) / 2.0)
end

"""
Doppler factor
"""
function δ(mps)
    β = sqrt(1.0 - 1.0/mps.Γ^2)
    μ_obs = cos(mps.θ)
    return(1.0 / (mps.Γ*(1 - μ_obs*β)))
end

"""
Volume Sphere
"""
function Vb(mps)
    return((4.0/3.0) * pi * mps.radius^3)
end

# Major functions
"""
Electron density distribution
"""
function dn_e(ϵ, mps)
    if (γ(ϵ, mps) < mps.γ_min) || (γ(ϵ, mps) > mps.γ_max)
        dn_e = 0.0
    else
        dn_e = mps.n_e0 * γ(ϵ, mps)^-mps.p
    end
    return dn_e
end
"""
Synchrotron Emissivity
"""
function j_syn(ϵ, mps)
    return(((c_cgs_f*σ_e_cgs_f*u_B(mps))/(6.0*pi*ϵ_B(mps))) * γ(ϵ, mps) * dn_e(ϵ, mps))
end
"""
Synchrotron Flux density
"""
function S_syn(ϵ, mps)
    return((δ(mps)^3 * (1.0+mps.z) * Vb(mps) * j_syn(ϵ*(1.0+mps.z)/δ(mps), mps)) / mps.dL^2)
end
"""
Synchrotron Spectral Power Flux
"""
function syncSpec(log_ν, mps)
    ν = 10.0.^log_ν
    ϵ = h_cgs_f*ν/(m_e_cgs_f*c_cgs_f^2)
    spec = zeros(length(log_ν))
    for ix in range(1, length(log_ν))
        spec[ix] = S_syn(ϵ[ix], mps)
    end
    return(spec)
end
"""
Synchrotron Self-Compton emissivity
"""
function j_ssc(ϵ, mps)
    Σ_c = min(ϵ^-1, ϵ/mps.γ_min^2, ϵ_B(mps)*mps.γ_max^2) / max(ϵ_B(mps)*mps.γ_min^2, ϵ/mps.γ_max^2)
    if Σ_c > 1.0
        return(((c_cgs_f * σ_e_cgs_f^2 * mps.n_e0^2 * u_B(mps) * mps.radius) / (9.0 * pi * ϵ_B(mps))) * (ϵ/ϵ_B(mps))^-α(mps) * log(Σ_c))
    else
        return(0.0)
    end
end
"""
Synchrotron Self-Compton Flux density
"""
function P_ssc(ϵ, mps)

    Σ_c = min(δ(mps)/(ϵ*(1.0+mps.z)), mps.γ_max^2 * ϵ_B(mps), ϵ*(1.0+mps.z)/(δ(mps)*mps.γ_min^2)) / max(mps.γ_min^2 * ϵ_B(mps), ϵ*(1+mps.z)/(δ(mps)*mps.γ_max^2))
    if Σ_c > 1.0
        return(((δ(mps))^(3.0+α(mps)) * (c_cgs_f*σ_e_cgs_f^2*mps.n_e0^2*u_B(mps)*mps.radius*Vb(mps)) * (1.0+mps.z)^(1.0-α(mps)) * (ϵ/ϵ_B(mps))^(1.0-α(mps)) * log(Σ_c)) / (9.0*pi*mps.dL^2))
    else
        return(0.0)
    end
end
"""
Synchrotron Self-Compton Spectral Power Flux
"""
function comptonSpec(log_ν, mps)
    ν = 10.0.^log_ν
    ϵ = h_cgs_f*ν/(m_e_cgs_f*c_cgs_f^2)
    spec = zeros(length(log_ν))
    for ix in range(1, length(log_ν))
        spec[ix] = P_ssc(ϵ[ix], mps) / ϵ[ix]
    end
    return(spec)
end

# the export keyword only really works when 
# developing with packages and modules are pre-compiled
# export S_syn, P_ssc
# end

# Main programme for the Sync+SSC Model
# using .Callmodel

# need mean function
# using StatsBase   # Commented this, not yet needed for this part of fitting

function evaluate_model(ν, mps)
    ϵ = h_cgs_f*ν/(m_e_cgs_f*c_cgs_f^2)

    # Get synchrotron
    syncFluxDensity = S_syn(ϵ, mps)
    # Get self-Compton
    comptonFluxDensity = P_ssc(ϵ, mps) / ϵ
    return (syncFluxDensity + comptonFluxDensity) * ϵ
end

function evaluate_model(νs::AbstractArray, mps)
    map(ν -> evaluate_model(ν, mps), νs)
end

@with_kw struct SSCModel{T,F} <: AbstractSpectralModel{T,Additive}
    @deftype T
    # need this field for additive models
    "Normalisation."
    K::T
    "Magnetic field strength in Gauss"
    B::T
    "Normalisation of electron density in cm^-3"
    n_e0::T
    "Power law index"
    p::T
    "Minimum Lorentz Factor"
    γ_min::T
    "Maximum Lorentz Factor"
    γ_max::T
    "Bulk Lorentz Factor"
    Γ::T
    "Radius of emitting Region"
    radius::T
    "Angle (radians) between the direction of the blob's motion and the direction to the observer"
    θ::T
    "Luminosity distance"
    dL::T
    "Redshift"
    z::T
end

function SSCModel(;
    K = FitParam(1.0),
    B = FitParam(3.3E-5),
    n_e0 = FitParam(5.3),
    radius = FitParam(7.7E20),
    Γ = FitParam(1.0),
    γ_min = FitParam(8.7E1),
    γ_max = FitParam(1.0E6),
    p = FitParam(2.48),
    dL = FitParam(4.752E26),
    θ = FitParam(23.0 * pi / 180.0),
    z = FitParam(0.035),
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

# ensure input length is same as output length
SpectralFitting.Δoutput_length(::Type{<:SSCModel}) = 0

# specify how the model should be invoked
# note: passing ν not log_10(ν)
# note: returning flux not log_10(flux) 
@inline function SpectralFitting.invoke!(flux, ν, model::SSCModel)
    @. flux = evaluate_model(ν, (model,))
end

νrange =  10 .^ collect(range(7, 26, 100))

model = SSCModel()
flux = invokemodel(νrange, model)

begin
    p = plot(νrange[1:end-1], flux[1:end-1], xscale=:log10, yscale=:log10, xlabel="ν (Hz)", ylabel="Flux (units)")
    display(p)
end

# read in the data values
begin
    lines = readlines(@__DIR__() * "/data_obs_sync_ssc.txt")
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
    "sample_data",
    10 .^ data[1, :],
    10 .^ data[2, :],
    x_units = SpectralFitting.SpectralUnits.u"Hz",
    yerr = 0.05 .* (10 .^ data[1, :]),
    xerr = 0.1 .* (10 .^ data[2, :]),
)

# Overplot dataset as large points, colored red, without lines
plot!(dataset.x, dataset.y, seriestype = :scatter, markersize = 3, markerstrokewidth = 0, mc=:red)

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
# print the result prettily
display(res)

plot(dataset, xscale = :identity, yscale = :identity, mc=:red, xrange=(1e7, 1e26), yrange=(1e-15, 1e-10))
# plot!(res, lc=:blue)

f = invokemodel(νrange, model, res.u)
plot!(νrange, f, lc=:blue)
