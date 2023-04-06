module SyncSSCModel

using Parameters
using PhysicalConstants
using Unitful
using UnitfulAstro
using SpectralFitting
# using Plots

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

# ensure input length is same as output length
SpectralFitting.Δoutput_length(::Type{<:SSCModel}) = 0

# specify how the model should be invoked
# note: passing ν not log_10(ν)
# note: returning flux not log_10(flux) 
@inline function SpectralFitting.invoke!(flux, ν, model::SSCModel)
    @. flux = evaluate_model(ν, (model,))
end

export SSCModel

end