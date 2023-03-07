using Parameters
using Plots

@with_kw struct MyParamModel
    @deftype Float64
    # Physical constants
    "Mass of electron in g"
    m_e = 9.1093897E-28
    "Charge of electron in esu - electrostatic system units"
    e = 4.8032068E-10
    "Speed of light in cm s^-1"
    c = 2.99792458E10
    "Planck's constant in ergs"
    h = 6.6260755E-27
    "Thomson cross section in cm^2"
    σ_T = 0.66524616E-24
    "Permeability of free space in G cm-1/2 g-1/2 s2"
    mu_o = 4.191696447656766E-10    
    # Input parameters for the physics of the problem
    "Magnetic field strength in Gauss"
    B
    "Normalisation of electron density in cm^-3"
    n_e0
    "Power law index"
    p
    "Minimum Lorentz Factor"
    γ_min
    "Maximum Lorentz Factor"
    γ_max
    "Bulk Lorentz Factor"
    Γ
    "Radius of emitting Region"
    radius
    "Angle (radians) between the direction of the blob's motion and the direction to the observer"
    θ
    "Luminosity distance"
    dL
    "Redshift"
    z   
    # Binning for problem
    "log_10 low frequency" 
    log_ν_low = 7.0
    "log_10 high frequency"
    log_ν_high = 26.0
    "number of logarithmic frequency bins"
    ν_n::Int = 100
end


module Callmodel

# no not need to import parameters within this module, as none of the
# functions use the Parameters module 

# Minor functions
"""
Cyclotron energy in umits of m_e*c^2
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
    return(((mps.c*mps.σ_T*u_B(mps))/(6.0*pi*ϵ_B(mps))) * γ(ϵ, mps) * dn_e(ϵ, mps))
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
    ϵ = mps.h*ν/(mps.m_e*mps.c^2)
    spec = zeros(mps.ν_n)
    for ix in range(1, mps.ν_n)
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
        return(((mps.c * mps.σ_T^2 * mps.n_e0^2 * u_B(mps) * mps.radius) / (9.0 * pi * ϵ_B(mps))) * (ϵ/ϵ_B(mps))^-α(mps) * log(Σ_c))
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
        return(((δ(mps))^(3.0+α(mps)) * (mps.c*mps.σ_T^2*mps.n_e0^2*u_B(mps)*mps.radius*Vb(mps)) * (1.0+mps.z)^(1.0-α(mps)) * (ϵ/ϵ_B(mps))^(1.0-α(mps)) * log(Σ_c)) / (9.0*pi*mps.dL^2))
    else
        return(0.0)
    end
end
"""
Synchrotron Self-Compton Spectral Power Flux
"""
function comptonSpec(log_ν, mps)
    ν = 10.0.^log_ν
    ϵ = mps.h*ν/(mps.m_e*mps.c^2)
    spec = zeros(mps.ν_n)
    for ix in range(1, mps.ν_n)
        spec[ix] = P_ssc(ϵ[ix], mps) / ϵ[ix]
    end
    return(spec)
end

# the export keyword only really works when 
# developing with packages and modules are pre-compiled
export MyParamModel,
 dn_e,
 γ,
 α,
 Vb,
 ϵ_B,
 u_B,
 j_syn,
 j_ssc,
 S_syn,
 syncSpec,
 P_ssc,
 comptonSpec
end

# Main programme for the Sync+SSC Model
using .Callmodel

# need mean function
# using StatsBase   # Commented this, not yet needed for this part of fitting

function S3Cmodel(ϵ, mps)
    all_logFluxDensity = zeros(100, length(mps))
    all_logνFluxDensity = zeros(100, length(mps))
    energy_ϵ = Any[]
    # all_ρ_ssc_syn = zeros(length(mps))    # Commented this, not yet needed for this part of fitting
    for x in eachindex(mps)
        # none of these variables _should_ be global
        log_ν = collect(range(mps[x].log_ν_low, stop=mps[x].log_ν_high, length=mps[x].ν_n))
        # Get synchrotron
        syncFluxDensity = Callmodel.syncSpec(log_ν, mps[x])
        # Get self-Compton
        comptonFluxDensity = Callmodel.comptonSpec(log_ν, mps[x])
        # Plot spectrum
        ix_low = mps[x].ν_n
        ix_high = 1
        for ix in range(1, mps[x].ν_n)
            if (syncFluxDensity[ix] > 0.0 || comptonFluxDensity[ix] > 0.0)
                ϵ = mps[x].h*10.0^log_ν[ix]/(mps[x].m_e*mps[x].c^2)
                append!(energy_ϵ, ϵ)    # Just to check at the end the energy values
                all_logFluxDensity[ix,x] = log10(syncFluxDensity[ix] + comptonFluxDensity[ix])
                all_logνFluxDensity[ix,x] = log10(ϵ*((syncFluxDensity[ix]) + (comptonFluxDensity[ix])))
                if ix < ix_low
                    ix_low = ix
                end
                if ix > ix_high
                    ix_high = ix
                end
            else
                all_logFluxDensity[ix,x] = -50.0
                all_logνFluxDensity[ix,x] = -50.0
            end
        end
        # all_ρ_ssc_syn[x] = mean(comptonFluxDensity) / mean(syncFluxDensity)    # Commented this, not yet needed for this part of fitting
    end    
    # return everything
    (all_logFluxDensity, all_logνFluxDensity) # all_ρ_ssc_syn)   # Commented this, not yet needed for this part of fitting
    print("Energy ϵ = ", energy_ϵ)
    print("\n Flux density Fν = ", all_logFluxDensity)
    print("\n Spectral Power Flux νFν = ", all_logνFluxDensity)
    plot(log_ν, all_logνFluxDensity, ylim=(-20,-10), xlims=(7, 26))   # To test the code
end

# Testing the model
# For PKS 0637-752 - Knot WK7.8
# mpsWK78 = [MyParamModel(B=1.25E-6, n_e0=19.0, p=2.6, γ_min=2.5E3, γ_max=4.0E6, Γ=2.0, radius=1.0E22, θ=60.0*pi/180.0, dL=1.26E28, z=0.654)]
# mps = [mpsWK78]

# For Pictor A - Western Hotspot
# For fitting, we may choose one of this model 
# and comment the one not needed, change in the mps = [mpsWHPicA_1, mpsWHPicA_2] applies.
mpsWHPicA_1 = MyParamStruct(B=3.3E-5, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=1.0E6, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=5.3, z=0.035, Ho=70)
mpsWHPicA_2 = MyParamStruct(B=5.3E-5, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=4.75E5, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=2.2, z=0.035, Ho=70)
mps = [mpsWHPicA_1, mpsWHPicA_2]

ϵ = 1.0
S3Cmodel(ϵ, mps)