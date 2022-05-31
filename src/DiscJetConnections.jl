module DiscJetConnections

using Parameters
using LaTeXStrings
using HCubature
using Plots
gr()

export dn_e
export j_syn
export S_syn
export P_syn
export syncPlot

export MyParamStruct

# Parameters and constants in cgs units
@with_kw struct MyParamStruct
    # set the default type of all parameters
    @deftype Float64

    "Mass of electron in g"
    m_e = 9.1093897E-28
    "Speed of light in cm s^-1"
    c = 2.99792458E10
    "Planck's constant in ergs"
    h = 6.6260755E-27
    "Thomson cross section in cm^2"
    σ_T = 0.66524616E-24
    
    "Magnetic field strength in Gauss"
    B = 1.5E-5  
    "Cyclotron energy in units of m_e*c^2"
    ϵ_B = B/4.414E13    # Dermer 1995 below eq.8
    "Magnetic field energy density"
    u_B = B^2/8.0*pi

    # PKS0637-752: Values taken from Schwartz et al. 2000 or Tavecchio et al. 2000 or Uchiyama et al. 2005
    "Normalisation of electron density"
    n_e0 = 6.0E-5 # This value should not be bigger (not 500.0) but less
    "Power law index of electron distribution function"
    p = 2.7             
    "Minimum Lorentz factor of electrons"
    γ_min = 1.78E3  # Calculated using the expression in section 3.1 Uchiyama et al. 2005
    "Maximum Lorentz factor of electrons"
    γ_max = 3.21E5  # Calculated using the expression in section 3.1 Uchiyama et al. 2005
    "Redshift"
    z = 0.651
    "Bulk Lorentz factor"
    Γ = 12.0
    "Angle (radians) between the direction of the blob's motion and the direction to the observer"
    θ = 3.5*pi/180.0
    "Radius of emitting region"
    radius = 3.086E21

    "Hubble parameter"
    ho = 0.67
    "Hublle constant in km s^-1 Mpc^-1"
    # Ho = 100.0*ho
    Ho = 71.0
    "Mass of BH in Solar Masses"
    M8 = 1.0E8
end

include("synchrotron.jl")

end # module
