module DiscJetConnections

using Parameters
using LaTeXStrings
using HCubature
using Plots
gr()   # default backend for Plots

# synchrotron.jl
export dn_e
export j_syn
export S_syn
export P_syn
export syncSpec
export syncPlot

# Compton.jl
export j_ssc
export P_ssc
export comptonSpec
export ssc_Plot

# Synchrotron_General.jl
export j_nu
export S_nu
export Sync_Plot

export MyParamStruct

# Parameters and constants in cgs units
@with_kw struct MyParamStruct
    # set the default type of all parameters
    @deftype Float64

    # Physical constants
    "Mass of electron in g"
    m_e = 9.1093897E-28
    "Speed of light in cm s^-1"
    c = 2.99792458E10
    "Planck's constant in ergs"
    h = 6.6260755E-27
    "Thomson cross section in cm^2"
    σ_T = 0.66524616E-24
    "Permeability of free space in G cm-1/2 g-1/2 s2"
    mu_o = 4.191696447656766E-10
    
    # Parameters specific to the physics of the problem
    "Magnetic field strength in Gauss"
    B = 1.5E-5 # PKS0637-752
    #B = 2.8E-6 # Pictor A  
    #B = [1E-6, 1E-4, 1E-2] # Just to test for a list of B-values
    # PKS0637-752: Values taken from Schwartz et al. 2000 or Tavecchio et al. 2000 or Uchiyama et al. 2005
    "Normalisation of electron density in cm^-3"
    n_e0 = 6.0E-5 # # PKS0637-752
    #n_e0 = 485 # Pictor A
    #n_e0 = [1E-5, 10, 450] # Just to test for a list of n_e0-values
    "Power law index of electron distribution function"
    p = 2.7 # PKS0637-752
    #p = 3.3 # Pictor A
    #p = [2, 2.5, 3] # Just to test for a list of B-values            
    "Minimum Lorentz factor of electrons"
    γ_min = 1.78E3  # PKS0637-752 # Calculated using the expression in section 3.1 Uchiyama et al. 2005
    #γ_min = 1.0E2 # Pictor A
    "Maximum Lorentz factor of electrons"
    γ_max = 3.21E5  # PKS0637-752 # Calculated using the expression in section 3.1 Uchiyama et al. 2005
    #γ_max = 4E6 # Pictor A
    "Redshift"
    z = 0.651 # PKS0637-752
    #z = 0.035 # Pictor A
    "Luminosity distance in cm"
    dL = 1.26E28
    "Bulk Lorentz factor"
    Γ = 12.0 # PKS0637-752
    #Γ = 1 # Pictor A
    #Γ = [1, 5, 10] # Just to test for a list of Γ-values
    "Angle (radians) between the direction of the blob's motion and the direction to the observer"
    θ = 3.5*pi/180.0 # PKS0637-752
    #θ = 23*pi/180.0 # Pictor A
    #θ = [0, 5, 10] # Just to test for a list of θ-values
    "Radius of emitting region"
    radius = 3.086E21 # PKS0637-752
    #radius = 7.7E20 # Pictor A
    "Mass of BH in Solar Masses"
    M_BH = 1.0E8

    "Hubble parameter"
    ho = 0.67
    "Hubble constant in km s^-1 Mpc^-1 in cgs units the Hubble constant should be in the units of s^-1"
    # Hubble Constant = (50-100) km s^-1 Mpc^-1 = (1.6-3.2) 10^-18 s^-1
    # Ho = 100.0*ho
    Ho = 71.0 # PKS0637-752
    #Ho = 50 # Pictor A

    # Calculated constants needed needed for calculations
    # I think we should move these to a function rather than have them as "parameters" because they could get out of sync
    "Cyclotron energy in units of m_e*c^2"
    ϵ_B = B/4.414E13    # Dermer 1995 below eq.8
    "Magnetic field energy density"
    u_B = B^2/(8.0*pi)

    # Binning for problem
    # log_10 low frequency 
    log_ν_low = 7.0
    # log_10 high frequency
    log_ν_high = 26.0
    # number of logarithmic frequency bins
    ν_n::Int = 100
end

include("synchrotron.jl")
include("Compton.jl")
include("Synchrotron_General.jl")

end # module
