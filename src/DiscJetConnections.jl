module DiscJetConnections

using LaTeXStrings
using Plots
gr()

# Global parameters and constants in cgs units
# There might be a better way of doing this

# Physical constants
m_e = 9.1093897E-28         # m_e : electron mass in g
c = 2.99792458E10           # c : speed of light in cm s^-1
h = 6.6260755E-27           # h : Planck's constant in ergs
σ_T = 0.66524616E-24        # σ_T : Thomson cross-section in cm^2

# Parameters
B = 3.0E-6                  # B : magnetic field strength in Gauss
ϵ_B = B/4.14E13             # ϵ_B : cyclotron energy in units of m_e*c^2
u_B = B^2/8.0*pi            # u_B : magnetic field energy density

n_e0 = 500.0                # n_e0 : normalisation of electron density
p = 3.0                     # p : power law index of electron distribution function
γ_min = 1.0E2               # γ_min : minimum Lorentz factor of electrons
γ_max = 1.0E7               # γ_max : maximum Lorentz factor of electrons

export dn_e
export j_syn
export syncPlot

include("synchrotron.jl")

end # module
