# In this part of the code I will try to write a script
# which can be used in the general description of
# synchrotron radiation using general equation

"""
    j_nu(mps, ν)

General equation of the volume emissivity of synchrotron radiation at frequency `ν` for parameters `mps`.

Equation from the book of [Longair (1994)](https://ui.adsabs.harvard.edu/abs/1994hea..book.....L/abstract)

```math
j_{\\nu} d\\nu = -\\left( \\frac{dE}{dt} \\right) N(E) dE \\propto \\kappa B^{\\alpha+1} \\nu^{-\\alpha} d\\nu
```
where, the average energy loss rate for ultra-relativistic electron is equal to
```math
-\\left( \\frac{dE}{dt} \\right) = \\frac{4}{3} \\sigma_T c u_B \\gamma^2 
```
N(E) is the electron distribution with power law index p equal to
```math
N(E) = \\kappa E^{-p}
```
E is the Einstein's mass equation and dE its derivative such that
```math
E = \\gamma m_e c^2
```
and
```math
\\alpha = \\frac{p-1}{2}
```
"""
function j_nu(mps, ν)
    C1 = (mps.σ_T*mps.c)/(3*mps.mu_o)
    C2 = mps.m_e * mps.c^2
    C3 = mps.m_e / (2 * pi * mps.m_e)
    α = (mps.p-1) \ 2
    j_nu = C1 * C2^(-2*α) * C3^(α-1) * mps.n_e0 * mps.B^(α+1) * ν^(-α)
    return(j_nu)
end

function S_nu(mps, ν)
    β = sqrt(1.0 - 1.0/mps.Γ^2) # Bulk Lorentz velocity
    μ_obs = cos(mps.θ)
    D = 1.0 / (mps.Γ*(1 - μ_obs*β)) # Doppler factor
    DL = (2.0*mps.c / (mps.Ho * 1.0E5 / 3.086E24 )) * (mps.z+1.0 - sqrt(mps.z+1.0)) # Luminosity distance in cm # Convert H_0 km/s/Mpc to s-1 (cgs)
    #Rg = 1.5E13 * mps.M8 * 1.989E33 # Rg in cm # Solar mass needs to be converted in g
    #R = 10^3 * Rg # Radius of the blob in cm
    Vb = (4.0/3.0) * pi * mps.radius^3 # Volume of the blob in cm^3
    S_nu = D^3 * (1.0+mps.z) * Vb * j_nu(mps, ν) / DL^2
    return(S_nu)
end

function Sync_Plot(mps)
    print("START SYNCHROTRON")

    # νg : Peak frequency (non-relativistic electron gyrofrequency)
    νg = mps.e*mps.B / (2*pi*mps.m_e)

    # Set up an array of photon frequencies
    log_nu_low = log10(mps.γ_min^2 * νg) # log_10 low frequency 
    log_nu_high = log10(mps.γ_max^2 * νg) # log_10 high frequency
    log_nu = range(log_nu_low, stop=log_nu_high, length=100) # create log frequency array
    nu_values = zeros(length(log_nu)) # Frequency array in normal scales (Hz)

    # j_nu_sync : synchrotron emissivity in ergs cm^-3 s^-1 sr^-1 epsilon^-1
    j_nu_sync = zeros(length(log_nu))

    # S_nu_sync : synchrotron flux density in ergs cm^-2 s^-1 sr^-1 epsilon^-1
    S_nu_sync = zeros(length(log_nu))

    # P_nu_sync : synchrotron spectral power flux in cgs units ergs cm^-2 s^-1 Hz^-1
    P_nu_sync = zeros(length(log_nu))

    for i in eachindex(j_nu_sync)
        ν = 10.0^log_nu[i]
        nu_values[i] = ν
        j_nu_sync[i] = j_nu(mps, ν)
        S_nu_sync[i] = S_nu(mps, ν)
        P_nu_sync[i] = ν * S_nu(mps, ν)
    end
    
    min_j_nu_sync = minimum(j_nu_sync)*10^(-1)
    max_j_nu_sync = maximum(j_nu_sync)*10

    min_S_nu_sync = minimum(S_nu_sync)*10^(-1)
    max_S_nu_sync = maximum(S_nu_sync)*10

    min_P_nu_sync = minimum(P_nu_sync)*10^(-1)
    max_P_nu_sync = maximum(P_nu_sync)*10

    print("\nSynchrotron Emissivity, j_nu_sync = ", j_nu_sync)
    print("\nSynchrotron Flux Density, S_nu_sync = ", S_nu_sync)
    print("\nSynchrotron Spectral Power Flux, P_nu_sync = ", P_nu_sync)

    # Plotting the synchrotron emissivity
    plot(nu_values, j_nu_sync, label = L"j_{syn} (\nu)", framestyle=:box, 
        title = "Synchrotron emissivity", titlefontsize = 10, fmt=:pdf,
        xminorticks = 2, yminorticks = 10, 
        xlims=(1E6, 1E27), ylims=(min_j_nu_sync, max_j_nu_sync),
        xaxis=:log10, yaxis=:log10)
    xlabel!(latexstring("\$\\nu\$ [Hz]"))
    ylabel!(latexstring("\$j_{\\nu}\$ [cgs]"))
    savefig("synchrotron_emissivity.png")

    # Plotting the synchrotron flux density
    plot(nu_values, S_nu_sync, label = L"S_{syn} (\nu)", framestyle=:box, 
        title = "Synchrotron flux density", titlefontsize = 10, fmt=:pdf,
        xminorticks = 2, yminorticks = 10, 
        xlims=(1E6, 1E27), ylims=(min_S_nu_sync, max_S_nu_sync),
        xaxis=:log10, yaxis=:log10)
    xlabel!(latexstring("\$\\nu\$ [Hz]"))
    ylabel!(latexstring("\$S_{\\nu}\$ [cgs]"))
    savefig("synchrotron_flux_density.png")

    # Plotting the synchrotron spectral power flux
    plot(nu_values, P_nu_sync, label = L"\nu S_{syn} (\nu)", framestyle=:box, 
        title = "Synchrotron spectral power flux", titlefontsize = 10, fmt=:pdf,
        xminorticks = 2, yminorticks = 10, 
        xlims=(1E6, 1E27), ylims=(min_P_nu_sync, max_P_nu_sync),
        xaxis=:log10, yaxis=:log10)
    xlabel!(latexstring("\$\\nu\$ [Hz]"))
    ylabel!(latexstring("\$\\nu S_{\\nu}\$ [cgs]"))
    savefig("synchrotron_spectral_power_flux.png")

    print("\nEND SYNCHROTRON")
end