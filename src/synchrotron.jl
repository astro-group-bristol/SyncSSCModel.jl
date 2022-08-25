"""
    dn_e(γ, mps)

Differential electron density at Lorentz factor `γ` for parameters `mps`.

The number of electrons with Lorentz factors in the range ``\\gamma`` to ``\\gamma + d\\gamma`` is given by ``dn_e``
```math
dn_e = n_{e_0} \\gamma^{-p} d\\gamma
```
where
```math
\\gamma_{min} \\le \\gamma \\le \\gamma_{max}
```
see also the equation 7.20 in [Rybicki and Lightman (1979)](https://ui.adsabs.harvard.edu/abs/1979rpa..book.....R/abstract)
"""
function dn_e(γ, mps)
    if (γ < mps.γ_min) || (γ > mps.γ_max)
        dn_e = 0.0
    else
        dn_e = mps.n_e0 * γ^-mps.p
    end
    return dn_e
end


"""
    j_syn(ϵ, mps)

Synchrotron emissivity at photon energy `ϵ` for parameters `mps`.

For an isotropic electron distribution in a randomly oriented magnetic field, the synchrotron emissivity is (see equation 13 of [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract); we've changed their ``H`` to our ``B``)

```math
j_{syn}(\\epsilon, \\Omega; x) = \\frac{c \\sigma_T u_B}{6 \\pi \\epsilon_B} \\left( \\frac{\\epsilon}{\\epsilon_B} \\right)^{\\tiny{1/2}} \\small{n_e} \\left[ \\left( \\frac{\\epsilon}{\\epsilon_B}\\right)^{\\tiny{1/2}};x \\right]
```
"""
function j_syn(ϵ, mps)
    # need to convert photon frequency to characteristic lorentz factor gamma
    # *** this needs to be checked and understood - might be wrong! ***
    γ = sqrt(ϵ / mps.ϵ_B)
    return((mps.c*mps.σ_T*mps.u_B)/(6.0*pi*mps.ϵ_B)) * sqrt(ϵ/mps.ϵ_B) * dn_e(γ, mps)
end


"""
    S_syn(ϵ, mps)
    
Synchrotron Flux Density at observed photon energy `ϵ`.

This impliments equation 3 of [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract).

```math
S_{syn}(\\epsilon, \\Omega; x) = \\frac{D^3 (1+z) V_b j_{syn}(\\frac{\\epsilon (1+z)}{D}, \\Omega; x)}{d_L^2}
```
"""
function S_syn(ϵ, mps)
    # doppler factor D(Γ,θ)
    # Γ is the Bulk Lorentz Factor
    # Θ is the angle between the direction of the blob's motion and the direction to the observer
    β = sqrt(1.0 - 1.0/mps.Γ^2)
    μ_obs = cos(mps.θ)
    D = 1.0 / (mps.Γ*(1 - μ_obs*β))
    
    # Luminosity Distance dL(z)
    # Convert H_0 km / s / Mpc to cm / s / cm (cgs)
    # dL = (2.0*mps.c / (mps.Ho * 1.0E5 / 3.086E24 )) * (mps.z+1.0 - sqrt(mps.z+1.0))
    dL = 1.26E28   # This is the Luminosity Distance value for PKS 0637-752
    
    # Volume of the blob (Sphere with radius in cm) Vb(R)
    # Rg = 1.5E13 * mps.M8      # Rg is the Gravitational radius
    # R = 10^3 Rg               # Size of blob calculation (see Dermer et al. 1997 paper)
    Vb = (4.0/3.0) * pi * mps.radius^3

    return((D^3 * (1.0+mps.z) * Vb * j_syn(ϵ*(1.0+mps.z)/D, mps)) / dL^2)
end


"""
    P_syn(ϵ, mps)
    
Synchrotron Spectral Power Flux at observed photon energy `ϵ`.

This is equivalent to \\nu F(\\nu) and presented in equation 23 of [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract).

```math
P_{syn}(\\epsilon, \\Omega; x) = \\epsilon  S_{syn}
```
"""
function P_syn(ϵ, mps)
    return(ϵ * S_syn(ϵ, mps))
end


"""
    syncPlot(mps)

Example synchotron plot for parameters `mps`.
"""
function syncPlot(mps)
    print("START SYN")
    # Set up an array of photon frequencies
    # log_10 low frequency 
    log_nu_low = 7.0
    # log_10 high frequency
    log_nu_high = 26.0
    # create log frequency array
    log_nu = range(log_nu_low, stop=log_nu_high, length=100)
    nu_values = zeros(length(log_nu))    # Frequency array in normal scales (Hz)

    # Define the array as function of frequency
    # Set up array for the emissivity and the spectral power flux
    # j_syn : synchrotron emissivity in ergs cm^-3 s^-1 sr^-1 epsilon^-1
    j_syn_values = zeros(length(log_nu))
    # S_syn : synchrotron flux density in ergs cm^-2 s^-1 sr^-1 epsilon^-1
    S_syn_values = zeros(length(log_nu))
    # P_syn : synchrotron spectral power flux in cgs units ergs cm^-2 s^-1 Hz^-1
    P_syn_values = zeros(length(log_nu))
    # Energy : E=hν (Define the x-axis to the energy in [ergs], try to produce fig2 in Uchiyama et al 2005)
    Energy = zeros(length(log_nu))
    
    # Logarithm scale array
    logj_syn_values = zeros(length(log_nu))
    logS_syn_values = zeros(length(log_nu))
    logP_syn_values = zeros(length(log_nu))

    # Calculate j_syn values
    # Calculate S_syn values
    # Calculate P_syn values
    for i in eachindex(j_syn_values)
        nu = 10.0^log_nu[i]
        nu_values[i] = (nu)     # Frequency values array (Hz) [This line might not really useful as it is already equal to log_nu define above]
        
        # need to convert photon frequency to epsilon
        ϵ = mps.h*nu/(mps.m_e*mps.c^2)
        
        # print("gamma = ", gamma, " and epsilon = ", epsilon, " and epsilon_B = ", epsilon_B, "\n")
        Energy[i] = (mps.h*nu / 1.602E-12) # energy in [eV] X-AXIS
        
        j_syn_values[i] = (j_syn(ϵ, mps))
        logj_syn_values[i] = log10(j_syn(ϵ, mps)) # IF NEEDED
        
        S_syn_values[i] = (S_syn(ϵ, mps))
        logS_syn_values[i] = log10(S_syn(ϵ, mps)) # IF NEEDED
        
        P_syn_values[i] = (ϵ*S_syn(ϵ, mps))
        logP_syn_values[i] = log10(ϵ*S_syn(ϵ, mps)) # IF NEEDED
    end

    print("\nEnergy = ", Energy)
    print("\n---------------------------")
    print("\nj_syn = ", j_syn_values)
    print("\n---------------------------")
    print("\nlog_j_syn = ", logj_syn_values)
    print("\n---------------------------")
    print("\nS_syn = ", S_syn_values)
    print("\n---------------------------")
    print("\nlog_S_syn = ", logS_syn_values)
    print("\n---------------------------")
    print("\nP_syn = ", P_syn_values)
    print("\n---------------------------")
    print("\nlog_P_syn = ", logP_syn_values)
    

    # Plotting the synchrotron emissivity
    plot(
        nu_values, j_syn_values, label = L"j_{syn} (\nu)",
        framestyle=:box,
        title = "Synchrotron emissivity", titlefontsize = 10,
        xminorticks = 2, yminorticks = 10, xlims=(1E6, 1E27), ylims=(1E-21, 1E-16), fmt=:jpg,
        xaxis=:log10, yaxis=:log10)
    xlabel!(latexstring("\$\\nu\$ [Hz]"))
    ylabel!(latexstring("\$j_{syn}\$ [cgs]"))
    savefig("syn_emissivity.png")

    # Plotting the synchrotron flux density
    plot(
        nu_values, S_syn_values, label = L"S_{syn} (\nu)",
        framestyle=:box,
        title = "Synchrotron Flux density", titlefontsize = 10,
        xminorticks = 2, yminorticks = 10, xlims=(1E6, 1E27), ylims=(1E-8, 1E-3), fmt=:jpg,
        xaxis=:log10, yaxis=:log10)
    xlabel!(latexstring("\$\\nu\$ [Hz]"))
    ylabel!(latexstring("\$S_{syn}\$ [cgs]"))
    savefig("syn_flux_density.png")

    # Plotting the synchrotron spectral power flux  ~ νF(ν) vs FREQUENCY
    plot(
        nu_values, P_syn_values, label = L"\nu S_{syn} (\nu)",
        framestyle=:box, 
        title = "Synchrotron spectral power flux", titlefontsize = 10,
        xminorticks = 2, yminorticks = 10, ylims=(1E-16, 1E-13), xlims=(1E6, 1E27), fmt=:jpg,
        xaxis=:log10, yaxis=:log10)
    xlabel!(latexstring("\$\\nu\$ [Hz]"))
    ylabel!(latexstring("\$\\nu S_{syn} \\left(\\nu \\right)\$ [cgs]"))
    savefig("syn_spectral_power_flux_FUNC_Frequency.png")

    # Plotting the synchrotron spectral power flux  ~ νF(ν) vs ENERGY (Uchiyama et al. 2005)
    plot(
        Energy, P_syn_values, label = L"\nu S_{syn} (\nu)",
        framestyle=:box,
        title = "Synchrotron spectral power flux", titlefontsize = 10, fmt=:jpg,
        xlims=(1.0E-7, 1.0E7), ylims=(1.0E-16, 1.0E-13), xminorticks = 2, yminorticks = 10,
        xaxis=:log10, yaxis=:log10)
    xlabel!(latexstring("\$Energy\$ [eV]"))
    ylabel!(latexstring("\$\\nu S_{syn} \\left(\\nu \\right)\$ [cgs]"))
    savefig("syn_spectral_power_flux_FUNC_Energy.png")

    print("\nEND SYN")

    # TRY TO USE SUBPLOT TO DISPLAY ALL THE PLOTS. JUST WONDERING WHAT IS THE BEST WHY IN JULIA?

end