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
    # B = 1 / β in the Dermer et al. (1997) paper
    β = sqrt(1.0 - 1.0/mps.Γ^2)
    μ_obs = cos(mps.θ)
    D = 1.0 / (mps.Γ*(1 - μ_obs/β))
    
    # Luminosity Distance dL(z)
    # Convert H_0 km / s / Mpc to cm / s / cm (cgs)
    dL = (2.0*mps.c / (mps.Ho * 1.0E5 / 3.086E24 )) * (mps.z+1.0 - sqrt(mps.z+1.0))

    # Volume of the blob (Sphere with radius in cm) Vb(R)
    # Rg = 1.5E13 * mps.M8
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
    # Set up an array of photon frequencies
    # log_10 low frequency 
    log_nu_low = 7.0
    # log_10 high frequency
    log_nu_high = 26.0
    # create log frequency array
    log_nu = range(log_nu_low, stop=log_nu_high, length=100)
    # Define the j_syn array as function of frequency

    # j_syn : synchrotron emissivity in ergs cm^-3 s^-1 sr^-1 epsilon^-1
    j_syn_values = zeros(length(log_nu))
    # S_syn : synchrotron flux density in ergs cm^-2 s^-1 sr^-1 epsilon^-1
    S_syn_values = zeros(length(log_nu))
    # P_syn : synchrotron spectral power flux in cgs units ergs cm^-2 s^-1 Hz^-1
    P_syn_values = zeros(length(log_nu))

    # Calculate j_syn values
    # Calculate S_syn values
    # Calculate P_syn values
    for i in eachindex(j_syn_values)
        nu = 10.0^log_nu[i]
        # need to convert photon frequency to epsilon
        ϵ = mps.h*nu/(mps.m_e*mps.c^2)
        
        # print("gamma = ", gamma, " and epsilon = ", epsilon, " and epsilon_B = ", epsilon_B, "\n")

        j_syn_values[i] = log10(j_syn(ϵ, mps))
        S_syn_values[i] = log10(S_syn(ϵ, mps))
        P_syn_values[i] = log10(ϵ*S_syn(ϵ, mps))
    end
    
    # Plotting the synchrotron emissivity
    # plot(log_nu, j_syn_values, label = L"j_{syn} (nu)", title = "Synchrotron emissivity", titlefontsize = 10, yaxis=:log10, xlims=(9, 19), ylims=(1.0E-20, 1.0E-12), fmt=:jpg)
    # xlabel!(L"\log_{10} (\nu)")
    # ylabel!(L"j_{syn}")
    # savefig("syn_emissivity.png")

    # Plotting the synchrotron flux density
    # plot(log_nu, S_syn_values, label = L"\epsilon S_{syn} (\nu)", title = "Synchrotron Flux density", titlefontsize = 10, ylims=(-50, 50), fmt=:jpg)
    # xlabel!(L"\log(\nu) - Hz")
    # ylabel!(L"\log(S_{syn}) - cgs")
    # savefig("syn_flux_density.png")

    # Plotting the synchrotron spectral power flux  ~ νF(ν)
    plot(
        log_nu, P_syn_values, label = L"\nu S_{syn} (\nu)", 
        title = "Synchrotron spectral power flux", titlefontsize = 10, 
        ylims=(-16, -5), fmt=:jpg
        )
    xlabel!(L"\log(\nu) [Hz]")
    ylabel!(L"\log(\nu S_{syn} (\nu)) [cgs]")
    #savefig("syn_spectral_power_flux.png")
end