"""
    dn_e(γ, mps)

Differential electron density at Lotentz factor `γ` for parameters `mps`.

The number of electrons with Lorentz factors in the range ``\\gamma`` to ``\\gamma + d\\gamma`` is given by ``dn_e``
```math
dn_e = n_{e_0} \\gamma^{-p} d\\gamma
```
where
```math
\\gamma_{min} \\le \\gamma \\le \\gamma_{max}
```
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

For an isotropic electron distribution in a randomly oriented magnetic field, the synchrotron emissivity is (see [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract) equation 13; we've changed their ``H`` to our ``B``)

```math
j_{syn}(\\epsilon, \\Omega; x) = \\frac{c \\sigma_T u_B}{6 \\pi \\epsilon_B} \\left( \\frac{\\epsilon}{\\epsilon_B} \\right)^{\\tiny{1/2}} \\small{n_e} \\left[ \\left( \\frac{\\epsilon}{\\epsilon_B}\\right)^{\\tiny{1/2}};x \\right]
```

(more documentation here)
"""
function j_syn(ϵ, mps)
    # need to convert photon frequency to characteristic lorentz factor gamma
    # *** this needs to be checked and understood - might be wrong! ***
    γ = sqrt(ϵ / mps.ϵ_B)
    return ((mps.c*mps.σ_T*mps.u_B)/(6.0*pi*mps.ϵ_B)) * sqrt(ϵ/mps.ϵ_B) * dn_e(γ, mps)
end

"""
    S_syn(ϵ, mps)
    
Synchrotron Flux Density at observed photon energy `ϵ`.

This impliments equation 3 of [dermer_nonthermal_1997](@cite).

```math
S_{syn}(\\epsilon, \\Omega; x) = \\frac{D^3 * (1+z) * Vb * j_{syn}}{dL^2}
```
"""
function S_syn(ϵ, mps)
    # doppler factor D(Γ,θ)
    # Γ is the Bulk Lorentz Factor
    # Θ is the angle between the direction of the blob's motion and the direction to the observer
    β = sqrt(1.0 - 1.0/mps.Γ^2)
    μ_obs = cos(mps.θ)
    D = 1 / (mps.Γ*(1-(β * μ_obs)))
    
    # Luminosity Distance dL(z)
    dL = (2.0*mps.c / mps.Ho) * (mps.z+1.0 - sqrt(mps.z+1.0))

    # Volume of the blob (Sphere with radius R) Vb(R)
    Rg = 1.5E13 * mps.M8
    R = 10^3 * Rg
    Vb = (4.0/3.0) * pi * R^3

    return((D^3 * (1.0+mps.z) * Vb * j_syn(ϵ*(1.0+mps.z)/D, mps)) / dL^2)
end

"""
    syncPlot(mps)

Example synchotron plot for parameters `mps`.
"""
function syncPlot(mps)
    # Set up an array of photon frequencies
    # log_10 low frequency 
    log_nu_low = 9.0
    # log_10 high frequency
    log_nu_high = 19.0
    # create log frequency array
    log_nu = range(log_nu_low, stop=log_nu_high, length=100)
    # Define the j_syn array as function of frequency

    # j_syn : synchrotron emissivity in ergs cm^-3 s^-1 sr^-1 epsilon^-1
    j_syn_values = zeros(length(log_nu))
    # S_syn : synchrotron flux density in erg cm^-2 s^-1 sr^-1 epsilon^-1
    S_syn_values = zeros(length(log_nu))

    # Calculate j_syn values
    # Calculate S_syn values
    for i in eachindex(j_syn_values)
        nu = 10.0^log_nu[i]
        # need to convert photon frequency to epsilon
        ϵ = mps.h*nu/(mps.m_e*mps.c^2)
        
        # print("gamma = ", gamma, " and epsilon = ", epsilon, " and epsilon_B = ", epsilon_B, "\n")

        j_syn_values[i] = j_syn(ϵ, mps)
        S_syn_values[i] = S_syn(ϵ, mps)
    end

    # Plotting the synchrotron emissivity or flux density
    # plot(log_nu, j_syn_values, label = L"j_{syn} (nu)", title = "Synchrotron emissivity", yaxis=:log10, xlims=(9, 19), ylims=(1.0E-20, 1.0E-12), fmt=:jpg)
    # xlabel!(L"\log_{10} (\nu)")
    # ylabel!(L"j_{syn}")

    plot(log_nu, S_syn_values, label = L"S_{syn} (\nu)", title = "Synchrotron flux density", yaxis=:log10, xlims=(9, 19), ylims=(1.0E15, 1.0E24), fmt=:jpg)
    xlabel!(L"\log_{10} (\nu)")
    ylabel!(L"S_{syn}")
end