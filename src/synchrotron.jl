"""
    dn_e(γ, n_e0, p, γ_min, γ_max)

Differential electron density at Lotentz factor `γ` assuming a power law distribution of Lorentz factors from `γ_min` to `γ_max` with a power law index `p` and normalisation n_e0.

The number of electrons with Lorentz factors in the range ``\\gamma`` to ``\\gamma + d\\gamma`` is given by ``dn_e``
```math
dn_e = n_{e_0} \\gamma^{-p} d\\gamma
```
where
```math
\\gamma_{min} \\le \\gamma \\le \\gamma_{max}
```
"""
function dn_e(γ, n_e0, p, γ_min, γ_max)
    if (γ < γ_min) || (γ > γ_max)
        dn_e = 0.0
    else
        dn_e = n_e0 * γ^-p
    end
    return dn_e
end

"""
    j_syn(ϵ)

Synchrotron emissivity at photon energy `ϵ`.

For an isotropic electron distribution in a randomly oriented magnetic field, the synchrotron emissivity is (see [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract) equation 13; we've changed their ``H`` to our ``B``)

```math
j_{syn}(\\epsilon, \\Omega; x) = \\frac{c \\sigma_T u_B}{6 \\pi \\epsilon_B} \\left( \\frac{\\epsilon}{\\epsilon_B} \\right)^{\\tiny{1/2}} \\small{n_e} \\left[ \\left( \\frac{\\epsilon}{\\epsilon_B}\\right)^{\\tiny{1/2}};x \\right]
```

(more documentation here)
"""
function j_syn(ϵ)
    # need to convert photon frequency to characteristic lorentz factor gamma
    # *** this needs to be checked and understood - might be wrong! ***
    γ = sqrt(ϵ / ϵ_B)
    return ((c*σ_T*u_B)/(6.0*pi*ϵ_B)) * sqrt(ϵ/ϵ_B) * dn_e(γ, n_e0, p, γ_min, γ_max)
end

"""
    syncPlot()

Example synchotron plot.
"""
function syncPlot()
    # Set up an array of photon frequencies
    # log_10 low frequency 
    log_nu_low = 9.0
    # log_10 high frequency
    log_nu_high = 19.0
    # create log frequency array
    log_nu = range(log_nu_low, stop=log_nu_high, length=100)
    # Define the j_syn array as function of frequency

    # j_syn : synchrotron emissivity in ergs cm-3 s-1 sr-1 epsilon-1
    j_syn_values = zeros(length(log_nu))

    # Calculate j_syn
    # The equation looks like below but I wonder if we still have to integrate it 
    # as J_sync is in function of the photon energy epsilon, photon direction Omega and the location x.

    for i in eachindex(j_syn_values)
        nu = 10.0^log_nu[i]
        # need to convert photon frequency to epsilon
        ϵ = h*nu/(m_e*c^2)
        
        # print("gamma = ", gamma, " and epsilon = ", epsilon, " and epsilon_B = ", epsilon_B, "\n")

        j_syn_values[i] = j_syn(ϵ)
    end

    return [j_syn_values, log_nu]

    # Plotting the synchrotron emissivity
    # The feature of the plot below may not be the required feature for the synchrotron emissivity
    # But I just tested how to plot with Julia
    # Values of some parameters still need to be fixed

    plot(log_nu, j_syn_values, label = L"j_{syn} (nu)", title = "Synchrotron emissivity", yaxis=:log10, xlims=(9, 19), ylims=(1.0E-20, 1.0E-12), fmt=:jpg)
    xlabel!(L"\log_{10} (\nu)")
    ylabel!(L"j_{syn}")

end


"""
    S_syn(j_syn, z)

Synchrotron Flux Density.

```math
S_{syn}(\\epsilon, \\Omega; x) = \\frac{D^3 * (1+z) * Vb * j_{syn}}{dL^2}
```
"""
function S_syn()
    # Redshift Values z
    #z = [z1, z2, z3]

    # synchrotron Flux Density
    j_syn_values = syncPlot()[1]
    S_syn_values = zeros(length(j_syn_values))
    
    # doppler factor D(Γ,θ)
    # Γ is the Bulk Lorentz Factor
    # Θ is the angle between the direction of the blob's motion and the direction to the observer
    β = 1 - (1 / 2*Γ^2)
    μ_obs = 1 - (θ^2 / 2)
    D = 1 / (Γ*(1-(β * μ_obs)))
    
    # Luminosity Distance dL(z)
    dL = (2*c / Ho) * (z+1 - sqrt(z+1))

    # Volume of the blob (Sphere with radius R) Vb(R)
    Rg = 1.5E13 * M8
    R = 10^3 * Rg    
    Vb = 4/3 * pi * R^3

    for i in eachindex(S_syn_values)
        S_syn_values[i] = (D^3 * (1+z) * Vb * j_syn_values[i]) / (dL^2)
    end
    
    #return S_syn_values

    log_nu = syncPlot()[2]

    #plot(log_nu, S_syn_values, label = L"S_{syn} (nu)", title = "Flux Density", yaxis=:log10, xlims=(9, 19), ylims=(1.0E-20, 1.0E-12), fmt=:jpg)
    plot(log_nu, S_syn_values, label = L"S_{syn} (\nu)", title = "Flux Density", fmt=:jpg)
    xlabel!(L"\log_{10} (\nu)")
    ylabel!(L"S_{syn}")

end