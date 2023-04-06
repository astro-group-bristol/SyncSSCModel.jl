"""
    j_ssc(ϵ, mps)

Synchrotron Self-Compton (SSC) emissivity at photon energy `ϵ` for parameters `mps`.

For an isotropic power-law electron distribution in a randomly oriented magnetic field, the first order SSC emissivity is given by the equation 20 of [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract)

```math
j_{SSC}(\\epsilon, \\Omega; x) = \\frac{c \\sigma_T^2 n_{e0}^2 u_B r_b}{9 \\pi \\epsilon_B} \\left( \\frac{\\epsilon}{\\epsilon_B} \\right)^{\\tiny{-\\alpha}} \\log{\\Sigma_C}
```
"""
function j_ssc(ϵ, mps)
    γ = sqrt(ϵ / mps.ϵ_B)
    # α is the Spectral Index as a function of the power law index
    α = (mps.p - 1.0) / 2.0

    # Σ_c is the Compton synchrotron logarithm (see [Gould 1979](https://ui.adsabs.harvard.edu/abs/1979A%26A....76..306G/abstract))
    Σ_c = min(ϵ^-1, ϵ/mps.γ_min^2, mps.ϵ_B*mps.γ_max^2) / max(mps.ϵ_B*mps.γ_min^2, ϵ/mps.γ_max^2)
    
    if Σ_c > 1.0
        return(((mps.c * mps.σ_T^2 * mps.n_e0^2 * mps.u_B * mps.radius) / (9.0 * pi * mps.ϵ_B)) * (ϵ/mps.ϵ_B)^-α * log(Σ_c))
    else
        return(0.0)
    end
end


"""
    P_ssc(ϵ, mps)
    
Synchrotron Self-Compton (SSC) Spectral Power Flux at observed photon energy `ϵ`.

This is equivalent to `νF(ν)` and presented in equation 24 of [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract).

```math
P_{ssc}(\\epsilon, \\Omega; x) = D^{(3+α)} \\frac{c \\sigma_T^2 n_{e0}^2 u_B r_b V_b}{9 \\pi d_L^2} \\left( 1+z \\right)^{1 - \\alpha} \\left( \\frac{\\epsilon}{\\epsilon_B} \\right) \\ln{\\overline{\\Sigma_c}}
```
where ``\\ln{\\overline{\\Sigma_c}}`` is the transformed Compton-Synchrotron logarithm in equation 25
"""
function P_ssc(ϵ, mps)
    # NOTE THAT THE CALCULATION OF THE SPECTRAL POWER FLUX HERE DOES NOT DEPEND ON THE EMISSIVITY 
    # I JUST USE DIRECTLY THE EQUATION 24 IN DERMER ET AL. 1997 PAPER    
    # D(Γ,θ) Doppler factor
    # Γ is the Bulk Lorentz Factor
    # Θ is the angle between the direction of the blob's motion and the direction to the observer
    β = sqrt(1.0 - 1.0/mps.Γ^2)
    μ_obs = cos(mps.θ)
    D = 1.0 / (mps.Γ*(1.0 - μ_obs*β))

    # Need to convert photon frequency to characteristic lorentz factor gamma
    γ = sqrt(ϵ / mps.ϵ_B)

    # α is the Spectral Index as a function of the power law index
    α = (mps.p - 1.0) / 2.0

    # Luminosity Distance dL(z)
    # Updated to use luminosity distance in parameter block mps.dL
    # Convert H_0 km / s / Mpc to cm / s / cm (cgs)
    # dL = (2.0*mps.c / (mps.Ho * 1.0E5 / 3.086E24 )) * (mps.z+1.0 - sqrt(mps.z+1.0))
    # dL = 1.26E28  # This is the Luminosity Distance value for PKS 0637-752

    # Volume of the blob (Sphere with radius in cm) Vb(R)
    # Rg = 1.5E13 * mps.M8      # Rg is the Gravitational radius
    # R = 10^3 Rg               # Size of blob calculation (see Dermer et al. 1997 paper) 
    Vb = (4.0/3.0) * pi * mps.radius^3

    # Σc is the transformed Compton synchrotron logarithm in EQUATION 25 of Dermer et al. 1997 paper
    Σ_c = min(D/(ϵ*(1.0+mps.z)), mps.γ_max^2 * mps.ϵ_B, ϵ*(1.0+mps.z)/(D*mps.γ_min^2)) / max(mps.γ_min^2 * mps.ϵ_B, ϵ*(1+mps.z)/(D*mps.γ_max^2))

    if Σ_c > 1.0
        return((D^(3.0+α) * (mps.c*mps.σ_T^2*mps.n_e0^2*mps.u_B*mps.radius*Vb) * (1.0+mps.z)^(1.0-α) * (ϵ/mps.ϵ_B)^(1.0-α) * log(Σ_c)) / (9.0*pi*mps.dL^2))
    else
        return(0.0)
    end
end


"""
    comptonSpec(log_ν, mps)

Populate the compton spectrum `flux_density` with frequency bins given by `log_ν` for parameters `mps`
"""
function comptonSpec(log_ν, mps)
    ν = 10.0.^log_ν
    ϵ = mps.h*ν/(mps.m_e*mps.c^2)
    spec = zeros(mps.ν_n)
    for ix in range(1, mps.ν_n)
        # Note divide by epsilon to get to flux density for comparison with synchrotron spectrum
        spec[ix] = P_ssc(ϵ[ix], mps) / ϵ[ix]
    end
    return(spec)
end

"""
    ssc_Plot(mps)

Example synchotron self-compton plot for parameters `mps`.
"""
function ssc_Plot(mps)
    print("START SSC")
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
    # j_ssc : synchrotron self-compton emissivity in ergs cm^-3 s^-1 sr^-1 epsilon^-1
    j_ssc_values = zeros(length(log_nu))
    # P_ssc : synchrotron spectral power flux in cgs units ergs cm^-2 s^-1 Hz^-1
    P_ssc_values = zeros(length(log_nu))
    
    # Calculate j_ssc values
    # Calculate P_ssc values
    for i in eachindex(j_ssc_values)
        nu = 10.0^log_nu[i]
        nu_values[i] = log10(nu)     # Frequency values array (Hz) [This line might not really useful as it is already equal to log_nu define above]

        # need to convert photon frequency to epsilon
        ϵ = mps.h*nu/(mps.m_e*mps.c^2)

        j_ssc_values[i] = (j_ssc(ϵ, mps))
        P_ssc_values[i] = (P_ssc(ϵ, mps))
    end

    # Display the complete values for the synchrotron self-compton emissivity and the spectral power flux
    print("\nj_ssc_1 = ", j_ssc_values)
    print("\n---------------------------")
    print("\nP_ssc_1 = ", P_ssc_values)

    # VISUALIZATION I (j_ssc vs ν) & (P_ssc vs ν) normal scales
    # Plotting the synchrotron self-Compton emissivity j_ssc
    p1a = plot(
        log_nu, j_ssc_values, label = L"j_{ssc} (\nu)", 
        framestyle=:box, 
        title = "Synchrotron self-compton emissivity", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-6.0E-22, 5.0E-23), fmt=:pdf)
    xlabel!(latexstring("\$\\log_{10}\\left(\\nu \\right)\$ [Hz]"))
    ylabel!(latexstring("\$\\log_{10}\\left[j_{ssc} \\left( \\nu \\right) \\right]\$ [cgs]"))
    #savefig("all_ssc_emissivity_all.png")
    # ----
    # Plotting the synchrotron self-Compton spectral power flux  ~ νF(ν) vs FREQUENCY
    p1b = plot(
        log_nu, P_ssc_values, label = L"\nu S_{ssc} (\nu)",
        framestyle=:box,     
        title = "Synchrotron self-compton Spectral power Flux", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-2.0E-19, 2.0E-19), fmt=:pdf)
    xlabel!(latexstring("\$\\log_{10}\\left(\\nu \\right)\$ [Hz]"))
    ylabel!(latexstring("\$\\log_{10}\\left[ \\nu S_{ssc} \\left( \\nu \\right) \\right]\$ [cgs]"))
    #savefig("all_syn_spectral_power_flux_FUNC_Frequency.png")

    # POSITIVE VALUES
    log_j_ssc_values_2 = zeros(length(log_nu))
    j_ssc_values_2 = j_ssc_values
    for i in eachindex(j_ssc_values_2)
        if j_ssc_values_2[i] <= 0
            j_ssc_values_2[i] = 0
        else
            j_ssc_values_2[i] = j_ssc_values_2[i]
        end
        log_j_ssc_values_2[i] = log10(j_ssc_values_2[i])
    end

    log_P_ssc_values_2 = zeros(length(log_nu))
    P_ssc_values_2 = P_ssc_values
    for i in eachindex(P_ssc_values_2)
        if P_ssc_values_2[i] <= 0
            P_ssc_values_2[i] = 0
        else
            P_ssc_values_2[i] = P_ssc_values_2[i]
        end
        log_P_ssc_values_2[i] = log10(P_ssc_values_2[i])
    end

    print("\n---------------------------")
    print("\nj_ssc_2 = ", j_ssc_values_2)
    print("\n---------------------------")
    print("\nP_ssc_2 = ", P_ssc_values_2)
    print("\n---------------------------")
    print("\nLog_j_ssc_2 = ", log_j_ssc_values_2)
    print("\n---------------------------")
    print("\nLog_P_ssc_2 = ", log_P_ssc_values_2)

    # VISUALIZATION II (j_ssc vs ν) & (P_ssc vs ν) log scales
    # Plotting the synchrotron self-Compton emissivity j_ssc
    p2a = plot(
        log_nu, log_j_ssc_values_2, label = L"j_{ssc} (\nu)", 
        framestyle=:box, 
        title = "Synchrotron self-compton emissivity", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-40, -25), fmt=:pdf)
    xlabel!(latexstring("\$\\log_{10}\\left(\\nu \\right)\$ [Hz]"))
    ylabel!(latexstring("\$\\log_{10}\\left[j_{ssc} \\left( \\nu \\right) \\right]\$ [cgs]"))
    #savefig("II_ssc_emissivity.png")
    # ----
    # Plotting the synchrotron self-Compton spectral power flux  ~ νF(ν) vs FREQUENCY
    p2b = plot(
        log_nu, log_P_ssc_values_2, label = L"\nu S_{ssc} (\nu)",
        framestyle=:box,     
        title = "Synchrotron self-compton Spectral power Flux", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-25, -15), fmt=:pdf)
    xlabel!(latexstring("\$\\log_{10}\\left(\\nu \\right)\$ [Hz]"))
    ylabel!(latexstring("\$\\log_{10}\\left[ \\nu S_{ssc} \\left( \\nu \\right) \\right]\$ [cgs]"))
    #savefig("II_ssc_spectral_power_flux_FUNC_Frequency.png")

    print("\nEND SSC")

    # TRY TO USE SUBPLOT TO DISPLAY ALL THE PLOTS. JUST WONDERING WHAT IS THE BEST WHY IN JULIA?
    #l = @layout [a b; c d] ; plot(p1a, p1b, p3a, p3b, layout = l)
    l = @layout [a; b]
    plot(p2a, p2b, layout = l)
    #savefig("Test_Ssc_subplot.pdf")
end
