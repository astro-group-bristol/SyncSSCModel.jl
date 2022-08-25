"""
    j_ssc(ϵ, mps)

Synchrotron Self-Compton (SSC) emissivity at photon energy `ϵ` for parameters `mps`.

For an isotropic power-law electron distribution in a randomly oriented magnetic field, the first order SSC emissivity is given by the equation 20 of [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract)

```math
j_{SSC}(\\epsilon, \\Omega; x) = \\frac{c \\sigma_T^2 n_{e0}^2 u_B r_b}{9 \\pi \\epsilon_B} \\left( \\frac{\\epsilon}{\\epsilon_B} \\right)^{\\tiny{-\\alpha}} \\log{\\Sigma_C}
```
"""
function j_ssc(ϵ, mps)
    # need to convert photon frequency to characteristic lorentz factor gamma
    # *** this needs to be checked and understood - might be wrong! ***
    γ = sqrt(ϵ / mps.ϵ_B)

    # α is the Spectral Index as a function of the power law index
    α = (mps.p - 1) / 2

    # Σ_c is the Compton synchrotron logarithm (see [Gould 1979](https://ui.adsabs.harvard.edu/abs/1979A%26A....76..306G/abstract))
    Σ_c = min(ϵ^-1, ϵ/mps.γ_min^2, mps.ϵ_B*mps.γ_max^2) / max(mps.ϵ_B*mps.γ_min^2, ϵ/mps.γ_max^2)
    
    return(((mps.c * mps.σ_T^2 * mps.n_e0^2 * mps.u_B * mps.radius) / (9.0 * pi * mps.ϵ_B)) * (ϵ/mps.ϵ_B)^-α * log(Σ_c))
end


"""
    P_ssc(ϵ, mps)
    
Synchrotron Self-Compton (SSC) Spectral Power Flux at observed photon energy `ϵ`.

This is equivalent to \\nu F(\\nu) and presented in equation 24 of [Dermer et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997ApJS..109..103D/abstract).

```math
P_{syn}(\\epsilon, \\Omega; x) = D^(3+α) * \\frac{c \\sigma_T^2 n_{e0}^2 u_B r_b V_b}{9 \\pi d_L^2} \\left( 1+z \\right)^{1 - \\alpha} * \\left( \\frac{\\epsilon}{\\epsilon_B} \\right) \\ln{\\overline{\\Sigma_c}}
```

where \\ln{\\overline{\\Sigma_c}} is the transformed Compton-Synchrotron logarithm in equation 25

"""
# NOTE THAT THE CALCULATION OF THE SPECTRAL POWER FLUX HERE DOES NOT DEPEND ON THE EMISSIVITY 
# I JUST USE DIRECTLY THE EQUATION 24 IN DERMER ET AL. 1997 PAPER
function P_ssc(ϵ, mps)
    # D(Γ,θ) Doppler factor
    # Γ is the Bulk Lorentz Factor
    # Θ is the angle between the direction of the blob's motion and the direction to the observer
    β = sqrt(1.0 - 1.0/mps.Γ^2)
    μ_obs = cos(mps.θ)
    D = 1.0 / (mps.Γ*(1 - μ_obs*β))

    # Need to convert photon frequency to characteristic lorentz factor gamma
    γ = sqrt(ϵ / mps.ϵ_B)

    # α is the Spectral Index as a function of the power law index
    α = (mps.p - 1) / 2

    # Luminosity Distance dL(z)
    # Convert H_0 km / s / Mpc to cm / s / cm (cgs)
    # dL = (2.0*mps.c / (mps.Ho * 1.0E5 / 3.086E24 )) * (mps.z+1.0 - sqrt(mps.z+1.0))
    dL = 1.26E28  # This is the Luminosity Distance value for PKS 0637-752

    # Volume of the blob (Sphere with radius in cm) Vb(R)
    # Rg = 1.5E13 * mps.M8      # Rg is the Gravitational radius
    # R = 10^3 Rg               # Size of blob calculation (see Dermer et al. 1997 paper) 
    Vb = (4.0/3.0) * pi * mps.radius^3

    # Σc is the transformed Compton synchrotron logarithm in EQUATION 25 of Dermer et al. 1997 paper
    Σc = min(D/(ϵ*(1+mps.z)), mps.γ_max^2 * mps.ϵ_B, ϵ*(1+mps.z)/(D*mps.γ_min^2)) / max(mps.γ_min^2 * mps.ϵ_B, ϵ*(1+mps.z)/(D*mps.γ_max^2))

    return((D^(3+α) * (mps.c*mps.σ_T^2*mps.n_e0^2*mps.u_B*mps.radius*Vb) * (1.0+mps.z)^(1-α) * (ϵ/mps.ϵ_B)^(1-α) * log(Σc)) / (9*pi*dL^2))
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
    plot(
        log_nu, j_ssc_values, label = L"j_{ssc} (\nu)", 
        framestyle=:box, 
        title = "Synchrotron self-compton emissivity", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-6.0E-22, 5.0E-23), fmt=:jpg)
    xlabel!(L"\log_{10} (\nu)")
    ylabel!(L"\log_{10} [j_{ssc} (\nu)]")
    savefig("all_ssc_emissivity_all.png")
    # ----
    # Plotting the synchrotron self-Compton spectral power flux  ~ νF(ν) vs FREQUENCY
    plot(
        log_nu, P_ssc_values, label = L"\nu S_{syn} (\nu)",
        framestyle=:box,     
        title = "Synchrotron self-compton Spectral power Flux", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-2.0E-19, 2.0E-19), fmt=:jpg)
    xlabel!(L"\log(\nu) [Hz]")
    ylabel!(L"\log(\nu S_{syn} (\nu)) [cgs]")
    savefig("all_syn_spectral_power_flux_FUNC_Frequency.png")

    # POSITIVE VALUES Method I
    # TAKE THE VALUES MANUALLY
    Pos_nu_j = nu_values[41:87]
    Pos_j_ssc_values = j_ssc_values[41:87]
    Pos_nu_P = nu_values[46:92]
    Pos_P_ssc_values = P_ssc_values[46:92]
    logPos_j_ssc_values = zeros(length(Pos_nu_j))
    logPos_P_ssc_values = zeros(length(Pos_nu_P))

    for i in eachindex(Pos_j_ssc_values)
        logPos_j_ssc_values[i] = log10(Pos_j_ssc_values[i])
        logPos_P_ssc_values[i] = log10(Pos_P_ssc_values[i])
    end
    print("\n------------------------POSITIVE-VALUES-METHODE-I---------------------------------")
    print("\nlogPos_j_ssc = ", logPos_j_ssc_values)
    print("\n---------------------------")
    print("\nlogPos_P_ssc = ", logPos_P_ssc_values)

    # VISUALIZATION II (j_ssc vs ν) & (P_ssc vs ν) log scales method 1
    # Plotting the synchrotron self-Compton emissivity j_ssc
    plot(
        Pos_nu_j, logPos_j_ssc_values, label = L"j_{ssc} (\nu)", 
        framestyle=:box, 
        title = "Synchrotron self-compton emissivity", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-37, -29), fmt=:jpg)
    xlabel!(L"\log_{10} (\nu)")
    ylabel!(L"\log_{10} [j_{ssc} (\nu)]")
    savefig("I_ssc_emissivity.png")
    # ----
    # Plotting the synchrotron self-Compton spectral power flux  ~ νF(ν) vs FREQUENCY
    plot(
        Pos_nu_P, logPos_P_ssc_values, label = L"\nu S_{syn} (\nu)",
        framestyle=:box,     
        title = "Synchrotron self-compton Spectral power Flux", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-22, -18), fmt=:jpg)
    xlabel!(L"\log(\nu) [Hz]")
    ylabel!(L"\log(\nu S_{syn} (\nu)) [cgs]")
    savefig("I_ssc_spectral_power_flux_FUNC_Frequency.png")


    # POSITIVE VALUES Method II
    log_j_ssc_values_2 = zeros(length(log_nu))
    j_ssc_values_2 = j_ssc_values
    for i in eachindex(j_ssc_values_2)
        if j_ssc_values_2[i] <= 0
            j_ssc_values_2[i] = 1
        else
            j_ssc_values_2[i] = j_ssc_values_2[i]
        end
        log_j_ssc_values_2[i] = log10(j_ssc_values_2[i])
    end

    log_P_ssc_values_2 = zeros(length(log_nu))
    P_ssc_values_2 = P_ssc_values
    for i in eachindex(P_ssc_values_2)
        if P_ssc_values_2[i] <= 0
            P_ssc_values_2[i] = 1
        else
            P_ssc_values_2[i] = P_ssc_values_2[i]
        end
        log_P_ssc_values_2[i] = log10(P_ssc_values_2[i])
    end

    print("\n-----------------------POSITIVE-VALUES-METHODE-II------------------------------")
    print("\nj_ssc_2 = ", j_ssc_values_2)
    print("\n---------------------------")
    print("\nP_ssc_2 = ", P_ssc_values_2)
    print("\n---------------------------")
    print("\nLog_j_ssc_2 = ", log_j_ssc_values_2)
    print("\n---------------------------")
    print("\nLog_P_ssc_2 = ", log_P_ssc_values_2)

    # VISUALIZATION III (j_ssc vs ν) & (P_ssc vs ν) log scales method 2
    # Plotting the synchrotron self-Compton emissivity j_ssc
    plot(
        log_nu, log_j_ssc_values_2, label = L"j_{ssc} (\nu)", 
        framestyle=:box, 
        title = "Synchrotron self-compton emissivity", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-37, 0), fmt=:jpg)
    xlabel!(L"\log_{10} (\nu)")
    ylabel!(L"\log_{10} [j_{ssc} (\nu)]")
    savefig("II_ssc_emissivity.png")
    # ----
    # Plotting the synchrotron self-Compton spectral power flux  ~ νF(ν) vs FREQUENCY
    plot(
        log_nu, log_P_ssc_values_2, label = L"\nu S_{syn} (\nu)",
        framestyle=:box,     
        title = "Synchrotron self-compton Spectral power Flux", titlefontsize = 10,
        xminorticks= 3, xlims=(6, 27), yminorticks=10, ylims=(-22, 0), fmt=:jpg)
    xlabel!(L"\log(\nu) [Hz]")
    ylabel!(L"\log(\nu S_{syn} (\nu)) [cgs]")
    savefig("II_ssc_spectral_power_flux_FUNC_Frequency.png")

    print("\nEND SSC")

    # TRY TO USE SUBPLOT TO DISPLAY ALL THE PLOTS. JUST WONDERING WHAT IS THE BEST WHY IN JULIA?

end