using DiscJetConnections
using Parameters
using Plots
gr()
print("START")

# Tavecchio et al. 2000 (SSC)
mps1 = MyParamStruct(B=1.25E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0, Ho=50)
# Harris and Krawczynski 2002
mps2 = MyParamStruct(B=195.0E-6, radius=1.0E22, Γ=9.9, γ_min=45, γ_max=1.4E6, p=2.62, dL=1.26E28, θ=5.7*pi/180.0, n_e0=6.0E-5, Ho=50)
# Uchiyama et al. 2005 (Synchrotron, IC+Sync)
# mps3 = MyParamStruct(B=1.5E-5, radius=3.086E21, Γ=12.0, γ_min=1.78E3, γ_max=2.35E5, p=2.6, dL=1.26E28, θ=3.5*pi/180.0, n_e0=6.0E-5, Ho=71)
mps3 = MyParamStruct(B=1.5E-5, radius=3.086E21, Γ=12.0, γ_min=20, γ_max=9.8E5, p=2.6, dL=1.26E28, θ=6.5*pi/180.0, n_e0=6.0E-5, Ho=71)
# Lucchini et al. 2017 (EC/CMB)
mps4 = MyParamStruct(B=14.5E-5, radius=1.0E22, Γ=10.0, γ_min=1.78E3, γ_max=1.1E6, p=2.65, dL=1.26E28, θ=9.0*pi/180.0, n_e0=44.0E-6, Ho=70)

# Pictor A Hotspot testing different value of n_e0
# mps1 = MyParamStruct(B=2.8E-6, radius=7.7E20, Γ=1.0, γ_min=1.78E3, γ_max=4.0E6, p=2.49, dL=4.752E26, θ=23.0*pi/180.0, n_e0=485.0, z=0.035, Ho=70)
# mps2 = MyParamStruct(B=2.8E-6, radius=7.7E20, Γ=1.0, γ_min=1.78E3, γ_max=4.0E6, p=2.49, dL=4.752E26, θ=23.0*pi/180.0, n_e0=19.0, z=0.035, Ho=70)
# mps3 = MyParamStruct(B=2.8E-6, radius=7.7E20, Γ=1.0, γ_min=1.78E3, γ_max=4.0E6, p=2.49, dL=4.752E26, θ=23.0*pi/180.0, n_e0=6.0E-5, z=0.035, Ho=70)

# Set an array of mps
mps = [mps1, mps2, mps3, mps4]

# Define a frequency grid and corresponding spectral flux density
# Create array
all_logFluxDensity = zeros(100, length(mps))
all_logνFluxDensity = zeros(100, length(mps))

for x in eachindex(mps)
    global log_ν = collect(range(mps[x].log_ν_low, stop=mps[x].log_ν_high, length=mps[x].ν_n))
    # Get synchrotron
    global syncFluxDensity = syncSpec(log_ν, mps[x])
    # Get self-Compton
    global comptonFluxDensity = comptonSpec(log_ν, mps[x])

    # Plot spectrum
    global ix_low = mps[x].ν_n
    global ix_high = 1
    # logFluxDensity = zeros(mps[x].ν_n)
    # logνFluxDensity = zeros(mps[x].ν_n)
    for ix in range(1, mps[x].ν_n)
        if (syncFluxDensity[ix] > 0.0 || comptonFluxDensity[ix] > 0.0)
            ϵ = mps[x].h*10.0^log_ν[ix]/(mps[x].m_e*mps[x].c^2)
            # logFluxDensity[ix] = log10(syncFluxDensity[ix] + comptonFluxDensity[ix])
            # logνFluxDensity[ix] = log10(ϵ*(syncFluxDensity[ix] + comptonFluxDensity[ix]))
            all_logFluxDensity[ix,x] = log10(syncFluxDensity[ix] + comptonFluxDensity[ix])
            all_logνFluxDensity[ix,x] = log10(ϵ*(syncFluxDensity[ix] + comptonFluxDensity[ix]))
            if ix < ix_low
                global ix_low = ix
            end
            if ix > ix_high
                global ix_high = ix
            end
        else
            # logFluxDensity[ix] = -50.0
            # logνFluxDensity[ix] = -50.0
            all_logFluxDensity[ix,x] = -50.0
            all_logνFluxDensity[ix,x] = -50.0
        end
    end
end    

print(all_logνFluxDensity)
print("\n", all_logFluxDensity)

# Plot Spectral Power Flux PKS0637-752 WK7.8
p0_WK78 = plot(log_ν, all_logνFluxDensity, ylim=(-22,-10), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["Tavecchio+00" "Harris+02" "Uchiyama+05 (IC+Syn)" "Lucchini+17"], framestyle=:box)
p0_WK78

# Plot Flux Density PKS0637-752 WK7.8
q0_WK78 = plot(log_ν, all_logFluxDensity, ylim=(-25,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = ["Tavecchio+00" "Harris+02" "Uchiyama+05 (IC+Syn)" "Lucchini+17"], framestyle=:box)
q0_WK78


# Plot Test
# p0_Γ = plot(log_ν, all_logνFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["Γ=1" "Γ=5" "Γ=10"], framestyle=:box)#, linecolor = "blue2", line = :dash)
# p0_p = plot(log_ν, all_logνFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["p=2" "p=2.5" "p=3"], framestyle=:box)#, linecolor = "blue2", line = :dash)
# p0_ne0 = plot(log_ν, all_logνFluxDensity, ylim=(-30,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["ne_0=485" "ne_0=19" "n_e0=6E-5"], framestyle=:box)#, linecolor = "blue2", line = :dash)
# p0_test = plot(log_ν, all_logνFluxDensity, ylim=(-30,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["PKS0637-752 p=2.6" "Pictor A p=2.49" "Pictor A p=3.3"], framestyle=:box)#, linecolor = "blue2", line = :dash)

# PKS0637-752
# q = plot(log_ν, logFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)
# p = plot(log_ν, logνFluxDensity, ylim=(-16,-12), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)

# PictorA
# q = plot(log_ν, logFluxDensity, ylim=(-26,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)
# p = plot(log_ν, logνFluxDensity, ylim=(-14,-8), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)

print("Done!\n")