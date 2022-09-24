using DiscJetConnections
using Parameters
using Plots
gr()

# Define a frequency grid and corresponding spectral flux density
mps = MyParamStruct(B=1.25E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0)
log_ν = collect(range(mps.log_ν_low, stop=mps.log_ν_high, length=mps.ν_n))
# Get synchrotron
syncFluxDensity = syncSpec(log_ν, mps)
# Get self-Compton
comptonFluxDensity = comptonSpec(log_ν, mps)

# Plot spectrum
ix_low = mps.ν_n
ix_high = 1
logFluxDensity = zeros(mps.ν_n)
logνFluxDensity = zeros(mps.ν_n)
for ix in range(1, mps.ν_n)
    if (syncFluxDensity[ix] > 0.0 || comptonFluxDensity[ix] > 0.0)
        ϵ = mps.h*10.0^log_ν[ix]/(mps.m_e*mps.c^2)
        logFluxDensity[ix] = log10(syncFluxDensity[ix] + comptonFluxDensity[ix])
        logνFluxDensity[ix] = log10(ϵ*(syncFluxDensity[ix] + comptonFluxDensity[ix]))
        if ix < ix_low
            global ix_low = ix
        end
        if ix > ix_high
            global ix_high = ix
        end
    else
        logFluxDensity[ix] = -50.0 # so the graph looks nice
        logνFluxDensity[ix] = -50.0
    end
end
print(logνFluxDensity)
print("\n", logFluxDensity)

q = plot(log_ν, logFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)
q
p = plot(log_ν, logνFluxDensity, ylim=(-16,-12), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)
p

print("Done!\n")