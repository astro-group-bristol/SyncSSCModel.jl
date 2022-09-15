using DiscJetConnections
using Parameters
using Plots
gr()

# Define a frequency grid and corresponding spectral flux density
mps = MyParamStruct(B=1.25E-6, radius=1.0E22, Γ=1.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=6.0*pi/180.0, n_e0=19.0)
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
        logFluxDensity[ix] = log10(syncFluxDensity[ix] + comptonFluxDensity[ix])
        logνFluxDensity[ix] = log10(10.0^log_ν[ix]*(syncFluxDensity[ix] + comptonFluxDensity[ix]))
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

# p = plot(log_ν, logFluxDensity, ylim=(-30,-2))
p = plot(log_ν, logνFluxDensity, ylim=(2,10), xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)")
p

print("Done!\n")