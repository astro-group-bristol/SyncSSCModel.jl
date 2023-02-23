using DiscJetConnections
using Parameters
using Statistics
using LaTeXStrings
using Plots
gr()
print("START")

# PARAMETERS FROM PAPERS WK7.8 ###############################################################
# Tavecchio et al. 2000 (SSC and EC/CMB)
mpstav_SSC = MyParamStruct(B=1.25E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0, Ho=50)
mpstav_EC = MyParamStruct(B=1.5E-5, radius=1.0E22, Γ=10.0, γ_min=10, γ_max=4.0E5, p=2.6, dL=1.26E28, θ=6.0*pi/180.0, n_e0=6.0E-5, Ho=50)
# Harris and Krawczynski 2002
mpshk = MyParamStruct(B=195.0E-6, radius=1.0E22, Γ=9.9, γ_min=45, γ_max=1.4E6, p=2.62, dL=1.26E28, θ=5.7*pi/180.0, n_e0=6.0E-5, Ho=50)
# Uchiyama et al. 2005 (Synchrotron, IC+Sync) Default n_e0=6.0E-5
mpsu1 = MyParamStruct(B=1.5E-5, radius=3.086E21, Γ=12.0, γ_min=1.78E3, γ_max=2.35E5, p=2.6, dL=1.26E28, θ=3.5*pi/180.0, n_e0=1.0E-3, Ho=71)
mpsu2 = MyParamStruct(B=1.5E-5, radius=3.086E21, Γ=12.0, γ_min=20, γ_max=9.8E5, p=2.6, dL=1.26E28, θ=6.5*pi/180.0, n_e0=1.0E-3, Ho=71)
# Lucchini et al. 2017 (EC/CMB)
mpslu = MyParamStruct(B=14.5E-6, radius=1.0E22, Γ=10.0, γ_min=20, γ_max=1.1E6, p=2.65, dL=1.26E28, θ=5.9*pi/180.0, n_e0=44.0E-6, Ho=70)
# Celotti et al. 2001 (SSC) TO VERIFY
mpscel1 = MyParamStruct(B=2.40E-5, radius=9.26E21, Γ=1.1, γ_min=30.0, γ_max=1.0E6, p=2.5, dL=1.26E28, θ=60.0*pi/180.0, n_e0=2.07E-3, Ho=50, z=0.654)
mpscel2 = MyParamStruct(B=2.4E-5, radius=9.257E21, Γ=14.0, γ_min=30, γ_max=1.0E6, p=2.5, dL=1.26E28, θ=5.0*pi/180.0, n_e0=6.0E-5, Ho=50)
# Atoyan & Dermer 2004  TO VERIFY
mpsad1 = MyParamStruct(B=1E-4, radius=3.09E21, Γ=5.0, γ_min=20.0, γ_max=1.8E5, p=2.6, dL=1.20E28, θ=10.0*pi/180.0, n_e0=1.0E-3, Ho=70)
mpsad2 = MyParamStruct(B=1.0E-4, radius=3.09E21, Γ=5.0, γ_min=2.0E7, γ_max=1.0E11, p=2.2, dL=1.20E28, θ=10.0*pi/180.0, n_e0=3.98E-3, Ho=70)
###############################################################

# TESTING DIFFERENT VALUE OF θ and B IN Tavecchio+00 parameter ###############################################################
# Small pitches angles θ=1°
mps1a1 = MyParamStruct(B=1E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=1.0*pi/180.0, n_e0=19.0, Ho=50)
mps1b1 = MyParamStruct(B=1E-5, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=1.0*pi/180.0, n_e0=19.0, Ho=50)
mps1c1 = MyParamStruct(B=1E-4, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=1.0*pi/180.0, n_e0=19.0, Ho=50)
mps1d1 = MyParamStruct(B=1E-3, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=1.0*pi/180.0, n_e0=19.0, Ho=50)
mps1e1 = MyParamStruct(B=1E-2, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=1.0*pi/180.0, n_e0=19.0, Ho=50)
# Small pitches angles θ=5°
mps1a5 = MyParamStruct(B=1E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=5.0*pi/180.0, n_e0=19.0, Ho=50)
mps1b5 = MyParamStruct(B=1E-5, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=5.0*pi/180.0, n_e0=19.0, Ho=50)
mps1c5 = MyParamStruct(B=1E-4, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=5.0*pi/180.0, n_e0=19.0, Ho=50)
mps1d5 = MyParamStruct(B=1E-3, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=5.0*pi/180.0, n_e0=19.0, Ho=50)
mps1e5 = MyParamStruct(B=1E-2, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=5.0*pi/180.0, n_e0=19.0, Ho=50)
# high pitches angles θ=10°
mps1a10 = MyParamStruct(B=1E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=10.0*pi/180.0, n_e0=19.0, Ho=50)
mps1b10 = MyParamStruct(B=1E-5, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=10.0*pi/180.0, n_e0=19.0, Ho=50)
mps1c10 = MyParamStruct(B=1E-4, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=10.0*pi/180.0, n_e0=19.0, Ho=50)
mps1d10 = MyParamStruct(B=1E-3, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=10.0*pi/180.0, n_e0=19.0, Ho=50)
mps1e10 = MyParamStruct(B=1E-2, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=10.0*pi/180.0, n_e0=19.0, Ho=50)
# high pitches angles θ=30°
mps1a30 = MyParamStruct(B=1E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=30.0*pi/180.0, n_e0=19.0, Ho=50)
mps1b30 = MyParamStruct(B=1E-5, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=30.0*pi/180.0, n_e0=19.0, Ho=50)
mps1c30 = MyParamStruct(B=1E-4, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=30.0*pi/180.0, n_e0=19.0, Ho=50)
mps1d30 = MyParamStruct(B=1E-3, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=30.0*pi/180.0, n_e0=19.0, Ho=50)
mps1e30 = MyParamStruct(B=1E-2, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=30.0*pi/180.0, n_e0=19.0, Ho=50)
# high pitches angles θ=45°
mps1a45 = MyParamStruct(B=1E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=45.0*pi/180.0, n_e0=19.0, Ho=50)
mps1b45 = MyParamStruct(B=1E-5, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=45.0*pi/180.0, n_e0=19.0, Ho=50)
mps1c45 = MyParamStruct(B=1E-4, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=45.0*pi/180.0, n_e0=19.0, Ho=50)
mps1d45 = MyParamStruct(B=1E-3, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=45.0*pi/180.0, n_e0=19.0, Ho=50)
mps1e45 = MyParamStruct(B=1E-2, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=45.0*pi/180.0, n_e0=19.0, Ho=50)
# high pitches angles θ=60°
mps1a60 = MyParamStruct(B=1E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0, Ho=50)
mps1b60 = MyParamStruct(B=1E-5, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0, Ho=50)
mps1c60 = MyParamStruct(B=1E-4, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0, Ho=50)
mps1d60 = MyParamStruct(B=1E-3, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0, Ho=50)
mps1e60 = MyParamStruct(B=1E-2, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0, Ho=50)
# high pitches angles θ=80°
mps1a80 = MyParamStruct(B=1E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=80.0*pi/180.0, n_e0=19.0, Ho=50)
mps1b80 = MyParamStruct(B=1E-5, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=80.0*pi/180.0, n_e0=19.0, Ho=50)
mps1c80 = MyParamStruct(B=1E-4, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=80.0*pi/180.0, n_e0=19.0, Ho=50)
mps1d80 = MyParamStruct(B=1E-3, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=80.0*pi/180.0, n_e0=19.0, Ho=50)
mps1e80 = MyParamStruct(B=1E-2, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=80.0*pi/180.0, n_e0=19.0, Ho=50)
# high pitches angles θ=90°
mps1a90 = MyParamStruct(B=1E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=90.0*pi/180.0, n_e0=19.0, Ho=50)
mps1b90 = MyParamStruct(B=1E-5, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=90.0*pi/180.0, n_e0=19.0, Ho=50)
mps1c90 = MyParamStruct(B=1E-4, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=90.0*pi/180.0, n_e0=19.0, Ho=50)
mps1d90 = MyParamStruct(B=1E-3, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=90.0*pi/180.0, n_e0=19.0, Ho=50)
mps1e90 = MyParamStruct(B=1E-2, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=90.0*pi/180.0, n_e0=19.0, Ho=50)
# mps = [mps1a1,mps1b1,mps1c1,mps1d1,mps1e1,mps1a5,mps1b5,mps1c5,mps1d5,mps1e5,mps1a10,mps1b10,mps1c10,mps1d10,mps1e10,mps1a30,mps1b30,mps1c30,mps1d30,mps1e30,mps1a45,mps1b45,mps1c45,mps1d45,mps1e45,mps1a60,mps1b60,mps1c60,mps1d60,mps1e60,mps1a80,mps1b80,mps1c80,mps1d80,mps1e80,mps1a90,mps1b90,mps1c90,mps1d90,mps1e90]
###############################################################

# PARAMETERS WH Pic A from Wilson et al. 2001 Table 4 ###############################################################
mpswi1 = MyParamStruct(B=3.3E-5, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=1.0E6, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=5.3, z=0.035, Ho=70)
mpswi2 = MyParamStruct(B=5.3E-5, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=4.75E5, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=2.2, z=0.035, Ho=70)
mpswi3 = MyParamStruct(B=2.8E-6, radius=7.7E20, Γ=1.0, γ_min=1.0E2, γ_max=9.0E3, p=3.30, dL=4.752E26, θ=23.0*pi/180.0, n_e0=485.0, z=0.035, Ho=70)
mpswi4a = MyParamStruct(B=4.7E-4, radius=7.7E20, Γ=1.0, γ_min=1.0E2, γ_max=5.0E4, p=2.14, dL=4.752E26, θ=23.0*pi/180.0, n_e0=1.0E-3, z=0.035, Ho=70)
mpswi4b = MyParamStruct(B=4.7E-4, radius=7.7E20, Γ=1.0, γ_min=1.0E2, γ_max=4.0E6, p=2.14, dL=4.752E26, θ=23.0*pi/180.0, n_e0=1.2E-5, z=0.035, Ho=70)
# mps = [mpswi1,mpswi2]
###############################################################

# WH PICTOR A testing different value of n_e0 ###############################################################
mpswi1a = MyParamStruct(B=3.3E-5, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=1.0E6, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=500.0, z=0.035, Ho=70)
mpswi1b = MyParamStruct(B=3.3E-5, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=1.0E6, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=5.3, z=0.035, Ho=70)
mpswi1c = MyParamStruct(B=3.3E-5, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=1.0E6, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=1.0E-5, z=0.035, Ho=70)
# mps = [mpswi1a,mpswi1b,mpswi1c]
mpswi2a = MyParamStruct(B=53.0E-6, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=4.75E5, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=500.0, z=0.035, Ho=70)
mpswi2b = MyParamStruct(B=53.0E-6, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=4.75E5, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=2.2, z=0.035, Ho=70)
mpswi2c = MyParamStruct(B=53.0E-6, radius=7.7E20, Γ=1.0, γ_min=8.7E1, γ_max=4.75E5, p=2.48, dL=4.752E26, θ=23.0*pi/180.0, n_e0=1.0E-5, z=0.035, Ho=70)
# mps = [mpswi2a,mpswi2b,mpswi2c]
###############################################################

# PKS0637-752 WK7.8 testing different value of n_e0 ###############################################################
mpsWK1 = MyParamStruct(B=1.25E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=500.0, Ho=50)
mpsWK2 = MyParamStruct(B=1.25E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=19.0, Ho=50)
mpsWK3 = MyParamStruct(B=1.25E-6, radius=1.0E22, Γ=2.0, γ_min=2.5E3, γ_max=4.0E6, p=2.6, dL=1.26E28, θ=60.0*pi/180.0, n_e0=1.0E-5, Ho=50)
# mps = [mpsWK1,mpsWK2,mpsWK3]
###############################################################

# SIMULATED DATA WITH RANDOM VALUE OF PARAMETERS with different value of θ ###############################################################
mpssim1 = MyParamStruct(B=1.0E-6, radius=1.0E22, Γ=2.0, γ_min=1.0E3, γ_max=1.0E6, p=2.5, dL=1E28, θ=5.0*pi/180.0, n_e0=1.0, Ho=50)
mpssim2 = MyParamStruct(B=1.0E-6, radius=1.0E22, Γ=2.0, γ_min=1.0E3, γ_max=1.0E6, p=2.5, dL=1E28, θ=10.0*pi/180.0, n_e0=1.0, Ho=50)
mpssim3 = MyParamStruct(B=1.0E-6, radius=1.0E22, Γ=2.0, γ_min=1.0E3, γ_max=1.0E6, p=2.5, dL=1E28, θ=30.0*pi/180.0, n_e0=1.0, Ho=50)
mpssim4 = MyParamStruct(B=1.0E-6, radius=1.0E22, Γ=2.0, γ_min=1.0E3, γ_max=1.0E6, p=2.5, dL=1E28, θ=80.0*pi/180.0, n_e0=1.0, Ho=50)
# mps = [mpssim1,mpssim2,mpssim3,mpssim4]
###############################################################

# Set an array of mps ###############################################################
mps = [mpstav_SSC] # Tavecchio et al. 2000 / Table 3.2 [1] in thesis
# mps = [mpscel1] # Celotti et al. 2001 / Table 3.2 [2] in thesis
# mps = [mpswi1,mpswi2] # Wilson et al. 2001 / Table 3.3 [3a] and [3b] in thesis
###############################################################

# MAIN PROGRAMME SYN+SSC MODEL ###############################################################
# Define a frequency grid and corresponding spectral flux density
# Create array
all_logFluxDensity = zeros(100, length(mps))
all_logνFluxDensity = zeros(100, length(mps))
all_ρ_ssc_syn = zeros(length(mps))
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
            global ϵ = mps[x].h*10.0^log_ν[ix]/(mps[x].m_e*mps[x].c^2)
            # logFluxDensity[ix] = log10(syncFluxDensity[ix] + comptonFluxDensity[ix])
            # logνFluxDensity[ix] = log10(ϵ*(syncFluxDensity[ix] + comptonFluxDensity[ix]))
            all_logFluxDensity[ix,x] = log10(syncFluxDensity[ix] + comptonFluxDensity[ix])
            all_logνFluxDensity[ix,x] = log10(ϵ*((syncFluxDensity[ix]) + (comptonFluxDensity[ix])))
            # all_logνFluxDensity[ix,x] = log10(ϵ*((2.09*syncFluxDensity[ix]) + (0.37*comptonFluxDensity[ix])))
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
    # all_ρ_ssc_syn[x] = log10(mean(comptonFluxDensity)) / log10(mean(syncFluxDensity))
    all_ρ_ssc_syn[x] = mean(comptonFluxDensity) / mean(syncFluxDensity)
end    
print("νFν = ", all_logνFluxDensity)
print("\n Fν = ", all_logFluxDensity)
print("\n ρ_ssc_syn = ", all_ρ_ssc_syn)
###############################################################

# Plot Ratio SSC Sync ρ(B,θ) ############################################################### 
# θ_test = [0, 5, 10, 15, 20, 25, 30, 45, 60, 90]
# p0_ρ = plot(θ_test, all_ρ_ssc_syn, yaxis=:log10)
###################### UNCOMMENT ME / Fig 3.4 in thesis
### see "TESTING DIFFERENT VALUE OF θ and B IN Tavecchio+00 parameter, mps = [mps1a1,mps1b1,mps1c1,mps1d1,mps1e1,mps1a5,mps1b5,mps1c5,mps1d5,mps1e5,mps1a10,mps1b10,mps1c10,mps1d10,mps1e10,mps1a30,mps1b30,mps1c30,mps1d30,mps1e30,mps1a45,mps1b45,mps1c45,mps1d45,mps1e45,mps1a60,mps1b60,mps1c60,mps1d60,mps1e60,mps1a80,mps1b80,mps1c80,mps1d80,mps1e80,mps1a90,mps1b90,mps1c90,mps1d90,mps1e90]"
# Bi_test = [1E-6,1E-5,1E-4,1E-3,1E-2]
# ρ_ssc_syn_1 = all_ρ_ssc_syn[1:5]
# ρ_ssc_syn_5 = all_ρ_ssc_syn[6:10]
# ρ_ssc_syn_10 = all_ρ_ssc_syn[11:15]
# ρ_ssc_syn_30 = all_ρ_ssc_syn[16:20]
# ρ_ssc_syn_45 = all_ρ_ssc_syn[21:25]
# ρ_ssc_syn_60 = all_ρ_ssc_syn[26:30]
# ρ_ssc_syn_80 = all_ρ_ssc_syn[31:35]
# ρ_ssc_syn_90 = all_ρ_ssc_syn[36:40]
# y_limit1 = 10^(-6.2)
# y_limit2 = 10^(-6.5)
# ρ_ssc_syn_deg = [ρ_ssc_syn_1,ρ_ssc_syn_5,ρ_ssc_syn_10,ρ_ssc_syn_30,ρ_ssc_syn_45,ρ_ssc_syn_60,ρ_ssc_syn_80,ρ_ssc_syn_90]
# p0_Bi = plot(Bi_test, ρ_ssc_syn_deg,xaxis=:log10, yaxis=:log10, ylim=(y_limit2,y_limit1), xlabel=L"B[G]",ylabel=L"ρ_{ssc/syn}",xminorticks=5, yminorticks=5,title = "Variation of the power ratio ρ",titlefontsize=11,label=["θ=1°" "θ=5°" "θ=10°" "θ=30°" "θ=45°" "θ=60°" "θ=80°" "θ=90°"], framestyle=:box, legend=:topright, size=(750,750), xguidefontsize=12, yguidefontsize=12, legendfontsize=11)
# p0_Bi
###############################################################

# n_e0 vs B ############################################################### UNCOMMENT ME FOR ONLY ONE VALUE OF mps / Fig 3.5 in thesis
# function Σ_c_test(B_test,δ_test,mean_ϵ,z_test,γ_min_test,γ_max_test)
#     ϵ_B_test = B_test/(4.414E13)
#     A = δ_test/(mean_ϵ*(1+z_test))
#     B = (mean_ϵ*(1+z_test))/(δ_test*γ_min_test^2)
#     C = γ_max_test^2*ϵ_B_test
#     min_test = min(A,B,C)
#     D = γ_min_test^2*ϵ_B_test
#     E = (mean_ϵ*(1+z_test))/(δ_test*γ_max_test^2)
#     max_test = max(D,E)
#     Σ_c_test = min_test/max_test
#     if Σ_c_test > 0
#     return(log(Σ_c_test))
#     else
#     return(0.0)
#     end
#     end
# function n_e0_test(B_test)
#     n_e0_test = (3*mean_Pssc)/(2*σ_T_test*R_b_test*mean_Psyn*Σ_c_test(B_test,δ_test,mean_ϵ,z_test,γ_min_test,γ_max_test))
#     return(n_e0_test)
#     end
# obj = ["tav_SSC"] ### DON'T FORGET TO CHANGE THIS TO THE APPROPRIATE NAME E.G. mps"obj" = mps"sim1"  or  mps"obj" = mps"tav" ...
# for x in eachindex(obj)
#     global β_test = sqrt(1.0 - 1.0/mps[x].Γ^2)
#     global μ_obs_test = cos(mps[x].θ)
#     global δ_test = 1.0 / (mps[x].Γ*(1 - μ_obs_test*β_test))
#     global σ_T_test = mps[x].σ_T
#     global R_b_test = mps[x].radius
#     global γ_max_test = mps[x].γ_max
#     global γ_min_test = mps[x].γ_min
#     global z_test = mps[x].z
#     global h_test = mps[x].h
#     global c_test = mps[x].c
#     global m_e_test = mps[x].m_e
# end
# log_ν_test = log_ν
# ϵ_test = zeros(length(log_ν_test))
# for x in eachindex(log_ν_test)
#     ϵ_test[x]=(h_test*10^(log_ν_test[x]))/(m_e_test*c_test^2)
# end
# mean_ϵ = mean(ϵ_test)
# syncFluxDensity_test = zeros(length(log_ν_test))
# comptonFluxDensity_test = zeros(length(log_ν_test))
# for x in eachindex(ϵ_test)
#     comptonFluxDensity_test[x]=ϵ_test[x]*comptonFluxDensity[x]
#     syncFluxDensity_test[x]=ϵ_test[x]*syncFluxDensity[x]
# end
# mean_Psyn = mean(syncFluxDensity_test)
# mean_Pssc = mean(comptonFluxDensity_test)
# B_test = [1.0E-6, 1.0E-5, 1.0E-4, 1.0E-3, 1.0E-2]
# n_e0_value_test = zeros(length(B_test))
# Σ_c_test_value = zeros(length(B_test))
# for x in eachindex(B_test)
#     Σ_c_test_value[x] = Σ_c_test(B_test[x],δ_test,mean_ϵ,z_test,γ_min_test,γ_max_test)
# end
# for x in range(1,5)
#     n_e0_value_test[x] = (3)/(2*σ_T_test*R_b_test*Σ_c_test_value[x])
#     # n_e0_value_test[x] = (3*mean_Pssc)/(2*σ_T_test*R_b_test*mean_Psyn*Σ_c_test_value[x])
# end
# p0_ne0_test = plot(B_test[1:5],n_e0_value_test[1:5],xlabel=L"\mathrm{B\, (G)}",ylabel=L"\mathrm{n_{e0}\, (cm^{-3})}",xaxis=:log10,xminorticks=10,yminorticks=5,framestyle=:box,label=L"n_{e0} \propto \frac{3}{2\sigma_T R_b \ln \Sigma_c (B)}",legendfontsize=16,size=(500,500))
###############################################################

# Plot νFν PKS0637-752 WK7.8 ###############################################################
# all_p0_WK78 = plot(log_ν, all_logνFluxDensity, ylim=(-22,-10), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["Tavecchio+00" "Harris+02" "Uchiyama+05 (IC+Syn)" "Lucchini+17" "Celotti+01"], framestyle=:box)
p0_WK78 = plot(log_ν, all_logνFluxDensity, ylim=(-50,5), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["Tavecchio+00" "Celotti+01" "Atoyan+04"], framestyle=:box)
p0_WK78_tav = plot(log_ν, all_logνFluxDensity, ylim=(-25,-12), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 12, label = "Tavecchio+00", framestyle=:box, legend=:topleft, size=(750,750), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)
p0_WK78_cel = plot(log_ν, all_logνFluxDensity, ylim=(-19,-14), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 12, label = "Celotti+01", framestyle=:box, legend=:topright, size=(500,350), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)
p0_WK78_ato = plot(log_ν, all_logνFluxDensity, ylim=(-20,-12), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 12, label = "Atoyan+04", framestyle=:box, legend=:topright, size=(500,350), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)
p0_WK78_uch = plot(log_ν, all_logνFluxDensity, ylim=(-20,-12), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 12, label = ["Sync" "IC+Sync" "Uchiyama+05"], framestyle=:box, legend=:topright, size=(500,350), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)
p0_WK78_lu = plot(log_ν, all_logνFluxDensity, ylim=(-20,-12), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 12, label = "Lucchini+17", framestyle=:box, legend=:topright, size=(500,350), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)

# Plotting the data fluxes from observation / Fig 3.1 in thesis
ν_Schwartz = [log10(2.42*10^17)] # log10[Hz]
νFν_Schwartz = [log10(1.52*10^(-14))] 
ν_Lovell = [log10(8.60*10^9)] # log10[Hz]
νFν_Lovell = [log10(3.78*10^(-15))]
ν_Uchiyama = [log10(5.17*10^(13)),log10(8.33*10^(13))] # log10[Hz]
νFν_Uchiyama = [log10(2.38*10^(-15)),log10(2.25*10^(-15))]
ν_Mehta = [log10(1.87*10^(14)),log10(4.3*10^(14)),log10(6.32*10^(14))] # log10[Hz]
νFν_Mehta = [log10(1.83*10^(-15)),log10(8.61*10^(-16)),log10(7.54*10^(-16))]
all_p0_WK78_tav = plot(log_ν, all_logνFluxDensity, ylim=(-16,-12), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux WK7.8 PKS0637-752", titlefontsize = 12, label = "Tavecchio+00", framestyle=:box, legend=:topleft, size=(500,500), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12) 
npoa = scatter!(ν_Schwartz,νFν_Schwartz, label="Schwartz+00", markershape=:diamond)
npob = scatter!(ν_Lovell,νFν_Lovell, label="Lovell+00", markershape=:dtriangle)
npoc = scatter!(ν_Uchiyama,νFν_Uchiyama, label="Uchiyama+05",markershape=:circle)
npod = scatter!(ν_Mehta,νFν_Mehta, label="Mehta+09", markershape=:+)

# Plot Fν PKS0637-752 WK7.8
# q0_WK78 = plot(log_ν, all_logFluxDensity, ylim=(-25,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = ["Tavecchio+00" "Celotti+01"], framestyle=:box)
q0_WK78_tav = plot(log_ν, all_logFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = "Tavecchio+00", framestyle=:box, legend=:topright, size=(500,500), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)
###############################################################

# Plot νFν WH Pic A ###############################################################
# Plotting the data fluxes from observation / Fig 3.3 in thesis
ν_Roser = [17.7] # log10[Hz]
νFν_Roser = [-12.35] 
ν_Meisen = [9.7] # log10[Hz]
νFν_Meisen = [-12.96]
ν_Tingay = [13.10,13.92] # log10[Hz]
νFν_Tingay = [-11.98,-11.9]
# all_p0_PicA_wi = plot(log_ν, all_logνFluxDensity, ylim=(-20,-10), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 12, label = ["model 1" "model 2" "model 3" "model 4a" "model 4b"], framestyle=:box, legend=:topleft, size=(750,750), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)
all_p0_PicA_wi = plot(log_ν, all_logνFluxDensity, ylim=(-14,-10), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux WH Pic A", titlefontsize = 12, label = ["Wilson+01 model 1" "Wilson+01 model 2"], framestyle=:box, legend=:topleft, size=(500,500), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=11) #, linestyle = [:dash :solid]
npo = scatter!(ν_Roser,νFν_Roser, label="Roser+87", mode="markers",sizemode="area") # yerror = 0.1,
npo1 = scatter!(ν_Meisen,νFν_Meisen, yerror=0.9E-2, label="Meisenheimer+97", mode="markers")
npo2 = scatter!(ν_Tingay,νFν_Tingay, yerror=[0.1, 0.05], label="Tingay+08")
p0_PicA_wi = plot(log_ν, all_logνFluxDensity, ylim=(-14,-10), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 12, label = "Wilson+01", framestyle=:box, legend=:topleft, size=(500,500), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)
q0_PicA_wi = plot(log_ν, all_logFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = "Wilson+01", framestyle=:box, legend=:topright, size=(500,500), xtickfontsize=12, ytickfontsize=12, xguidefontsize=12, yguidefontsize=12, legendfontsize=12)
###############################################################

# OTHER PLOT ############################################################### IF NEEDED OR IF YOU JUST WANT TO TEST SOME OTHER PLOTS
# Plot Test Γ, p, ne0
# p0_Γ = plot(log_ν, all_logνFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["Γ=1" "Γ=5" "Γ=10"], framestyle=:box)#, linecolor = "blue2", line = :dash)
# p0_p = plot(log_ν, all_logνFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["p=2" "p=2.5" "p=3"], framestyle=:box)#, linecolor = "blue2", line = :dash)
p0_ne0_picA = plot(log_ν, all_logνFluxDensity, ylim=(-30,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux WH Pic A", titlefontsize = 11, label = ["n_e0=500" "n_e0=5.3" "n_e0=1E-5"], framestyle=:box)#, legend=:topright, size=(500,500), xtickfontsize=11, ytickfontsize=11, xguidefontsize=11, yguidefontsize=11, legendfontsize=11linecolor = "blue2", line = :dash)
p0_ne0_WK = plot(log_ν, all_logνFluxDensity, ylim=(-30,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux WK7.8 PKS0637-752", titlefontsize = 11, label = ["n_e0=500" "n_e0=19.0" "n_e0=1E-5"], framestyle=:box)#, legend=:topright, size=(500,500), xtickfontsize=11, ytickfontsize=11, xguidefontsize=11, yguidefontsize=11, legendfontsize=11linecolor = "blue2", line = :dash)
# p0_test = plot(log_ν, all_logνFluxDensity, ylim=(-30,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = ["PKS0637-752 p=2.6" "Pictor A p=2.49" "Pictor A p=3.3"], framestyle=:box)#, linecolor = "blue2", line = :dash)
##########
# PKS0637-752 OLD version
# q = plot(log_ν, logFluxDensity, ylim=(-20,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)
# p = plot(log_ν, logνFluxDensity, ylim=(-16,-12), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)
##########
# PictorA OLD version
# q = plot(log_ν, logFluxDensity, ylim=(-26,0), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log F(ν) (cgs)", title = "Flux Density", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)
# p = plot(log_ν, logνFluxDensity, ylim=(-14,-8), xlims=(7, 26), xminorticks=5, yminorticks=5, xlabel="log ν (Hz)", ylabel="log νF(ν) (cgs)", title = "Spectral Power Flux", titlefontsize = 11, label = false, framestyle=:box)#, linecolor = "blue2", line = :dash)
###############################################################
print("Done!\n")