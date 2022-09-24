using DiscJetConnections
using Parameters

# Set up parameters, changing mass to 10^7 M_Sun
mps = MyParamStruct(M8=1.0E8)

# Test Synchrotron plot
syncPlot(mps)

# Test Synchrotron Self-Compton plot
ssc_Plot(mps)

# Test General Synchrotron plot
Sync_Plot(mps)

print("Done!\n")
