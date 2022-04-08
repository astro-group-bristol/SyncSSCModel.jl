using DiscJetConnections
using Parameters

# Set up parameters, changing mass to 10^7 M_Sun
mps = MyParamStruct(M8=0.1)

# Test Synchrotron plot
syncPlot(mps)
