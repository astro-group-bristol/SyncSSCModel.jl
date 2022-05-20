var documenterSearchIndex = {"docs":
[{"location":"#synchrotron.jl-Documentation","page":"Home","title":"synchrotron.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"dn_e(γ, mps)\nj_syn(ϵ, mps)\nS_syn(ϵ, mps)\nP_syn(ϵ, mps)\nsyncPlot(mps)","category":"page"},{"location":"#DiscJetConnections.dn_e-Tuple{Any, Any}","page":"Home","title":"DiscJetConnections.dn_e","text":"dn_e(γ, mps)\n\nDifferential electron density at Lorentz factor γ for parameters mps.\n\nThe number of electrons with Lorentz factors in the range gamma to gamma + dgamma is given by dn_e\n\ndn_e = n_e_0 gamma^-p dgamma\n\nwhere\n\ngamma_min le gamma le gamma_max\n\nsee also the equation 7.20 in Rybicki and Lightman (1979)\n\n\n\n\n\n","category":"method"},{"location":"#DiscJetConnections.j_syn-Tuple{Any, Any}","page":"Home","title":"DiscJetConnections.j_syn","text":"j_syn(ϵ, mps)\n\nSynchrotron emissivity at photon energy ϵ for parameters mps.\n\nFor an isotropic electron distribution in a randomly oriented magnetic field, the synchrotron emissivity is (see equation 13 of Dermer et al. (1997); we've changed their H to our B)\n\nj_syn(epsilon Omega x) = fracc sigma_T u_B6 pi epsilon_B left( fracepsilonepsilon_B right)^tiny12 smalln_e left left( fracepsilonepsilon_Bright)^tiny12x right\n\n\n\n\n\n","category":"method"},{"location":"#DiscJetConnections.S_syn-Tuple{Any, Any}","page":"Home","title":"DiscJetConnections.S_syn","text":"S_syn(ϵ, mps)\n\nSynchrotron Flux Density at observed photon energy ϵ.\n\nThis impliments equation 3 of Dermer et al. (1997).\n\nS_syn(epsilon Omega x) = fracD^3 (1+z) V_b j_syn(fracepsilon (1+z)D Omega x)d_L^2\n\n\n\n\n\n","category":"method"},{"location":"#DiscJetConnections.P_syn-Tuple{Any, Any}","page":"Home","title":"DiscJetConnections.P_syn","text":"P_syn(ϵ, mps)\n\nSynchrotron Spectral Power Flux at observed photon energy ϵ.\n\nThis is equivalent to \\nu F(\\nu) and presented in equation 23 of Dermer et al. (1997).\n\nP_syn(epsilon Omega x) = epsilon  S_syn\n\n\n\n\n\n","category":"method"},{"location":"#DiscJetConnections.syncPlot-Tuple{Any}","page":"Home","title":"DiscJetConnections.syncPlot","text":"syncPlot(mps)\n\nExample synchotron plot for parameters mps.\n\n\n\n\n\n","category":"method"}]
}
