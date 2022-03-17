push!(LOAD_PATH,"src")

using Documenter
using DiscJetConnections

makedocs(sitename="Synchrotron Documentation")
makedocs(
    modules=[DiscJetConnections],
    clean=true,
    sitename="DiscJetConnections.jl"
)
