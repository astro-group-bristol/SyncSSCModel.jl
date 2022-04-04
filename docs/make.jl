push!(LOAD_PATH,"src")

using Documenter
using DiscJetConnections

makedocs(
    modules=[DiscJetConnections],
    clean=false,
    sitename="Documentation",

    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "github.com/phajy/DiscJetConnections.jl.git",
)
