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

# Update the following to have the correct user name
# for your local version or the astro-group-bristol version
deploydocs(
    repo = "github.com/phajy/DiscJetConnections.jl.git"
    # repo = "github.com/astro-group-bristol/DiscJetConnections.jl.git"
)
