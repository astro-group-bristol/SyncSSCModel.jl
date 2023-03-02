push!(LOAD_PATH,"src")

using Documenter
using SyncSSCModel

makedocs(
    modules=[SyncSSCModel],
    clean=false,
    sitename="Documentation",

    pages = [
        "Home" => "index.md"
    ]
)

# Update the following to have the correct user name
# for your local version or the astro-group-bristol version
deploydocs(
# repo = "github.com/phajy/SyncSSCModel.jl.git"
repo = "github.com/astro-group-bristol/SyncSSCModel.jl.git"
)
