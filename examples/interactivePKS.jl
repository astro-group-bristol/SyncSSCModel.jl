### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4a510b9e-213a-452e-b14b-3a8a12383c19
using Pkg; Pkg.activate("../.")

# ╔═╡ d2304800-21d3-495c-ad26-0f994e7328e9
begin
	using SyncSSCModel
	using SpectralFitting
	using Plots
	gr()
	using PlutoUI
end

# ╔═╡ cbf7004c-e681-464d-ba56-677e312bf3b5
md"""
# Testing interactive models in Pluto.jl
"""

# ╔═╡ 14dff0ad-3a72-48ab-ab62-a1386a8530f5
md"""
Interactive exploration of the SSC model for PKS 0637-752 WK7.8
"""

# ╔═╡ bd0e25a0-b015-40e3-942f-219286a6ef90
md"""
Define SSC model
"""

# ╔═╡ 0ac10751-444f-423b-a4e7-5839e29be045
function PKSSSCModel(;
    K = FitParam(1.0, lower_limit = 0.5, upper_limit = 2.0),
    log_B = FitParam(log10(1.25E-5), lower_limit = -6.0, upper_limit = -3.0),
    n_e0 = FitParam(19.0, lower_limit = 0.0, upper_limit = 1.0E2),
    log_radius = FitParam(log10(1.0E22), lower_limit = 21.0, upper_limit = 23.0),
    Γ = FitParam(2.0),
    log_γ_min = FitParam(log10(2.5E3), lower_limit = 1.0),
    log_γ_max = FitParam(log10(4.0E6), upper_limit = 8.0),
    p = FitParam(2.6, lower_limit = 2.0, upper_limit = 3.5),
    log_dL = FitParam(log10(1.26E28), lower_limit=23.0, upper_limit=29.0),
    θ = FitParam(60.0 * pi / 180.0),
    z = FitParam(0.651),
)
    SSCModel{
        typeof(K),
        SpectralFitting.FreeParameters{(:log_B, :n_e0, :θ)},
    }(
        K, log_B, n_e0, p, log_γ_min, log_γ_max, Γ, log_radius, θ, log_dL, z,
    )
end

# ╔═╡ 7c509230-b97e-4ed3-b46e-8151625096d6
md"""
Define PKS 0637-752 dataset of flux densities at different frequencies
"""

# ╔═╡ ec2cfccf-c976-4363-8de4-dac4b4f442ae
begin
    lines = readlines(@__DIR__() * "/PKS0637-752.txt")
    number_expr = r"(-?\d*\.\d*)"
    search_expr = r"^" * number_expr * r"\s+" * number_expr
    data_stacked = map(filter(!isnothing, match.(search_expr, lines))) do m
        parse.(Float64, m.captures)
    end
    # sort just for coherence
    sort!(data_stacked)
    # flatten array
    data = reduce(hcat, data_stacked)
end

# ╔═╡ a980821a-e07d-4983-a60c-1a47874f52db
dataset = SimpleDataset(
    "PKS0637data",
    10 .^ data[1, :],
    10 .^ data[2, :],
    x_units = SpectralFitting.SpectralUnits.u"Hz",
    x_err = 0.05 .* (10 .^ data[1, :]),
    y_err = 0.1 .* (10 .^ data[2, :]),
)

# ╔═╡ 53a26f8a-89de-4e1d-bfc1-9a0ea676bcbb
md"""
Adjustable model parameters
"""

# ╔═╡ 2b095fe7-f16d-4d30-a0db-764e087e897f
@bind var_n_e0 Slider(1:0.1:20.0, default=5.0)

# ╔═╡ 2728a610-48dd-40a0-9f19-306d9b33e80f
@bind var_log_B Slider(-7.0:0.01:-4.0, default=-5.13)

# ╔═╡ f3b1babf-c7ad-4f66-91be-4b14b3c24ba8
@bind var_θ Slider(1:1:90, default=74)

# ╔═╡ a022f9ca-9a3d-46c4-af36-053a1ecd3588
print("Density n_e0 = ", var_n_e0, " Magnetic field B = ", var_log_B, " θ = ", var_θ, " degrees")

# ╔═╡ 1175735a-22f1-4662-8029-8983ad747318
md"""
Evaluate model, plot model, overplot data
"""

# ╔═╡ 7efde37e-52ee-4239-846d-08b33fc9d5fb
begin
	νrange =  10 .^ collect(range(7, 26, 100))
	model = PKSSSCModel()
	model_free_pars = [var_log_B, var_n_e0, var_θ*π/180.0]
	flux = invokemodel(νrange, model, model_free_pars)
	nonzero = findall(x -> x > 0.0, flux)
	plot(νrange[nonzero], flux[nonzero], xrange=(1e7, 1e26), yrange=(1E-17, 5E-13), xscale=:log10, yscale=:log10, xlabel="ν (Hz)", ylabel="Flux (units)", label = "Model", legend = :topleft)
	plot!(dataset.x, dataset.y, seriestype = :scatter, markersize = 3, markerstrokewidth = 0, mc=:red, label = "Data")
end

# ╔═╡ Cell order:
# ╟─cbf7004c-e681-464d-ba56-677e312bf3b5
# ╟─14dff0ad-3a72-48ab-ab62-a1386a8530f5
# ╠═4a510b9e-213a-452e-b14b-3a8a12383c19
# ╠═d2304800-21d3-495c-ad26-0f994e7328e9
# ╟─bd0e25a0-b015-40e3-942f-219286a6ef90
# ╠═0ac10751-444f-423b-a4e7-5839e29be045
# ╟─7c509230-b97e-4ed3-b46e-8151625096d6
# ╠═ec2cfccf-c976-4363-8de4-dac4b4f442ae
# ╠═a980821a-e07d-4983-a60c-1a47874f52db
# ╟─53a26f8a-89de-4e1d-bfc1-9a0ea676bcbb
# ╠═2b095fe7-f16d-4d30-a0db-764e087e897f
# ╠═2728a610-48dd-40a0-9f19-306d9b33e80f
# ╠═f3b1babf-c7ad-4f66-91be-4b14b3c24ba8
# ╠═a022f9ca-9a3d-46c4-af36-053a1ecd3588
# ╟─1175735a-22f1-4662-8029-8983ad747318
# ╠═7efde37e-52ee-4239-846d-08b33fc9d5fb
