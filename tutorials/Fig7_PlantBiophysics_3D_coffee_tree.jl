### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 56757555-725c-49b8-8ad9-84d99052a0c8
begin
	import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="PlantBiophysics", version="0.9.0"),
		Pkg.PackageSpec(name="PlantGeom"),
		Pkg.PackageSpec(name="BenchmarkTools"),
		Pkg.PackageSpec(name="PlantSimEngine"),
		Pkg.PackageSpec(name="PlantMeteo"),
		Pkg.PackageSpec(name="Dates"),
		Pkg.PackageSpec(name="FLoops"),
		Pkg.PackageSpec(name="MultiScaleTreeGraph"),
    ])
	using PlantBiophysics
	using PlantGeom, PlantMeteo, PlantSimEngine
	using BenchmarkTools
	using Dates
	using FLoops
	using MultiScaleTreeGraph
	nothing
end

# ╔═╡ 65ec1544-c81f-11ed-0a29-174eeb8e92e0
md"""
# PlantBiophysics.jl 3D global tree simulation

This Pluto notebook presents the computation of Fig. 7 from the scientific article. It displays leaf temperature on a 3D coffee tree simulated by PlantBiophysics.jl. Non-Pluto Julia script is also available (see [here](https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/Fig7_PlantBiophysics_3D_coffee_tree_noPluto.jl)).

## Importing the dependencies:

Loading the Julia packages:

"""

# ╔═╡ 8c41fd8a-648c-4d50-b39b-d6087b264d95
md"""
## Reading data 
#### MTG file
"""

# ╔═╡ 552692ae-2f8e-4ae1-88a0-5565edabdf3f
mtg = read_opf(joinpath(dirname(dirname(pathof(PlantBiophysics))), "test", "inputs", "scene", "opf", "coffee.opf"))

# ╔═╡ 22f4ac9d-0181-476b-a3dd-277b2c5e31a2
md"""
#### Meteorological data
"""

# ╔═╡ 2a106a14-3ab9-4ebc-87bd-8c86453a3166
weather = PlantMeteo.read_weather(
    joinpath(dirname(dirname(pathof(PlantMeteo))), "test", "data", "meteo.csv"),
    :temperature => :T,
    :relativeHumidity => (x -> x ./ 100) => :Rh,
    :wind => :Wind,
    :atmosphereCO2_ppm => :Cₐ,
    date_format = DateFormat("yyyy/mm/dd")
)

# ╔═╡ c13c9dd7-3a8c-436a-9520-c638a76b7135
md"""
#### Models list
"""

# ╔═╡ 00877758-9349-4aa6-a2c8-ac2e50e150a0
begin
	file = joinpath(dirname(dirname(pathof(PlantBiophysics))), "test", "inputs", "models", "plant_coffee.yml")
	models = read_model(file)

	to_initialize(models)
end


# ╔═╡ e9293038-6179-481e-93c9-a619bcbce8ca
md"""
## Adding light interception data: incident radiation
"""

# ╔═╡ f02f400c-cc1f-4798-8045-e6fb93ec5647
transform!(
    mtg,
    [:Ra_PAR_f, :Ra_NIR_f] => ((x, y) -> x + y * 1.2) => :Rᵢ, # This would be the incident radiation
    [:Ra_PAR_f, :Ra_NIR_f] => ((x, y) -> x + y) => :Rₛ,
    :Ra_PAR_f => (x -> x * 4.57) => :PPFD,
    (x -> 0.3) => :d,
    ignore_nothing = true
)

# ╔═╡ 0ce68c8e-13a2-499d-86c1-59ea31c69f42
md"""
## Running the simulation
"""

# ╔═╡ 37f50a97-2842-4b12-ab70-ab56de09af0b
run!(mtg, models, weather)

# ╔═╡ 5376ccaf-a41a-4ad2-85a8-8b4badd4d938
md"""
## Quickly benchmarking

For sake of simplicity, we benchmark here in Sequential mode (i.e. no parallelization).
"""

# ╔═╡ b22d0e7c-aeaf-4c28-9cf9-5368a4384cf3
begin
	B = @benchmark run!($mtg, $models, $weather,executor=SequentialEx())
	B
end

# ╔═╡ 840c1748-0502-4d33-ad81-bf8047b58037
md"""
## Plotting the result in 3D
"""

# ╔═╡ 2fcdc6ba-0b5a-450a-ac4d-b013c9dfa2fd
begin
	transform!(
    mtg,
    :Tₗ => (x -> x[1]) => :Tₗ_1,
    ignore_nothing = true
)

	f, ax, p = viz(mtg, color = :Tₗ_1)
	colorbar(f[1, 2], p)
	f
end

# ╔═╡ Cell order:
# ╟─65ec1544-c81f-11ed-0a29-174eeb8e92e0
# ╟─56757555-725c-49b8-8ad9-84d99052a0c8
# ╟─8c41fd8a-648c-4d50-b39b-d6087b264d95
# ╠═552692ae-2f8e-4ae1-88a0-5565edabdf3f
# ╟─22f4ac9d-0181-476b-a3dd-277b2c5e31a2
# ╠═2a106a14-3ab9-4ebc-87bd-8c86453a3166
# ╟─c13c9dd7-3a8c-436a-9520-c638a76b7135
# ╠═00877758-9349-4aa6-a2c8-ac2e50e150a0
# ╟─e9293038-6179-481e-93c9-a619bcbce8ca
# ╟─f02f400c-cc1f-4798-8045-e6fb93ec5647
# ╟─0ce68c8e-13a2-499d-86c1-59ea31c69f42
# ╠═37f50a97-2842-4b12-ab70-ab56de09af0b
# ╟─5376ccaf-a41a-4ad2-85a8-8b4badd4d938
# ╠═b22d0e7c-aeaf-4c28-9cf9-5368a4384cf3
# ╟─840c1748-0502-4d33-ad81-bf8047b58037
# ╠═2fcdc6ba-0b5a-450a-ac4d-b013c9dfa2fd
