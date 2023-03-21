using PlantGeom
using PlantBiophysics
using BenchmarkTools
using PlantSimEngine, PlantMeteo
using Dates
using FLoops
using MultiScaleTreeGraph

# Reading MTG file
mtg = read_opf(joinpath(dirname(dirname(pathof(PlantBiophysics))), "test", "inputs", "scene", "opf", "coffee.opf"))

# Reading meteorological data
weather = PlantMeteo.read_weather(
    joinpath(dirname(dirname(pathof(PlantMeteo))), "test", "data", "meteo.csv"),
    :temperature => :T,
    :relativeHumidity => (x -> x ./ 100) => :Rh,
    :wind => :Wind,
    :atmosphereCO2_ppm => :Cₐ,
    date_format = DateFormat("yyyy/mm/dd")
)

# Reading models
file = joinpath(dirname(dirname(pathof(PlantBiophysics))), "test", "inputs", "models", "plant_coffee.yml")
models = read_model(file)

to_initialize(models)

# Add incident radiation
transform!(
    mtg,
    [:Ra_PAR_f, :Ra_NIR_f] => ((x, y) -> x + y * 1.2) => :Rᵢ, # This would be the incident radiation
    [:Ra_PAR_f, :Ra_NIR_f] => ((x, y) -> x + y) => :Rₛ,
    :Ra_PAR_f => (x -> x * 4.57) => :PPFD,
    (x -> 0.3) => :d,
    ignore_nothing = true
)

# Run the models
run!(mtg, weather)

# Benchmark it in sequential execution
B = @benchmark run!($mtg, $weather,executor=SequentialEx())

# Plot
transform!(
    mtg,
    :Tₗ => (x -> x[1]) => :Tₗ_1,
    ignore_nothing = true
)

f, ax, p = viz(mtg, color = :Tₗ_1)
colorbar(f[1, 2], p)
f