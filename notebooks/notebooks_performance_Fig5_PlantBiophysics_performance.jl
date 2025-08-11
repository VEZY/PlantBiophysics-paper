### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 211f1a10-061d-452a-b061-aa5eb5f07428
begin
    using Statistics
    using DataFrames
    using CSV
    using Random
    using PlantBiophysics
    using Downloads
end

# ╔═╡ ca26f490-bc68-4489-97cc-95a4c39097e5
HTML("   <div style=\"\n   position: absolute;\n   width: calc(100% - 30px);\n   border: 50vw solid #282936;\n   border-top: 500px solid #282936;\n   border-bottom: none;\n   box-sizing: content-box;\n   left: calc(-50vw + 15px);\n   top: -500px;\n   height: 180px;\n   pointer-events: none;\n   \"></div>\n\n   <div style=\"\n   height: 180px;\n   width: 100%;\n   background: #282936;\n   color: #fff;\n   padding-top: 68px;\n   \">\n   <span style=\"\n   font-family: Vollkorn, serif;\n   font-weight: 700;\n   font-feature-settings: 'lnum', 'pnum';\n   \"> \n\n   <p style=\"text-align: center; font-size: 2rem;\">\n   <em>Benchmark</em>\n   </p>\n\n\n\n   <style>\n   body {\n   overflow-x: hidden;\n   }\n   </style>")

# ╔═╡ 994a4062-39cf-49b7-b19e-63ff90e9a574
md"""

# PlantBiophysics.jl benchmark

The main objective of this notebook is to compare the computational times of `PlantBiophysics.jl` against the [plantecophys](https://github.com/RemkoDuursma/plantecophys) R package and the [LeafGasExchange.jl](https://github.com/cropbox/LeafGasExchange.jl) Julia package from the [Cropbox.jl](https://github.com/cropbox/Cropbox.jl) framework. The comparison follows three steps:
- create an N-large basis of random conditions.
- benchmark the computational time of the three packages via similar functions (_i.e._ photosynthesis-stomatal conductance-energy balance coupled model for C3 leaves): `energy_balance`, `photosynEB` and `simulate` with `ModelC3MD`.
- compare the results with plots and statistics.

This notebook does not perform the benchmark by itself for the obvious reason that it takes forever to run, and because there is an overhead cause by Pluto (benchmark and reactive are not a good mix). Instead, it shows the outputs of [a script from this repository](https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/Fig5_PlantBiophysics_performance_noPluto.jl) that implements the code shown here. If you want to perform the benchmark by yourself, you can run this script from the command line. The versions used for the above dependencies are available in the [Project.toml](https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/Project.toml) of the repository.


## Importing the dependencies:

###### Note
Make sure to have R installed on your computer first.

Loading the Julia packages:

```julia
begin
    using CairoMakie
	using BenchmarkTools
	using PlantBiophysics
    using Cropbox
    using LeafGasExchange
    using RCall
end
```
	
"""

# ╔═╡ 9c3f432b-4d1f-4765-9989-17e983cee5f7
md"""

## Parameters

### Benchmark parameters

You'll find below the main parameters of the benchmark. In a few words, each package runs a simulation for `N` different time-steps `microbenchmark_steps` times repeated `microbenchmark_evals` times. We make `N` different simulations because the simulation duration can vary depending on the inputs due to iterative computations in the code, *i.e.* different initial conditions can make the algorithms converge more or less rapidly.
"""

# ╔═╡ e7e23262-1813-4fe5-80fb-2be20cb54f89
begin
    Random.seed!(1) # Set random seed
    microbenchmark_steps = 100 # Number of times the microbenchmark is run
    microbenchmark_evals = 1 # N. times each sample is run to be sure of the output
    N = 100 # Number of timesteps simulated for each microbenchmark step
end

# ╔═╡ 65d00bf6-8f9a-40b4-947a-fd46ddc89753
md"""

### Random input simulation dataset

We create possible ranges for input parameters. These ranges where chosen so all of the three packages don't return errors during computation (`plantecophys` has issues with low temperatures).

- Ta: air temperature ($°C$)
- Wind: wind speed ($m.s^{-1}$)
- P: ambient pressure ($kPa$)
- Rh: relative humidity (between 0 and 1)
- Ca: air CO₂ concentration ($ppm$)
- Jmax: potential rate of electron transport ($\mu mol_{CO2}.m^{-2}.s^{-1}$)
- Vmax: maximum rate of Rubisco activity ($\mu mol_{CO2}.m^{-2}.s^{-1}$)
- Rd: mitochondrial respiration in the light at reference temperature ($\mu mol_{CO2}.m^{-2}.s^{-1}$)
- TPU: triose phosphate utilization-limited photosynthesis rate ($\mu mol_{CO2}.m^{-2}.s^{-1}$)
- Rs: short-wave net radiation ($W.m^{-1}$)
- skyF: Sun-visible fraction of the leaf (between 0 and 1)
- d: characteristic length ($m$)
- g0: residual stomatal conductance ($mol_{CO2}.m^{-2}.s^{-1}$)
- g1: slope of the stomatal conductance relationship.
"""

# ╔═╡ acc685f9-7ae8-483d-982e-ff45ccd9e860
begin
    # Create the ranges of input parameters
    length_range = 10000
    Rs = range(10, 500, length=length_range)
    Ta = range(18, 40, length=length_range)
    Wind = range(0.5, 20, length=length_range)
    P = range(90, 101, length=length_range)
    Rh = range(0.1, 0.98, length=length_range)
    Ca = range(360, 900, length=length_range)
    skyF = range(0.0, 1.0, length=length_range)
    d = range(0.001, 0.5, length=length_range)
    Jmax = range(200.0, 300.0, length=length_range)
    Vmax = range(150.0, 250.0, length=length_range)
    Rd = range(0.3, 2.0, length=length_range)
    TPU = range(5.0, 20.0, length=length_range)
    g0 = range(0.001, 2.0, length=length_range)
    g1 = range(0.5, 15.0, length=length_range)
    vars = hcat([Ta, Wind, P, Rh, Ca, Jmax, Vmax, Rd, Rs, skyF, d, TPU, g0, g1])
    nothing
end

# ╔═╡ b8126c39-8c6c-490d-9992-8922e88f8857
md"We then sample `N` conditions from the given ranges:"

# ╔═╡ a0476d0d-64d9-4457-b78d-41519e38e859
begin
    set = [rand.(vars) for i = 1:N]
    set = reshape(vcat(set...), (length(set[1]), length(set)))'
    name = [
        "T",
        "Wind",
        "P",
        "Rh",
        "Ca",
        "JMaxRef",
        "VcMaxRef",
        "RdRef",
        "Rs",
        "sky_fraction",
        "d",
        "TPURef",
        "g0",
        "g1",
    ]
    set = DataFrame(set, name)
    @. set[!, :vpd] = e_sat(set.T) - vapor_pressure(set.T, set.Rh)
    @. set[!, :PPFD] = set.Rs * 0.48 * 4.57
    set
end

# ╔═╡ ae00e44e-d864-4d04-b2e8-fac510ce44bb
md"
## Benchmarking
"

# ╔═╡ 8419c830-8e35-493b-91b2-adfc97e30a61
md"""
##### plantecophys

Preparing R to make the benchmark:

```julia
R\"\"\"
if(!require("plantecophys")){
    install.packages("plantecophys", repos = "https://cloud.r-project.org")
}
if(!require("microbenchmark")){
    install.packages("microbenchmark", repos = "https://cloud.r-project.org")
}
\"\"\"

# Make variables available to the R session
@rput set N microbenchmark_steps

```

Making the benchmark:

```julia
R\"\"\"
# Define the function call in a function that takes a list as input to limit DataFrame overhead
function_EB <- function(input) {
    PhotosynEB(
        Tair = input$Tair, VPD = input$VPD, Wind = input$Wind,
        Wleaf = input$Wleaf,Ca = input$Ca,  StomatalRatio = 1,
        LeafAbs = input$LeafAbs, gsmodel = "BBOpti", g0 = input$g0, g1 = input$g1,
        alpha = 0.24, theta = 0.7, Jmax = input$Jmax,
        Vcmax = input$Vcmax, TPU = input$TPU, Rd = input$Rd,
        RH = input$RH, PPFD=input$PPFD, Patm = input$Patm
    )
}

time_PE = c()
for(i in seq_len(N)){
    # Put the inputs into a vector to limit dataframe overhead:
    input = list(
        Tair = set$T[i], VPD = set$vpd[i], Wind = set$Wind[i], Wleaf = set$d[i],
        Ca = set$Ca[i], LeafAbs = set$sky_fraction[i], g0 = set$g0[i], g1 = set$g1[i],
        Jmax = set$JMaxRef[i], Vcmax = set$VcMaxRef[i], TPU = set$TPURef[i],
        Rd = set$RdRef[i], RH = set$Rh[i]*100, PPFD=set$PPFD[i],Patm = set$P[i]
    )

    m = microbenchmark(function_EB(input), times = microbenchmark_steps)

    time_PE = append(time_PE,m$time * 10e-9) # transform in seconds
}
\"\"\"

@rget time_PE
```

"""

# ╔═╡ 4dad9b75-3a03-4da1-9aa8-fea02fb48a97
md"""

##### LeafGasExchange.jl

Note that we benchmark `LeafGasExchange.jl` with the `nounit` flag to compute a fair comparison with `PlantBiophysics.jl` in case computing units takes time (it shouldn't much).

```julia
time_LG = []
n_lg = fill(0, N)
for i = 1:N
    config =
        :Weather => (
            PFD=set.PPFD[i],
            CO2=set.Ca[i],
            RH=set.Rh[i] * 100,
            T_air=set.T[i],
            wind=set.Wind[i],
            P_air=set.P[i],
            g0=set.g0[i],
            g1=set.g1[i],
            Vcmax=set.VcMaxRef[i],
            Jmax=set.JMaxRef[i],
            Rd=set.RdRef[i],
            TPU=set.TPURef[i],
        )
    b_LG =
        @benchmark simulate($ModelC3MD; config=$config) evals = microbenchmark_evals samples =
            microbenchmark_steps
    append!(time_LG, b_LG.times .* 1e-9) # transform in seconds
    n_lg[i] = 1
end
```

"""

# ╔═╡ 60f769b9-bf8c-461c-9739-f31ca43b27b8
md"""
##### PlantBiophysics.jl

Benchmarking `PlantBiophysics.jl`:

```julia
constants = Constants()
time_PB = []
for i = 1:N
    leaf = ModelList(
        energy_balance=Monteith(),
        photosynthesis=Fvcb(
            VcMaxRef=set.VcMaxRef[i],
            JMaxRef=set.JMaxRef[i],
            RdRef=set.RdRef[i],
            TPURef=set.TPURef[i],
        ),
        stomatal_conductance=Medlyn(set.g0[i], set.g1[i]),
        status=(
            Rₛ=set.Rs[i],
            sky_fraction=set.sky_fraction[i],
            aPPFD=set.PPFD[i],
            d=set.d[i],
        ),
    )
    meteo = Atmosphere(T=set.T[i], Wind=set.Wind[i], P=set.P[i], Rh=set.Rh[i], Cₐ=set.Ca[i])
    b_PB = @benchmark run!($leaf, $meteo, $constants, nothing) evals = microbenchmark_evals samples = microbenchmark_steps
    append!(time_PB, b_PB.times .* 1e-9) # transform in seconds
end
```
"""

# ╔═╡ 87cc52dd-9878-422d-a10e-d491db2a04c6
md"
## Comparison
"

# ╔═╡ 8294a065-5285-4221-8553-e38ca23257f2
md"""
##### Statistics

We compute here basic statistics, _i.e._ mean, median, min, max, standard deviation. 

```julia
statsPB = basic_stat(time_PB)
statsPE = basic_stat(time_PE)
statsLG = basic_stat(time_LG)

factorPE = mean(time_PE) / mean(time_PB)
factorLG = mean(time_LG) / mean(time_PB)

# Write overall timings:
df = DataFrame(
	[getfield(j, i) for i in fieldnames(StatResults), j in [statsPB, statsPE, statsLG]],
	["PlantBiophysics", "plantecophys", "LeafGasExchange"]
)
insertcols!(df, 1, :Stat => [fieldnames(StatResults)...])
CSV.write("benchmark.csv", df)

# Write timing for each sample:
CSV.write("benchmark_full.csv",
	DataFrame(
		"package" => vcat(
			[
				repeat([i.first], length(i.second)) for i in [
					"PlantBiophysics" => time_PB,
					"plantecophys" => time_PE,
					"LeafGasExchange" => time_LG
				]
			]...
		),
		"sample_time" => vcat(time_PB, time_PE, time_LG)
	)
)
```
"""

# ╔═╡ 59ce5f30-3b74-402d-919d-b1b3d1c76788
df_res = CSV.read(Downloads.download("https://raw.githubusercontent.com/VEZY/PlantBiophysics-paper/main/notebooks/performance/benchmark.csv"), DataFrame)

# ╔═╡ 12497dfb-87ff-4978-801a-c78cd118b032
md"
##### Histogram plotting
"

# ╔═╡ 84a4b3b4-23d8-4b12-9204-110d3ba4724b
Markdown.parse("""
We here display the computational time histogram of each package on the same scale in order to compare them: `PlantBiophysics.jl` (a), `plantecophys` (b) and `LeafGasExchange.jl` (c). The y-axis represents the density (_i.e._ reaching 0.3 means that 30% of the computed times are in this bar). Orange zone represents the interval [mean - standard deviation; mean + standard deviation]. Red dashed line represents the mean. Note that the x-axis is logarithmic.

```julia
    fig = plot_benchmark_Makie(statsPB, statsPE, statsLG, time_PB, time_PE, time_LG)
    save("benchmark_each_time_steps.png", fig, px_per_unit=3)
```

![](https://raw.githubusercontent.com/VEZY/PlantBiophysics-paper/main/notebooks/performance/benchmark_each_time_steps.png)

!!! note
	PlantBiophysics.jl is $(Int(round(df_res[1,:plantecophys] / df_res[1,:PlantBiophysics]))) times faster than plantecophys, and  $(Int(round(df_res[1,:LeafGasExchange] / df_res[1,:PlantBiophysics]))) times faster than LeafGasExchange.jl.

!!! warning
	This is the plot from the latest commit on <https://github.com/VEZY/PlantBiophysics-paper/>. If you want to make your own benchmark, run the [script used to perform it](https://github.com/VEZY/PlantBiophysics-paper/blob/main/notebooks/performance/Fig5_PlantBiophysics_performance_noPluto.jl), but be careful, it takes several minutes to run!
""")

# ╔═╡ b2958ece-3a4c-497e-88e3-afc45d8f3688
md"""
# References
"""

# ╔═╡ 22ec1978-aaee-4fa1-b739-1246a6aeda20
"""
	StatResults(
		mean::AbstractFloat
    	median::AbstractFloat
    	stddev::AbstractFloat
    	min::AbstractFloat
    	max::AbstractFloat
	)

Structure to hold basic statistics of model performance.
"""
struct StatResults
    mean::AbstractFloat
    median::AbstractFloat
    stddev::AbstractFloat
    min::AbstractFloat
    max::AbstractFloat
end

# ╔═╡ e56f829e-4b83-457a-8e90-a08e071cd5f3
begin
    function Base.show(io::IO, ::MIME"text/plain", m::StatResults)
        print(
            io,
            "Benchmark:",
            "\nMean time -> ",
            m.mean,
            " ± ",
            m.stddev,
            "\nMedian time -> ",
            m.median,
            "\nMinimum time -> ",
            m.min,
            "\nMaximum time -> ",
            m.max,
        )
    end
    Base.show(io::IO, m::StatResults) = print(io, m.mean, "(±", m.stddev, ')')


    md"""
    **Base.show**

    	Base.show(io::IO, m::StatResults)
    	Base.show(io::IO, ::MIME"text/plain", m::StatResults)

    Add a show method for our `StatResults` type.
    """
end

# ╔═╡ 30a52bf1-ff77-4f9b-a87a-e9dd8ad0c1c4
"""
	basic_stat(df)

Compute basic statistics from the benchmarking
"""
function basic_stat(df)
    m = mean(df)
    med = median(df)
    std = Statistics.std(df)
    min = findmin(df)[1]
    max = findmax(df)[1]
    return StatResults(m, med, std, min, max)
end

# ╔═╡ ac330756-f1ac-431c-9e60-7d346e06aa1f
function plot_benchmark_Makie(statsPB, statsPE, statsLG, time_PB, time_PE, time_LG)
    size_inches = (6.7, 5)
    size_pt = 72 .* size_inches
    bins = 220
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    fig = Figure(
        backgroundcolor=RGBf(1, 1, 1),
        size=size_pt,
        font=noto_sans,
        fontsize=10,
    )
    ep = 1e-9
    extr = extrema(vcat(time_PB, time_PE, time_LG))
    interval = (extr[1] * 1e-1, extr[2])

    ga = fig[1, 1] = GridLayout()

    axa = Axis(
        ga[1, 1],
        title="(a) PlantBiophysics.jl",
        xscale=log10,
        titlealign=:left,
        titlesize=10,
    )
    stddevi = poly!(
        axa,
        Rect(max(ep, statsPB.mean - statsPB.stddev), 0.0, 2 * statsPB.stddev, 1),
        color=(:orange, 0.3),
        yautolimits=false,
    )
    moy = vlines!(axa, statsPB.mean; color=:red, linewidth=3, linestyle=:dot)
    hist!(axa, time_PB, normalization=:probability, bins=bins)
    # h = axa.finallimits[].widths[2]
    axislegend(
        axa,
        [stddevi, moy],
        ["95% confidence interval", "Mean"],
        "",
        position=:rb,
        orientation=:vertical,
        labelsize=8,
        framevisible=false,
    )
    xlims!(axa, interval)

    axb = Axis(
        ga[2, 1],
        title="(b) plantecophys",
        xscale=log10,
        ylabel="Density",
        titlealign=:left,
        titlesize=10,
    )
    stddevi = poly!(
        axb,
        Rect(statsPE.mean - statsPE.stddev, 0.0, 2 * statsPE.stddev, 1),
        color=(:orange, 0.3),
        yautolimits=false,
    )
    vlines!(axb, statsPE.mean; color=:red, linewidth=3, linestyle=:dot)
    hist!(axb, time_PE, normalization=:probability, bins=bins)
    xlims!(axb, interval)

    # axc = Axis(ga[3, 1], title="(c) LeafGasExchange.jl", yminorticks=IntervalsBetween(10),
    #     xscale=log10, xminorticks=IntervalsBetween(10), yminorgridvisible=true, yminorticksvisible=true,
    #     xminorgridvisible=true, xminorticksvisible=true, xlabel="Time (s)")
    axc = Axis(
        ga[3, 1],
        title="(c) LeafGasExchange.jl",
        xscale=log10,
        xlabel="Time (s)",
        titlealign=:left,
        titlesize=10,
    )
    stddevi = poly!(
        axc,
        Rect(statsLG.mean - statsLG.stddev, 0.0, 2 * statsLG.stddev, 1),
        color=(:orange, 0.3),
        yautolimits=false,
    )
    vlines!(axc, statsLG.mean; color=:red, linewidth=3, linestyle=:dot)
    hist!(axc, time_LG, normalization=:probability, bins=bins)
    xlims!(axc, interval)

    rowgap!(ga, 7)
    hidexdecorations!(axa, grid=false)
    hidexdecorations!(axb, grid=false)
    fig
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
PlantBiophysics = "7ae8fcfa-76ad-4ec6-9ea7-5f8f5e2d6ec9"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.3.4"
PlantBiophysics = "~0.3.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "798794481d8bc7d1ef5cd7a2590726f9c0888ed2"

[[deps.ADTypes]]
git-tree-sha1 = "7927b9af540ee964cc5d1b73293f1eb0b761a3a1"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.16.0"

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

    [deps.ADTypes.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "9606d7832795cbef89e06a550475be300364a8aa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.19.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.BSON]]
git-tree-sha1 = "4c3e506685c527ac6a54ccc0c8c76fd6f91b42fb"
uuid = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
version = "0.3.9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "0037835448781bb46feb39866934e243886d756a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataDeps]]
deps = ["HTTP", "Libdl", "Reexport", "SHA", "Scratch", "p7zip_jll"]
git-tree-sha1 = "8ae085b71c462c2cb1cfedcb10c3c877ec6cf03f"
uuid = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
version = "0.7.13"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "db2a9cb664fcea7836da4b414c3278d71dd602d2"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.6"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "54d7b8c74408048aea9c055ac8573b2b5c5ec11f"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.4"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

    [deps.Distances.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3e6d038b77f22791b8e3472b7c633acea1ecac06"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.120"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "680a39c9aadce7c721b68d979e66dc65d2021aa6"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.2.2"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "b66970a70db13f45b7e57fbda1736e1cf72174ea"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.17.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "f089ab1f834470c525562030c8cfde4025d5e915"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.27.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Formatting]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "fb409abab2caf118986fc597ba84b50cbaf00b87"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.3"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "c5abfa0ae0aaee162a3fbb053c13ecda39be545b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Impute]]
deps = ["BSON", "CSV", "DataDeps", "Distances", "IterTools", "LinearAlgebra", "Missings", "NamedDims", "NearestNeighbors", "Random", "Statistics", "StatsBase", "TableOperations", "Tables"]
git-tree-sha1 = "ec5ca28d28bf8c807f1f8f9db885037ae660c89d"
uuid = "f7bf1975-0170-51b9-8c5f-a992d46b9575"
version = "0.6.12"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "8594fac023c5ce1ef78260f24d1ad18b4327b420"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.4"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "PrecompileTools", "Requires", "TranscodingStreams"]
git-tree-sha1 = "89e1e5c3d43078d42eed2306cab2a11b13e5c6ae"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.54"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "91aa1442e63a77f101aff01dec5a821a17f43922"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.12.1"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.MetaGraphsNext]]
deps = ["Graphs", "JLD2"]
git-tree-sha1 = "f8e0351036278130f6bf966cb903e5fddb93778c"
uuid = "fa8bd995-216d-47f1-8a91-f3b68fbeb377"
version = "0.2.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultiScaleTreeGraph]]
deps = ["AbstractTrees", "DataFrames", "DelimitedFiles", "Graphs", "MetaGraphsNext", "MutableNamedTuples", "OrderedCollections", "Printf", "SHA", "XLSX"]
git-tree-sha1 = "d0dc30c70217c03211d690a6e16f77552aa58fb0"
uuid = "dd4a991b-8a45-4075-bede-262ee62d5583"
version = "0.5.0"

[[deps.MutableNamedTuples]]
git-tree-sha1 = "0faaabea6ebbfde9a5a01455f851009fb2603aac"
uuid = "af6c499f-54b4-48cc-bbd2-094bba7533c7"
version = "0.1.3"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NamedDims]]
deps = ["LinearAlgebra", "Statistics"]
git-tree-sha1 = "f9e4a49ecd1ea2eccfb749a506fa882c094152b4"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "1.2.3"

    [deps.NamedDims.extensions]
    AbstractFFTsExt = "AbstractFFTs"
    ChainRulesCoreExt = "ChainRulesCore"
    CovarianceEstimationExt = "CovarianceEstimation"
    TrackerExt = "Tracker"

    [deps.NamedDims.weakdeps]
    AbstractFFTs = "621f4979-c628-5d54-868e-fcf4e3e8185c"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    CovarianceEstimation = "587fd27a-f159-11e8-2dae-1979310e6154"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ca7e18198a166a1f3eb92a3650d53d94ed8ca8a1"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.22"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "87510f7292a2b21aeff97912b0898f9553cc5c2c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlantBiophysics]]
deps = ["CSV", "DataFrames", "Dates", "Impute", "LsqFit", "MultiScaleTreeGraph", "MutableNamedTuples", "OrderedCollections", "RecipesBase", "Statistics", "Test", "YAML"]
git-tree-sha1 = "c8891016722498da99581b983f753c4d44f86ab6"
uuid = "7ae8fcfa-76ad-4ec6-9ea7-5f8f5e2d6ec9"
version = "0.3.0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "cbea8a6bd7bed51b1619658dec70035e07b8502f"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.14"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "8e45cecc66f3b42633b8ce14d431e8e57a3e242e"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "b765e46ba27ecf6b44faf70df40c57aa3a547dcb"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.7"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XLSX]]
deps = ["Dates", "EzXML", "Printf", "Tables", "ZipFile"]
git-tree-sha1 = "7fa8618da5c27fdab2ceebdff1da8918c8cd8b5d"
uuid = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"
version = "0.7.10"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "c217bad8fccb3bbfef7d7902326eacfbd0d702ad"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.14.4+0"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "2f58ac39f64b41fb812340347525be3b590cce3b"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.14"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "3593e69e469d2111389a9bd06bac1f3d730ac6de"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.4"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─ca26f490-bc68-4489-97cc-95a4c39097e5
# ╟─994a4062-39cf-49b7-b19e-63ff90e9a574
# ╠═211f1a10-061d-452a-b061-aa5eb5f07428
# ╟─9c3f432b-4d1f-4765-9989-17e983cee5f7
# ╠═e7e23262-1813-4fe5-80fb-2be20cb54f89
# ╟─65d00bf6-8f9a-40b4-947a-fd46ddc89753
# ╠═acc685f9-7ae8-483d-982e-ff45ccd9e860
# ╟─b8126c39-8c6c-490d-9992-8922e88f8857
# ╟─a0476d0d-64d9-4457-b78d-41519e38e859
# ╟─ae00e44e-d864-4d04-b2e8-fac510ce44bb
# ╟─8419c830-8e35-493b-91b2-adfc97e30a61
# ╟─4dad9b75-3a03-4da1-9aa8-fea02fb48a97
# ╟─60f769b9-bf8c-461c-9739-f31ca43b27b8
# ╟─87cc52dd-9878-422d-a10e-d491db2a04c6
# ╟─8294a065-5285-4221-8553-e38ca23257f2
# ╟─59ce5f30-3b74-402d-919d-b1b3d1c76788
# ╟─12497dfb-87ff-4978-801a-c78cd118b032
# ╟─84a4b3b4-23d8-4b12-9204-110d3ba4724b
# ╟─b2958ece-3a4c-497e-88e3-afc45d8f3688
# ╟─22ec1978-aaee-4fa1-b739-1246a6aeda20
# ╟─e56f829e-4b83-457a-8e90-a08e071cd5f3
# ╟─30a52bf1-ff77-4f9b-a87a-e9dd8ad0c1c4
# ╟─ac330756-f1ac-431c-9e60-7d346e06aa1f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
