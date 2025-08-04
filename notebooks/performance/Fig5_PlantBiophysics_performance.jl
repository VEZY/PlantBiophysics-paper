### A Pluto.jl notebook ###
# v0.19.23

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
            PPFD=set.PPFD[i],
            d=set.d[i],
        ),
    )
    deps = PlantSimEngine.dep(leaf)
    meteo = Atmosphere(T=set.T[i], Wind=set.Wind[i], P=set.P[i], Rh=set.Rh[i], Cₐ=set.Ca[i])
    st = PlantMeteo.row_struct(leaf.status[1])
    b_PB = @benchmark run!($leaf, $deps, $st, $meteo, $constants, nothing) evals =
        microbenchmark_evals samples = microbenchmark_steps
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
	This is the plot from the latest commit on <https://github.com/VEZY/PlantBiophysics-paper/>. If you want to make your own benchmarking, run the script that was used to perform it, but careful, it takes a long time to perform!
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
        resolution=size_pt,
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

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "9573fa42e6a51a723efc4b333934638735eee6ec"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "81f0cb60dc994ca17f68d9fb7c942a5ae70d9ee4"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "5.0.8"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BSON]]
git-tree-sha1 = "306bb5574b0c1c56d7e1207581516c557d105cad"
uuid = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
version = "0.3.5"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CovarianceEstimation]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "a3e070133acab996660d31dcf479ea42849e368f"
uuid = "587fd27a-f159-11e8-2dae-1979310e6154"
version = "0.2.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataDeps]]
deps = ["BinaryProvider", "HTTP", "Libdl", "Reexport", "SHA", "p7zip_jll"]
git-tree-sha1 = "4f0e41ff461d42cfc62ff0de4f1cd44c6e6b3771"
uuid = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
version = "0.7.7"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "28d605d9a0ac17118fe2c5e9ce0fbb76c3ceb120"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "8a6b49396a4058771c5c072239b2e0a76e2e898c"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.58"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "0fa3b52a04a4e210aeb1626def9c90df3ae65268"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.1.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "51c8f36c81badaa0e9ec405dcbabaf345ed18c84"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.11.1"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "89cc49bf5819f0a10a7a3c38885e7c7ee048de57"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.29"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Impute]]
deps = ["BSON", "CSV", "DataDeps", "Distances", "IterTools", "LinearAlgebra", "Missings", "NamedDims", "NearestNeighbors", "Random", "Statistics", "StatsBase", "TableOperations", "Tables"]
git-tree-sha1 = "d3fb6342d94030706ad31f05c23514962c29296c"
uuid = "f7bf1975-0170-51b9-8c5f-a992d46b9575"
version = "0.6.8"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "76c987446e8d555677f064aaac1145c4c17662f8"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.14"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "91aa1442e63a77f101aff01dec5a821a17f43922"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.12.1"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MetaGraphsNext]]
deps = ["Graphs", "JLD2"]
git-tree-sha1 = "f8e0351036278130f6bf966cb903e5fddb93778c"
uuid = "fa8bd995-216d-47f1-8a91-f3b68fbeb377"
version = "0.2.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MultiScaleTreeGraph]]
deps = ["AbstractTrees", "DataFrames", "DelimitedFiles", "Graphs", "MetaGraphsNext", "MutableNamedTuples", "OrderedCollections", "Printf", "SHA", "XLSX"]
git-tree-sha1 = "d0dc30c70217c03211d690a6e16f77552aa58fb0"
uuid = "dd4a991b-8a45-4075-bede-262ee62d5583"
version = "0.5.0"

[[deps.MutableNamedTuples]]
git-tree-sha1 = "bd551bdbe872d8f132bb0c5fb6e6dec236fa01f9"
uuid = "af6c499f-54b4-48cc-bbd2-094bba7533c7"
version = "0.1.1"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NamedDims]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "0856b62716585eb90cc1dada226ac9eab5f69aa5"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "0.2.47"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ded92de95031d4a8c61dfb6ba9adb6f1d8016ddd"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "3114946c67ef9925204cc024a73c9e679cebe0d7"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.8"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlantBiophysics]]
deps = ["CSV", "DataFrames", "Dates", "Impute", "LsqFit", "MultiScaleTreeGraph", "MutableNamedTuples", "OrderedCollections", "RecipesBase", "Statistics", "Test", "YAML"]
git-tree-sha1 = "c8891016722498da99581b983f753c4d44f86ab6"
uuid = "7ae8fcfa-76ad-4ec6-9ea7-5f8f5e2d6ec9"
version = "0.3.0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "5309da1cdef03e95b73cd3251ac3a39f887da53e"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "cd56bf18ed715e8b09f06ef8c6b781e6cdc49911"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c82aaa13b44ea00134f8c9c89819477bd3986ecd"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.3.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "50ccd5ddb00d19392577902f0079267a72c5ab04"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.XLSX]]
deps = ["Dates", "EzXML", "Printf", "Tables", "ZipFile"]
git-tree-sha1 = "7fa8618da5c27fdab2ceebdff1da8918c8cd8b5d"
uuid = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"
version = "0.7.10"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "3c6e8b9f5cdaaa21340f841653942e1a6b6561e5"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.7"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "3593e69e469d2111389a9bd06bac1f3d730ac6de"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.4"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
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
