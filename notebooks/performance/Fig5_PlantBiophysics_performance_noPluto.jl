# To run this script, activate the environment in tutorials:
# ]activate .

using Random, Statistics, DataFrames, BenchmarkTools, CairoMakie, CSV
using Cropbox, LeafGasExchange, RCall, PlantBiophysics, PlantMeteo, PlantSimEngine

# Parameters :
Random.seed!(1)
microbenchmark_steps = 100 # Number of times the microbenchmark is run
microbenchmark_evals = 1 # Number of times each sample is run to be sure of the output
N = 100        # Number of timesteps simulated for each microbenchmark step

# Create the ranges of input parameters
function make_meteo(N; length_range=10000)
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

    # Create the randomly drawn input set
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

    return set
end

set = make_meteo(N)
######################################################################################################
# BENCHMARKING PLANTECOPHYS
######################################################################################################

# Import dependencies, install if not installed (but please prefer installing from R directly):
R"""
if(!require("plantecophys")){
    install.packages("plantecophys", repos = "https://cloud.r-project.org")
}
if(!require("microbenchmark")){
    install.packages("microbenchmark", repos = "https://cloud.r-project.org")
}
"""

# Make variables available to the R session
@rput set N microbenchmark_steps

# Benchmark:
R"""
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
"""

@rget time_PE


######################################################################################################
# BENCHMARKING LEAFGASEXCHANGE.JL
#############1#########################################################################################

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

######################################################################################################
# BENCHMARKING PLANTBIOPHYSICS.JL
######################################################################################################

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
    # run!(leaf, meteo, constants, nothing)
    b_PB = @benchmark run!($leaf, $meteo, $constants, nothing) evals = microbenchmark_evals samples = microbenchmark_steps
    append!(time_PB, b_PB.times .* 1e-9) # transform in seconds
end

######################################################################################################
# STATISTICS
######################################################################################################

# Defining a structure to hold basic statistics
struct StatResults
    mean::AbstractFloat
    median::AbstractFloat
    stddev::AbstractFloat
    min::AbstractFloat
    max::AbstractFloat
end

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

# Function computing the basic statistics
function basic_stat(df)
    m = mean(df)
    med = median(df)
    std = Statistics.std(df)
    min = findmin(df)[1]
    max = findmax(df)[1]
    return StatResults(m, med, std, min, max)
end

statsPB = basic_stat(time_PB)
statsPE = basic_stat(time_PE)
statsLG = basic_stat(time_LG)

factorPE = mean(time_PE) / mean(time_PB)
factorLG = mean(time_LG) / mean(time_PB)

######################################################################################################
# PLOTTING
######################################################################################################

function plot_benchmark_Makie(
    statsPB,
    statsPE,
    statsLG,
    time_PB,
    time_PE,
    time_LG;
    bins=220,
    backgroundcolor=:transparent
)

    size_inches = (6.7, 5)
    size_pt = 72 .* size_inches
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    fig = Figure(
        backgroundcolor=backgroundcolor,
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
        backgroundcolor=backgroundcolor,
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
        " ",
        position=:rb,
        orientation=:vertical,
        labelsize=8,
        framevisible=false,
    )
    CairoMakie.xlims!(axa, interval)

    axb = Axis(
        ga[2, 1],
        title="(b) plantecophys",
        xscale=log10,
        ylabel="Density",
        titlealign=:left,
        titlesize=10,
        backgroundcolor=backgroundcolor,
    )
    stddevi = poly!(
        axb,
        Rect(statsPE.mean - statsPE.stddev, 0.0, 2 * statsPE.stddev, 1),
        color=(:orange, 0.3),
        yautolimits=false,
    )
    vlines!(axb, statsPE.mean; color=:red, linewidth=3, linestyle=:dot)
    hist!(axb, time_PE, normalization=:probability, bins=bins)
    CairoMakie.xlims!(axb, interval)

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
        backgroundcolor=backgroundcolor,
    )
    stddevi = poly!(
        axc,
        Rect(statsLG.mean - statsLG.stddev, 0.0, 2 * statsLG.stddev, 1),
        color=(:orange, 0.3),
        yautolimits=false,
    )
    vlines!(axc, statsLG.mean; color=:red, linewidth=3, linestyle=:dot)
    hist!(axc, time_LG, normalization=:probability, bins=bins)
    CairoMakie.xlims!(axc, interval)

    rowgap!(ga, 7)
    hidexdecorations!(axa, grid=false)
    hidexdecorations!(axb, grid=false)
    fig
end

######################################################################################################
# SAVING
######################################################################################################

# Save the figure:
fig = plot_benchmark_Makie(
    statsPB,
    statsPE,
    statsLG,
    time_PB,
    time_PE,
    time_LG,
    bins=1000,
    backgroundcolor=:white,
)

save("benchmark_each_time_steps.png", fig, px_per_unit=6)

# Write overall timings:
df = DataFrame(
    [getfield(j, i) for i in fieldnames(StatResults), j in [statsPB, statsPE, statsLG]],
    ["PlantBiophysics", "plantecophys", "LeafGasExchange"],
)
insertcols!(df, 1, :Stat => [fieldnames(StatResults)...])
CSV.write("benchmark.csv", df)

# Write timing for each sample:
CSV.write(
    "benchmark_full.csv",
    DataFrame(
        "package" => vcat(
            [
                repeat([i.first], length(i.second)) for i in [
                    "PlantBiophysics" => time_PB,
                    "plantecophys" => time_PE,
                    "LeafGasExchange" => time_LG,
                ]
            ]...,
        ),
        "sample_time" => vcat(time_PB, time_PE, time_LG),
    ),
)


######################################################################################################
# Benchmark with limiting overhead
######################################################################################################

# Benchmark but passing all time-steps at once and computing the average simulation time.
# This is helpfull to avoid the overhead of calling the wrapper function, the downside is that
# we don't have the variability between simulations

# Parameters :
Random.seed!(1)
microbenchmark_steps2 = 100 # Number of times the microbenchmark is run
microbenchmark_evals2 = 1 # Number of times each sample is run to be sure of the output
N2 = 100        # Number of timesteps simulated for each microbenchmark step
set2 = make_meteo(N2)

# Make variables available to the R session
@rput set2 N2 microbenchmark_steps2
R"""
m2 = microbenchmark(
    PhotosynEB(
        Tair = set2$T, VPD = set2$vpd, Wind = set2$Wind,
        Wleaf = set2$d, Ca = set2$Ca,  StomatalRatio = 1,
        LeafAbs = set2$sky_fraction, gsmodel = "BBOpti", g0 = set2$g0, g1 = set2$g1,
        alpha = 0.24, theta = 0.7, Jmax = set2$JMaxRef,
        Vcmax = set2$VcMaxRef, TPU = set2$TPURef, Rd = set2$RdRef,
        RH = set2$Rh, PPFD=set2$PPFD, Patm = set2$P
    ),
    times = microbenchmark_steps2
)

    time_PE2 = m2$time / N2 * 10e-9 # transform in seconds
"""
@rget time_PE2


# Configuration:
leaves = [
    ModelList(
        Monteith(), Medlyn(set2.g0[i], set2.g1[i]),
        Fvcb(VcMaxRef=set2.VcMaxRef[i], JMaxRef=set2.JMaxRef[i], RdRef=set2.RdRef[i], TPURef=set2.TPURef[i]),
        status=(
            Rₛ=set2.Rs[i], sky_fraction=set2.sky_fraction[i],
            aPPFD=set2.PPFD[i], d=set2.d[i]
        )
    ) for i in 1:size(set2, 1)
]

meteos = Weather([Atmosphere(T=set2.T[i], Wind=set2.Wind[i], P=set2.P[i], Rh=set2.Rh[i], Cₐ=set2.Ca[i]) for i in 1:size(set2, 1)])
constants = Constants()
# b_PB2 = @benchmark eval_PB($leaves, $meteos) evals = microbenchmark_evals2 samples = microbenchmark_steps2
extra = nothing
b_PB2 = @benchmark run!($leaves, $meteos, $constants, $extra) evals = microbenchmark_evals2 samples = microbenchmark_steps2
time_PB2 = b_PB2.times ./ N2 .* 1e-9 # transform in seconds
mean(time_PB2) * 1e6 # 600 μs per timestep
#  Range (min … max):  34.426 ms … 149.017 ms  ┊ GC (min … max): 49.40% … 83.39%
#  Time  (median):     44.679 ms               ┊ GC (median):    55.17%
#  Time  (mean ± σ):   51.614 ms ±  21.050 ms  ┊ GC (mean ± σ):  60.09% ±  7.45%

#     ▁▃█▇▆  ▁                                                    
#   ▅▇██████▅█▅▅▇▅▁▁▇█▅▁▅▅▁▁▁▁▁▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▁▁▅▁▁▁▁▁▅▁▁▁▁▅ ▁
#   34.4 ms       Histogram: log(frequency) by time       144 ms <

#  Memory estimate: 646.03 MiB, allocs estimate: 740361.
# For LeafGasExchange:
# This is useful to know the overhead cost from Julia to Python and back.

configs = [
    begin
        :Weather => (
            PFD=set2.PPFD[i],
            CO2=set2.Ca[i],
            RH=set2.Rh[i] * 100,
            T_air=set2.T[i],
            wind=set2.Wind[i],
            P_air=set2.P[i],
            g0=set2.g0[i],
            g1=set2.g1[i],
            Vcmax=set2.VcMaxRef[i],
            Jmax=set2.JMaxRef[i],
            Rd=set2.RdRef[i],
            TPU=set2.TPURef[i]
        )
    end for i in 1:size(set2, 1)
]
b_LG2 = @benchmark simulate($ModelC3MD; configs=$configs) evals = microbenchmark_evals2 samples = microbenchmark_steps2
time_LG2 = b_LG2.times ./ N2 .* 1e-9 # transform in seconds


######################################################################################################
# SAVING
######################################################################################################

statsPB2 = basic_stat(time_PB2)
statsPE2 = basic_stat(time_PE2)
statsLG2 = basic_stat(time_LG2)

plot_benchmark(statsPB2, statsPE2, statsLG2)
savefig("benchmark_all_steps.png")





# Benchmarking PlantBiophysics on either one leaf and 1000 timesteps, or 1000 leaves and one timestep.

microbenchmark_steps3 = 100 # Number of times the microbenchmark is run
microbenchmark_evals3 = 1 # Number of times each sample is run to be sure of the output

# One time-step, 1000 leaves:
N_ts = 1
N_leaves = 1000
set3 = make_meteo(N_ts)
i = 1
leaves = [
    ModelList(
        Monteith(), Medlyn(set3.g0[i], set3.g1[i]),
        Fvcb(VcMaxRef=set3.VcMaxRef[i], JMaxRef=set3.JMaxRef[i], RdRef=set3.RdRef[i], TPURef=set3.TPURef[i]),
        status=(
            Rₛ=set3.Rs[i], sky_fraction=set3.sky_fraction[i],
            aPPFD=set3.PPFD[i], d=set3.d[i]
        )
    ) for l in 1:1000
]

meteos = Atmosphere(T=set3.T[i], Wind=set3.Wind[i], P=set3.P[i], Rh=set3.Rh[i], Cₐ=set3.Ca[i])
constants = Constants()
extra = nothing
b_PB3 = @benchmark run!($leaves, $meteos, $constants, $extra) evals = microbenchmark_evals3 samples = microbenchmark_steps3
time_PB_one_ts_1000_leaves_average = mean(b_PB3.times) .* 1e-6 # transform in ms for 1 time-step, 1000 leaves
time_PB_one_ts_1000_leaves = time_PB_one_ts_1000_leaves_average / N_ts / N_leaves # in ms per timestep per leaf
time_PB_one_ts_1000_leaves * 1e3 # 87.9 μs per timestep per leaf

# 1000 time-steps, 1 leaf:
N_ts = 1000
set4 = make_meteo(N_ts)
i = 1
leaves =
    ModelList(
        Monteith(), Medlyn(set4.g0[i], set4.g1[i]),
        Fvcb(VcMaxRef=set4.VcMaxRef[i], JMaxRef=set4.JMaxRef[i], RdRef=set4.RdRef[i], TPURef=set4.TPURef[i]),
        status=(
            Rₛ=set4.Rs[i], sky_fraction=set4.sky_fraction[i],
            aPPFD=set4.PPFD[i], d=set4.d[i]
        )
    )
meteos = Weather([Atmosphere(T=set4.T[i], Wind=set4.Wind[i], P=set4.P[i], Rh=set4.Rh[i], Cₐ=set4.Ca[i]) for i in 1:N_ts])
constants = Constants()
# b_PB2 = @benchmark eval_PB($leaves, $meteos) evals = microbenchmark_evals2 samples = microbenchmark_steps2
extra = nothing
b_PB4 = @benchmark run!($leaves, $meteos, $constants, $extra) evals = microbenchmark_evals3 samples = microbenchmark_steps3
time_PB_1000_ts_one_leaf_mean = mean(b_PB4.times) .* 1e-6 # in ms
time_PB_1000_ts_one_leaf = time_PB_1000_ts_one_leaf_mean / N_ts # transform in seconds
mean(time_PB_1000_ts_one_leaf) * 1e3 # 26 μs per timestep

# 1000 time-steps, 1000 leaves:
N_ts = 1000
set5 = make_meteo(N_ts)
i = 1
leaves = [
    ModelList(
        Monteith(), Medlyn(set5.g0[i], set5.g1[i]),
        Fvcb(VcMaxRef=set5.VcMaxRef[i], JMaxRef=set5.JMaxRef[i], RdRef=set5.RdRef[i], TPURef=set5.TPURef[i]),
        status=(
            Rₛ=set5.Rs[i], sky_fraction=set5.sky_fraction[i],
            aPPFD=set5.PPFD[i], d=set5.d[i]
        )
    ) for l in 1:1000
]

meteos = Weather([Atmosphere(T=set4.T[i], Wind=set4.Wind[i], P=set4.P[i], Rh=set4.Rh[i], Cₐ=set4.Ca[i]) for i in 1:N_ts])
constants = Constants()
extra = nothing

b_PB5 = @benchmark run!($leaves, $meteos, $constants, $extra) evals = microbenchmark_evals3 samples = microbenchmark_steps3
time_PB_1000_ts_1000_leaves_average = mean(b_PB5.times) .* 1e-9 # in s
time_PB_1000_ts_1000_leaves = time_PB_1000_ts_1000_leaves_average / N_ts / N_leaves * 1e3 # in ms per timestep per leaf

# 1 time-step, 1 leaf:
N_ts = 1
set6 = make_meteo(N_ts)
i = 1
leaves =
    ModelList(
        Monteith(), Medlyn(set6.g0[i], set6.g1[i]),
        Fvcb(VcMaxRef=set6.VcMaxRef[i], JMaxRef=set6.JMaxRef[i], RdRef=set6.RdRef[i], TPURef=set6.TPURef[i]),
        status=(
            Rₛ=set6.Rs[i], sky_fraction=set6.sky_fraction[i],
            aPPFD=set6.PPFD[i], d=set6.d[i]
        )
    )
meteos = Atmosphere(T=set6.T[i], Wind=set6.Wind[i], P=set6.P[i], Rh=set6.Rh[i], Cₐ=set6.Ca[i])
constants = Constants()
# b_PB2 = @benchmark eval_PB($leaves, $meteos) evals = microbenchmark_evals2 samples = microbenchmark_steps2
extra = nothing
b_PB6 = @benchmark run!($leaves, $meteos, $constants, $extra) evals = microbenchmark_evals3 samples = microbenchmark_steps3
time_PB_one_ts_one_leaf_mean = mean(b_PB6.times) .* 1e-6 # in ms
time_PB_one_ts_one_leaf = time_PB_one_ts_one_leaf_mean / N_ts
mean(time_PB_one_ts_one_leaf) * 1e3 # 26 μs per timestep
