using Random, Statistics, LaTeXStrings, DataFrames, BenchmarkTools, Plots, CSV
using Cropbox, LeafGasExchange
using RCall
using PlantBiophysics

# Parameters :
Random.seed!(1)
microbenchmark_steps = 100 # Number of times the microbenchmark is run
microbenchmark_evals = 1 # Number of times each sample is run to be sure of the output
N = 100        # Number of timesteps simulated for each microbenchmark step

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

# Create the randomly drawn input set
set = [rand.(vars) for i in 1:N]
set = reshape(vcat(set...), (length(set[1]), length(set)))'
name = ["T", "Wind", "P", "Rh", "Ca", "JMaxRef", "VcMaxRef", "RdRef", "Rs", "sky_fraction", "d", "TPURef", "g0", "g1"]
set = DataFrame(set, name)
@. set[!, :vpd] = e_sat(set.T) - vapor_pressure(set.T, set.Rh)
@. set[!, :PPFD] = set.Rs * 0.48 * 4.57

######################################################################################################
# BENCHMARKING PLANTECOPHYS
######################################################################################################

# Import dependencies, install if not installed:
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
######################################################################################################

time_LG = []
n_lg = fill(0, N)
for i in 1:N
    config = :Weather => (
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
        TPU=set.TPURef[i]
    )
    b_LG = @benchmark simulate($ModelC3MD; config=$config) evals = microbenchmark_evals samples = microbenchmark_steps
    append!(time_LG, b_LG.times .* 1e-9) # transform in seconds
    n_lg[i] = 1
end

######################################################################################################
# BENCHMARKING PLANTBIOPHYSICS.JL
######################################################################################################

time_PB = []
for i in 1:N
    leaf = LeafModels(
        energy=Monteith(),
        photosynthesis=Fvcb(
            VcMaxRef=set.VcMaxRef[i],
            JMaxRef=set.JMaxRef[i],
            RdRef=set.RdRef[i],
            TPURef=set.TPURef[i]
        ),
        stomatal_conductance=Medlyn(set.g0[i], set.g1[i]),
        Rₛ=set.Rs[i], sky_fraction=set.sky_fraction[i], PPFD=set.PPFD[i], d=set.d[i]
    )
    meteo = Atmosphere(T=set.T[i], Wind=set.Wind[i], P=set.P[i], Rh=set.Rh[i], Cₐ=set.Ca[i])

    b_PB = @benchmark energy_balance!($leaf, $meteo) evals = microbenchmark_evals samples = microbenchmark_steps
    append!(time_PB, b_PB.times .* 1e-9) # transform in seconds
end


######################################################################################################
# STATISTICS
######################################################################################################

# Defining a structure to hold basic statistics
struct stat_res
    mean::AbstractFloat
    median::AbstractFloat
    stddev::AbstractFloat
    min::AbstractFloat
    max::AbstractFloat
end

# Function computing the basic statistics
function basic_stat(df)
    m = mean(df)
    med = median(df)
    std = Statistics.std(df)
    min = findmin(df)[1]
    max = findmax(df)[1]
    return stat_res(m, med, std, min, max)
end

statsPB = basic_stat(time_PB)
statsPE = basic_stat(time_PE)
statsLG = basic_stat(time_LG)

factorPE = mean(time_PE) / mean(time_PB)
factorLG = mean(time_LG) / mean(time_PB)

######################################################################################################
# PLOTTING
######################################################################################################

function plot_benchmark(statsPB, statsPE, statsLG)
    Plots.plot(layout=grid(3, 1), xminorgrid=true) # Using `Plots.plot` as `plot` function is defined in Cropbox
    interval = (1e-6, 0.0) # x-axis
    # Ribbon plots (i.e. plot the standard deviation)
    eps = 1e-9 # limits for PlantBiophysics.jl standard deviation: sometimes the lower limit of the interval is negative and there are problems with logscale
    Plots.plot!(
        [
            (max(eps, statsPB.mean - statsPB.stddev), 0.0),
            (statsPB.mean + statsPB.stddev, 0.0),
            (statsPB.mean + statsPB.stddev, 0.9),
            (max(eps, statsPB.mean - statsPB.stddev), 0.9),
            (max(eps, statsPB.mean - statsPB.stddev), 0.0)
        ],
        seriestype=:shape, fillcolor=:orange, alpha=0.3, linecolor=:blue, linewidth=0.0, sp=1, label="95% confidence interval"
    )
    Plots.plot!(
        [
            (statsPE.mean - statsPE.stddev, 0.0),
            (statsPE.mean + statsPE.stddev, 0.0),
            (statsPE.mean + statsPE.stddev, 1.0),
            (statsPE.mean - statsPE.stddev, 1.0),
            (statsPE.mean - statsPE.stddev, 0.0)
        ],
        seriestype=:shape, fillcolor=:orange, alpha=0.3, linecolor=:blue, linewidth=0.0, sp=2, label=""
    )
    Plots.plot!(
        [
            (statsLG.mean - statsLG.stddev, 0.0),
            (statsLG.mean + statsLG.stddev, 0.0),
            (statsLG.mean + statsLG.stddev, 1.0),
            (statsLG.mean - statsLG.stddev, 1.0),
            (statsLG.mean - statsLG.stddev, 0.0)
        ],
        seriestype=:shape, fillcolor=:orange,
        alpha=0.3, linecolor=:blue, linewidth=0.0, sp=3, label="")

    # Histogram densities plots
    histogram!(time_PB[time_PB.<1e-4], sp=1, xaxis=(:log10, interval), normalize=:probability, legend=false, bins=100, dpi=300, color=:lightblue, label="")
    histogram!(time_PE[time_PE.<1e-0], sp=2, xaxis=(:log10, interval), normalize=:probability, bins=30, dpi=300, color=:lightblue, label="")
    histogram!(time_LG[time_LG.<1e-1], sp=3, xaxis=(:log10, interval), normalize=:probability, bins=30, dpi=300, color=:lightblue, label="")

    # Mean plotting
    Plots.plot!(statsPB.mean * [1, 1], [0.0, 0.9], sp=1, linewidth=2, linestyle=:dot, color=:red, label="Mean", legend=:topright)
    Plots.plot!(statsPE.mean * [1, 1], [0.0, 1.0], sp=2, linewidth=2, linestyle=:dot, color=:red, label="")
    Plots.plot!(statsLG.mean * [1, 1], [0.0, 1.0], sp=3, linewidth=2, linestyle=:dot, color=:red, label="")

    # Annotating (a), (b), (c)
    annotate!(1.5e-6, 1.0, Plots.text("(a) PlantBiophysics.jl", :black, :left, 9), sp=1)
    annotate!(1.5e-6, 1.0, Plots.text("(b) plantecophys", :black, :left, 9), sp=2)
    annotate!(1.5e-6, 1.0, Plots.text("(c) LeafGasExchange.jl", :black, :left, 9), sp=3)
    Plots.xlabel!("time (s)", xguidefontsize=10, sp=3)
    Plots.ylabel!("density", xguidefontsize=10, sp=2)
    Plots.xticks!([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
    Plots.ylims!(0.0, 1.0)
end


function plot_benchmark_Makie(statsPB,statsPE,statsLG)
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    fig = Makie.Figure(backgroundcolor = RGBf(1, 1, 1),resolution = (1000, 800),font=noto_sans,dpi=600,size=(1000,1000),xminorgridstyle=true)
    ep=1e-9
    interval = (3e-6, 2e-1) # x-axis

    axa = Axis(fig[1,1:3],title="(a) PlantBiophysics.jl",yminorticks= IntervalsBetween(10),xscale=log10,xminorticks= IntervalsBetween(10),yminorgridvisible=true,yminorticksvisible = true,xminorgridvisible=true,xminorticksvisible = true)
    stddevi = Makie.poly!(axa,Rect(max(ep, statsPB.mean - statsPB.stddev), 0., 2*statsPB.stddev, 0.2),color=(:orange, 0.3))
    Makie.hist!(axa,time_PB[time_PB.<1e-4], normalization=:probability, bins=300)
    moy = Makie.lines!(axa,statsPB.mean * [1, 1], [0.0, 0.2], linewidth=3, linestyle=:dot, color=:red)
    axislegend(axa, [stddevi,moy], ["95% confidence interval","Mean"], "", position = :rb,
        orientation = :vertical,labelsize=13,framevisible = false)
    Makie.xlims!(axa,interval)

    axb = Axis(fig[2,1:3],title="(b) plantecophys",yminorticks= IntervalsBetween(10),xscale=log10,xminorticks= IntervalsBetween(10),yminorgridvisible=true,yminorticksvisible = true,xminorgridvisible=true,xminorticksvisible = true)
    stddevi = Makie.poly!(axb,Rect(statsPE.mean - statsPE.stddev, 0., 2*statsPE.stddev, 0.2),color=(:orange, 0.3))
    Makie.hist!(axb,time_PE[time_PE.<1e-0],normalization=:probability, bins=100)
    moy = Makie.lines!(axb,statsPE.mean * [1, 1], [0.0, 0.2], linewidth=3, linestyle=:dot, color=:red)
    #axislegend(axb, [stddevi,moy], ["95% confidence interval","Mean"], "", position = :rb,
    #    orientation = :vertical,labelsize=13,framevisible = false)
    Makie.xlims!(axb,interval)

    axc = Axis(fig[3,1:3],title="(c) LeafGasExchange.jl",yminorticks= IntervalsBetween(10),xscale=log10,xminorticks= IntervalsBetween(10),yminorgridvisible=true,yminorticksvisible = true,xminorgridvisible=true,xminorticksvisible = true)
    stddevi = Makie.poly!(axc,Rect(statsLG.mean - statsLG.stddev, 0., 2*statsLG.stddev, 0.2),color=(:orange, 0.3))
    Makie.hist!(axc,time_LG[time_LG.<1e-1],normalization=:probability, bins=80)
    moy = Makie.lines!(axc,statsLG.mean * [1, 1], [0.0, 0.2], linewidth=3, linestyle=:dot, color=:red)
    #axislegend(axb, [stddevi,moy], ["95% confidence interval","Mean"], "", position = :rb,
    #    orientation = :vertical,labelsize=13,framevisible = false)
    Makie.xlims!(axc,interval)

    Label(fig[1:3, 0], "Density", rotation = pi/2,labelsize=25)
    Label(fig[4,1:4], "Time (s)")

    fig
end

######################################################################################################
# SAVING
######################################################################################################

plot_benchmark(statsPB, statsPE, statsLG)
savefig("benchmark_each_time_steps.png")

fig = plot_benchmark_Makie(statsPB, statsPE, statsLG)
save("/Users/simon/Desktop/perf.pdf", fig, px_per_unit = 2) # output size = 1600 x 1200 pixels
# perfs = DataFrame("TimePB" => time_PB, "TimeLGE" => time_LG, "TimePE" => time_PE)
# perfs = hcat(set, perfs)
# CSV.write("Simulation_times_" * string(N) * ".csv", perfs)


######################################################################################################
# Benchmark with limiting overhead
######################################################################################################

# Benchmark but passing all time-steps at once and computing the average simulation time.
# This is helpfull to avoid the overhead of calling the wrapper function, the downside is that
# we don't have the variability between simulations

# Parameters :
Random.seed!(1)
microbenchmark_steps2 = 1000 # Number of times the microbenchmark is run
microbenchmark_evals2 = 1 # Number of times each sample is run to be sure of the output
N2 = 100        # Number of timesteps simulated for each microbenchmark step

# Create the randomly drawn input set
set2 = [rand.(vars) for i in 1:N2]
set2 = reshape(vcat(set2...), (length(set2[1]), length(set2)))'
name = ["T", "Wind", "P", "Rh", "Ca", "JMaxRef", "VcMaxRef", "RdRef", "Rs", "sky_fraction", "d", "TPURef", "g0", "g1"]
set2 = DataFrame(set2, name)
@. set2[!, :vpd] = e_sat(set2.T) - vapor_pressure(set2.T, set2.Rh)
@. set2[!, :PPFD] = set2.Rs * 0.48 * 4.57

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
    begin
        LeafModels(
            energy=Monteith(),
            photosynthesis=Fvcb(VcMaxRef=set2.VcMaxRef[i], JMaxRef=set2.JMaxRef[i], RdRef=set2.RdRef[i], TPURef=set2.TPURef[i]),
            stomatal_conductance=Medlyn(set2.g0[i], set2.g1[i]),
            Rₛ=set2.Rs[i], sky_fraction=set2.sky_fraction[i],
            PPFD=set2.PPFD[i], d=set2.d[i]
        )
    end for i in 1:size(set2, 1)
]

meteos = [
    begin
        Atmosphere(T=set2.T[i], Wind=set2.Wind[i], P=set2.P[i], Rh=set2.Rh[i], Cₐ=set2.Ca[i])
    end for i in 1:size(set2, 1)
]

function eval_PB(leaves, meteos)
    for i in 1:N2
        energy_balance!(leaves[i], meteos[i])
    end
end

b_PB2 = @benchmark eval_PB($leaves, $meteos) evals = microbenchmark_evals2 samples = microbenchmark_steps2
time_PB2 = b_PB2.times ./ N2 .* 1e-9 # transform in seconds


# For LeafGasExchange:
# This is usefull to know the overhead cost from Julia to Python and back.
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

# perfs2 = DataFrame("TimePB" => time_PB2, "TimeLGE" => time_LG2, "TimePE" => time_PE2)
# perfs2 = hcat(set2, perfs2)
# CSV.write("Simulation_times_" * string(N2) * "_all.csv", perfs2)
