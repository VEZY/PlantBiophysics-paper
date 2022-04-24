using Statistics
using LaTeXStrings
using Revise
using Cropbox
using LeafGasExchange
using RCall
using DataFrames
using Plots
using CSV

# Parameters : 
N    = 100000        # Length of the random imput set
PPFD = 1500.      # Absorbed Photosynthetic Photon Flux Density (μmol m-2 s-1)
@rput PPFD

# Create the ranges of input parameters
Rn   = range(0,500,length=100000)
Ta   = range(-10,50,length=100000)
Wind = range(0,20,length=100000)
P    = range(90,101,length=100000)
Rh   = range(0.1,1.,length=100000)
Ca   = range(360,900,length=100000)
skyF = range(0.,1.,length=100000)
d    = range(0.001,0.5,length=100000)
Jmax = range(200.,300.,length=100000)
Vmax = range(150.,250.,length=100000)
Rd   = range(0.3,2.,length=100000)
TPU  = range(5.,20.,length=100000)
g0   = range(0.001, 2., length=100000)
g1   = range(0.5, 15., length=100000)
vars = hcat([Ta,Wind,P,Rh,Ca,Jmax,Vmax,Rd,Rn,skyF,d,TPU,g0,g1])
print("Ranges defined.")

# Create the randomly drawn input set
set   = [rand.(vars) for i in 1:N]
set   = reshape(vcat(set...),(length(set[1]),length(set)))'
name = ["T","Wind","P","Rh","Ca","JMaxRef","VcMaxRef","RdRef","Rs","sky_fraction","d","TPURef","g0","g1"]
set   = DataFrame(set,name)

######################################################################################################
# BENCHMARKING PLANTECOPHYS
######################################################################################################

R"""
library(plantecophys)
library(dplyr)
library(readr)
library(microbenchmark)
time_PE = c()
"""

for i in 1:N
    vpd = 0.61375 * exp((17.502 * 25.) / (25. + 240.97)) - set.Rh[i] * 0.61375 * exp((17.502 * set.T[i]) / (set.T[i] + 240.97))
    @rput vpd
    
    df=DataFrame(set[i,:])
    @rput df
    R"""
    function_EB <- function(df,vpd) {
    res = PhotosynEB(Tair = df$T, VPD = vpd,Wind = df$Wind,
                    Wleaf = df$d,Ca = df$Ca,  StomatalRatio = 1,
                    LeafAbs = df$sky_fraction,gsmodel = "BBOpti",g0 = df$g0, g1 = df$g1,
                    alpha = 0.24,theta = 0.7, Jmax = df$JMaxRef, 
                    Vcmax = df$VcMaxRef, TPU = df$TPURef,Rd = df$RdRef,
                    RH = df$Rh*100,PPFD=PPFD,
                    Patm = df$P)
        return(res)
    }
    """
    
    R"""
    resBenchmark <- microbenchmark(function_EB(df,VPD), times=1)
    time_PE = append(time_PE,resBenchmark$time)
    """
end
@rget time_PE 
time_PE .*= 10e-9

######################################################################################################
# BENCHMARKING LEAFGASEXCHANGE.JL
######################################################################################################

time_LG = []
for i in 1:N
    config = :Weather => (
        PFD = PPFD,
        CO2 = set.Ca[i],
        RH = set.Rh[i] * 100,
        T_air = set.T[i],
        wind = set.Wind[i],
        P_air = set.P[i],
        g0 = set.g0[i],
        g1 = set.g1[i],
        Vcmax = set.VcMaxRef[i],
        Jmax = set.JMaxRef[i],
        Rd = set.RdRef[i],
        TPU = set.TPURef[i]
    )
    function LGtest()
        simulate(ModelC3MD;config=(config),nounit=true)
    end
    resLG = @elapsed LGtest()
    push!(time_LG,resLG)
end


######################################################################################################
# BENCHMARKING PLANTBIOPHYSICS.JL
######################################################################################################

time_PB = []
for i in 1:N
    leaf = LeafModels(energy = Monteith(),
                photosynthesis = Fvcb(VcMaxRef=set.VcMaxRef[i],JMaxRef=set.JMaxRef[i],RdRef=set.RdRef[i],TPURef=set.TPURef[i]),
                stomatal_conductance = Medlyn(set.g0[i], set.g1[i]),
                Rₛ = set.Rs[i], sky_fraction =set.sky_fraction[i], 
                PPFD = PPFD, d = set.d[i])
    meteo = Atmosphere(T = set.T[i], Wind = set.Wind[i], P = set.P[i], Rh = set.Rh[i],Cₐ= set.Ca[i])
    function PBtest()
        energy_balance!(leaf,meteo)
    end
    resPB = @elapsed PBtest()
    push!(time_PB,resPB)
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
    m   = mean(df)
    med = median(df)
    std = Statistics.std(df)
    min = findmin(df)[1]
    max = findmax(df)[1]
    return stat_res(m,med,std,min,max)
end

statsPB = basic_stat(time_PB)
statsPE = basic_stat(time_PE)
statsLG = basic_stat(time_LG)

factorPE = mean(time_PE ./ time_PB)
factorLG = mean(time_LG ./ time_PB)

######################################################################################################
# PLOTTING
######################################################################################################

Plots.plot(layout=grid(3,1),xminorgrid=true) # Using `Plots.plot` as `plot` function is defined in Cropbox
interval=(1e-6,0.5) # x-axis

# Ribbon plots (i.e. plot the standard deviation)
eps=1e-6 # limits for PlantBiophysics.jl standard deviation: sometimes the lower limit of the interval is negative and there are problems with logscale
Plots.plot!([(max(eps,statsPB.mean-statsPB.stddev),0.),(statsPB.mean+statsPB.stddev,0.),(statsPB.mean+statsPB.stddev,0.6),(max(eps,statsPB.mean-statsPB.stddev),0.6),(max(eps,statsPB.mean-statsPB.stddev),0.)],seriestype=:shape,fillcolor=:orange,alpha=0.3,linecolor=:blue,linewidth=0.,sp=1,label="[μ-σ, μ+σ]")
Plots.plot!([(statsPE.mean-statsPE.stddev,0.),(statsPE.mean+statsPE.stddev,0.),(statsPE.mean+statsPE.stddev,0.6),(statsPE.mean-statsPE.stddev,0.6),(statsPE.mean-statsPE.stddev,0.)],seriestype=:shape,fillcolor=:orange,alpha=0.3,linecolor=:blue,linewidth=0.,sp=2,label="")
Plots.plot!([(statsLG.mean-statsLG.stddev,0.),(statsLG.mean+statsLG.stddev,0.),(statsLG.mean+statsLG.stddev,0.6),(statsLG.mean-statsLG.stddev,0.6),(statsLG.mean-statsLG.stddev,0.)],seriestype=:shape,fillcolor=:orange,alpha=0.3,linecolor=:blue,linewidth=0.,sp=3,label="")

# Histogram densities plots
histogram!(time_PB[time_PB.<1e-4],sp=1,xaxis=(:log10, interval),normalize=:probability,legend=false,bins=100,dpi=300,color=:lightblue,label="")
histogram!(time_PE[time_PE.<1e-0],sp=2,xaxis=(:log10, interval),normalize=:probability,bins=30,dpi=300,color=:lightblue,label="")
histogram!(time_LG[time_LG.<1e-1],sp=3,xaxis=(:log10, interval),normalize=:probability,bins=30,dpi=300,color=:lightblue,label="")

# Mean plotting
Plots.plot!(statsPB.mean*[1,1],[0.,0.6],sp=1,linewidth=2,linestyle=:dot,color=:red,label="Mean μ",legend=:bottomright)
Plots.plot!(statsPE.mean*[1,1],[0.,0.6],sp=2,linewidth=2,linestyle=:dot,color=:red,label="")
Plots.plot!(statsLG.mean*[1,1],[0.,0.6],sp=3,linewidth=2,linestyle=:dot,color=:red,label="")

# Annotating (a), (b), (c)
annotate!(0.9*interval[2], 0.62, Plots.text("(a) PlantBiophysics.jl", :black, :right, 9),sp=1)
annotate!(0.9*interval[2], 0.68, Plots.text("(b) plantecophys", :black, :right, 9),sp=2)
annotate!(0.9*interval[2], 0.68, Plots.text("(c) LeafGasExchange.jl", :black, :right, 9),sp=3)
Plots.xlabel!("time (s)",xguidefontsize=10,sp=3)
Plots.ylabel!("density",xguidefontsize=10,sp=2)
Plots.xticks!([1e-6,1e-5,1e-4,1e-3,1e-2,1e-1])

######################################################################################################
# SAVING
######################################################################################################
savefig("fig")

perfs = DataFrame("TimePB" => time_PB, "TimeLGE" => time_LG,"TimePE" => time_PE)
perfs = hcat(set,perfs)
CSV.write("Simulations_temps_"*string(N)*".csv",perfs)