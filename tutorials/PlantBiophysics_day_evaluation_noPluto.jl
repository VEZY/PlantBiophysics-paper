using PlantBiophysics
using DataFrames
using CSV
using Statistics
using Plots
using MonteCarloMeasurements
using Makie
using CairoMakie
using Dates
unsafe_comparisons(true)
constants = Constants()
CairoMakie.activate!()

# Defining data path and choosing day and tree displayed
datapath="/Users/simon/Code/PlantBiophysics-paper/data/"
date="14/11/2001"
tree=3

# Reading file
file = datapath*"/TumbarumbaGasex_Spot_Medlyn.csv"
df = read_licor6400(file)
day = df[(df.Date.==date)&(day.Tree.==tree),:]

# Fitting paramaters (photosynthesis and stomatal conductance)
day.Asim .= day.Esim .= day.Gssim .= day.Dlsim .= day.Tlsim .= 0. ± 0.
VcMaxRef,JMaxRef,RdRef,TPURef = PlantBiophysics.fit(Fvcb,day)
g0,g1 = PlantBiophysics.fit(Medlyn,day)

# Simulating
for i in 1:size(day,1)
    meteo = Atmosphere(T = day.T[i] ± 0.1, Wind = 40. ± 10., P = day.P[i] ± 0.001*day.P[i], Rh = day.Rh[i] ± 0.01 ,Cₐ= day.Cₐ[i] ± 10.)
    leaf = LeafModels(energy = Monteith(maxiter=100),
        photosynthesis = Fvcb(VcMaxRef=VcMaxRef,JMaxRef=JMaxRef,RdRef=RdRef,TPURef=TPURef,θ=0.9),
    stomatal_conductance = Medlyn(g0, g1),
    Rₛ = (day.PPFD[i] ± 0.1*day.PPFD[i])/(4.57), sky_fraction =1. ± 0.1, 
    PPFD = day.PPFD[i] ± 0.1*day.PPFD[i], d = 0.01 .. 0.10)
    energy_balance!(leaf,meteo)
    day.Asim[i] = leaf[:A]
    day.Esim[i] = leaf[:λE]/(meteo.λ * constants.Mₕ₂ₒ)*1000.
    day.Gssim[i] = leaf[:Gₛ]
    day.Dlsim[i] = leaf[:Dₗ]
    day.Tlsim[i] = leaf[:Tₗ] 
end

# Plotting
noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
fig = Makie.Figure(backgroundcolor = RGBf(1, 1, 1),resolution = (1000, 500),font=noto_sans,dpi=600,size=(1000,1000),xminorgridstyle=true)
#sideinfo1 = Label(fig[2,0], "Simulations", rotation = pi/2,labelsize=25)
#sideinfo2 = Label(fig[2,0], "Observations")
ax = Axis(fig[1,1],title="Leaf temperature in day "*date*" and tree "*string(tree),xlabel="Time (hh:mm:ss)",ylabel="Leaf temperature (°C)")
ax.xticks = ([1,2,3,4,5,6],string.(day.Time))
dat = Makie.scatter!(ax,day.Tₗ,color=:black,marker='∆',markersize=20)
err = errorbars!(ax,[1,2,3,4,5,6],pmean.(day.Tlsim), pstd.(day.Tlsim), pstd.(day.Tlsim),whiskerwidth = 15)
sim = lines!(ax,[1,2,3,4,5,6],pmean.(day.Tlsim),color=:black)
Makie.ylims!(ax,(14.5,27.))
axislegend(ax, [dat,sim,err], ["data","simulations","95% confidence interval"], position = :rb,
    orientation = :vertical,labelsize=13)
fig

save(path*"../tutorials/figure_day.pdf", fig, px_per_unit = 2) # output size = 1600 x 1200 pixels