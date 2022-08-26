########################################################################
# Global simulation of ACis curves from Medlyn data
# Simon Treillou, 2022
########################################################################
using DataFrames
using Plots
using RCall
using CSV
using PlantBiophysics
using Statistics
using Revise
constants = Constants()
using CairoMakie
using Dates
using MonteCarloMeasurements
unsafe_comparisons(true)
using LeafGasExchange
using Cropbox
using Makie

# Add absolute path for PlantBiophysics-ComplCode
path = "/Users/simon/Code/PlantBiophysics-paper/data/"

# Did you already have fitted the parameters ?
fitted = false

# Do you want to save the simulations ?
saving_simulations = false

# Do you want to save the figure ?
saving_figure = false

function RMSE(obs, sim)
    return sqrt(sum((obs .- sim) .^ 2) / length(obs))
end

function EF(obs, sim)
    SSres = sum((obs - sim) .^ 2)
    SStot = sum((obs .- mean(obs)) .^ 2)
    return 1 - SSres / SStot
end

# Loading R packages
R"""
library(plantecophys)
library(dplyr)
library(readr)
library(microbenchmark)
"""

########################################################################
# Parameters : fit them or download already fitted parameters
########################################################################
if !fitted
    df = read_licor6400(path * "TumbarumbaGasex_ACis_Medlyn.csv")
    df.Date = Date.(df.Date, Dates.DateFormat("yyyy/mm/dd"))

    df.VcMaxRef .= df.JMaxRef .= df.RdRef .= df.TPURef .= df.g0 .= df.g1 .= 0.0
    df.VcMaxRefPE .= df.JMaxRefPE .= df.RdRefPE .= df.TPURefPE .= df.g0PE .= df.g1PE .= 0.0
    for i in unique(df.Curve)
        dfi = filter(x -> x.Curve == i, df)
        sort!(dfi, :Cᵢ)

        g0, g1 = PlantBiophysics.fit(Medlyn, dfi)
        df.g0[df.Curve.==i] .= g0
        df.g1[df.Curve.==i] .= g1

        filter!(x -> x.PPFD > 1400.0, dfi)

        VcMaxRef, JMaxRef, RdRef, TPURef = PlantBiophysics.fit(Fvcb, dfi)
        df.VcMaxRef[df.Curve.==i] .= VcMaxRef
        df.JMaxRef[df.Curve.==i] .= JMaxRef
        df.RdRef[df.Curve.==i] .= RdRef
        df.TPURef[df.Curve.==i] .= TPURef

        dfiPE = dfi[:, 3:end]
        P = mean(dfi.P)
        rename!(dfiPE, :Tₗ => :Tleaf, :Cᵢ => :Ci, :Cₐ => :Ca, :Gₛ => :gs)
        @rput P
        @rput dfiPE
        R"""
        fit = fitaci(
            dfiPE,
            varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PPFD"),
            Tcorrect = TRUE,
            Patm = P,
            citransition = NULL,
            quiet = TRUE,
            startValgrid = TRUE,
            fitmethod = "bilinear",
            algorithm = "nls",
            fitTPU = TRUE,
            alphag = 0,
            useRd = FALSE,
            PPFD = dfiPE$PPFD,
            Tleaf = dfiPE$Tleaf,
            alpha = 0.24,
            theta = 0.7,
            gmeso = NULL,
            EaV = 58550.0,
            EdVC = 2e+05,
            delsC = 629.26,
            EaJ = 29680.0,
            EdVJ = 2e+05,
            delsJ = 631.88,
            #GammaStar = NULL,
            #Km = NULL,
            id = NULL)
        coefs = coef(fit)
        """
        @rget coefs
        dfi = filter(x -> x.Curve == i, df)
        dfiPE = dfi[:, 3:end]
        rename!(dfiPE, :Tₗ => :Tleaf, :Cᵢ => :Ci, :Cₐ => :Ca, :Gₛ => :gs, :Dₗ => :VpdL)
        filter!(x -> x.A ./ (x.Ca) > 0.0, dfiPE)
        @rput dfiPE
        R"""
        dfiPE$Rh = dfiPE$Rh*100
        myfit = fitBB(
            dfiPE,
            varnames = list(ALEAF = "A", GS = "gs", VPD = "VpdL", Ca = "Ca", RH = "Rh"),
            gsmodel = c("BBOpti"),
            fitg0 = TRUE,
            D0 = NULL
            )
        gP = coef(myfit)
        """
        @rget gP
        df.VcMaxRefPE[df.Curve.==i] .= coefs[1]
        df.JMaxRefPE[df.Curve.==i] .= coefs[2]
        df.RdRefPE[df.Curve.==i] .= coefs[3]
        if !ismissing(coefs[4])
            df.TPURefPE[df.Curve.==i] .= coefs[4]
        else
            df.TPURefPE[df.Curve.==i] .= 1000000.0
        end
        df.g0PE[df.Curve.==i] .= gP[1]
        df.g1PE[df.Curve.==i] .= gP[2]
    end
    CSV.write(path * "../tutorials/Medlyn_ACis_param.csv", df)

else
    df = CSV.File(path * "../tutorials/Medlyn_ACis_param.csv", dateformat="yyyy/mm/dd") |> DataFrame
end
########################################################################
# Simulation
########################################################################
# Define characteristic length (here we have the Licor6400 area, we deduce a
# characteristic length in m
d = sqrt(df.Area[1]) / 100#sqrt(df.Area[1]/pi) / 100
Wind = 20.0
@rput Wind
@rput d

########################################################################
# PlantBiophysics.jl

df.AsimPB .= df.EsimPB .= df.TlsimPB .= df.GssimPB .= 0.0
for i in unique(df.Curve)
    dfi = filter(x -> x.Curve == i, df)

    cols = fieldnames(Atmosphere)
    dfiMeteo = select(dfi, names(dfi, x -> Symbol(x) in cols))
    dfiMeteo.Wind .= Wind
    # Note that as we only use A-Ci curves, there is no NIR in the Licor6400
    dfiMeteo.Ri_SW_f .= dfi.PPFD ./ (4.57)
    meteo = Weather(dfiMeteo)

    leaf = LeafModels(
        energy=Monteith(aₛₕ=2, aₛᵥ=1, ε=0.955, maxiter=100),
        photosynthesis=Fvcb(VcMaxRef=dfi.VcMaxRef[1], JMaxRef=dfi.JMaxRef[1], RdRef=dfi.RdRef[1], TPURef=dfi.TPURef[1]),
        stomatal_conductance=Medlyn(dfi.g0[1], dfi.g1[1]),
        Rₛ=meteo[:Ri_SW_f],
        sky_fraction=1.0,
        PPFD=dfi.PPFD,
        d=d
    )

    energy_balance!(leaf, meteo)
    df.AsimPB[df.Curve.==i, :] = DataFrame(leaf).A
    df.EsimPB[df.Curve.==i, :] = DataFrame(leaf).λE ./ (meteo[:λ] * constants.Mₕ₂ₒ) * 1000
    df.TlsimPB[df.Curve.==i, :] = DataFrame(leaf).Tₗ
    df.GssimPB[df.Curve.==i, :] = DataFrame(leaf).Gₛ
end

########################################################################
# LeafGasExchange.jl

df.AsimLG .= df.EsimLG .= df.TlsimLG .= df.GssimLG .= 0.0
for i in unique(df.Curve)
    dfi = filter(x -> x.Curve == i, df)
    #dfi = filter(x->x.PPFD > 1400.,dfi)
    #dfi = filter(x->x.Cₐ > 150.,dfi)

    configs = []
    for i in 1:length(dfi.T)
        config = :Weather => (
            PFD=dfi.PPFD[i],
            CO2=dfi.Cₐ[i],
            RH=dfi.Rh[i] * 100,
            T_air=dfi.T[i],
            wind=Wind,
            P_air=dfi.P[i],
            g0=dfi.g0PE[i],
            g1=1.57 * dfi.g1PE[i],
            Vc25=dfi.VcMaxRefPE[i],
            Jm25=dfi.JMaxRefPE[i],
            Rd25=dfi.RdRefPE[i],
            Tp25=dfi.TPURefPE[i],
            Ear=46.39,
            Haj=29.68,
            w=d,
            d=d,
            EaVc=58.55,
            ϵ=0.955,
            Dh=21.5
        )
        push!(configs, config)
    end
    res = simulate(ModelC3MD; configs=configs, nounit=true)
    df.AsimLG[df.Curve.==i, :] = res.A_net
    df.EsimLG[df.Curve.==i, :] = res.E
    df.TlsimLG[df.Curve.==i, :] = res.T
    df.GssimLG[df.Curve.==i, :] = res.gsc
end

########################################################################
# plantecophys

R"""
A_sim = c()
E_sim = c()
Tl_sim = c()
Gs_sim = c()
failed =c()
"""

df.AsimPE .= df.EsimPE .= df.TlsimPE .= df.GssimPE .= df.PEfailed .= 0.0
for i in unique(df.Curve)
    dfi = filter(x -> x.Curve == i, df)
    dfi = dfi[:, 3:end]
    @rput dfi
    Ca = dfi.Cₐ
    @rput Ca
    R"""
    VcMaxRef = dfi$VcMaxRefPE[1]
    JMaxRef = dfi$JMaxRefPE[1]
    TPURef = dfi$TPURefPE[1]
    RdRef = dfi$RdRefPE[1]
    g0 = dfi$g0PE[1]
    g1 = dfi$g1PE[1]
    res = PhotosynEB(Tair = dfi$T,Wind = Wind,VPD=dfi$VPD,
             Wleaf = d,Ca = Ca,  StomatalRatio = 1,
             LeafAbs = 1.,gsmodel = "BBOpti",g0 = g0, g1 = g1,
             EaV = 58550.0,EdVC = 2e+05, delsC = 629.26,
             EaJ = 29680.0,EdVJ = 2e+05,delsJ = 631.88,
             alpha = 0.24,theta = 0.7, Jmax = JMaxRef,
             Vcmax = VcMaxRef, TPU = TPURef,Rd = RdRef,
             RH = dfi$Rh*100,PPFD=dfi$PPFD,
             Patm = dfi$P,gk=0.,
             Tcorrect = FALSE)
    A_sim = append(A_sim,res$ALEAF)
    Gs_sim = append(Gs_sim,res$GS)
    failed = append(failed,res$failed)
    Tl_sim = append(Tl_sim,res$Tleaf)
    E_sim = append(E_sim,res$ELEAF)
    """
end
@rget A_sim

@rget failed
@rget Gs_sim
@rget E_sim
@rget Tl_sim

df.AsimPE .= A_sim
df.EsimPE .= E_sim
df.TlsimPE .= Tl_sim
df.GssimPE .= Gs_sim
df.PEfailed .= failed

if saving_simulations
    CSV.write(path * "../tutorials/Medlyn_ACis_simulations.csv", df)
end


########################################################################
# Plotting with Makie.jl
########################################################################
CairoMakie.activate!(type="svg")
noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
fig = Makie.Figure(backgroundcolor=RGBf(1, 1, 1), resolution=(1000, 1000), font=noto_sans, dpi=600, size=(1000, 1000), xminorgridstyle=true)
sideinfo1 = Label(fig[1:2, 0], "Simulations", rotation=pi / 2, labelsize=25)
sideinfo2 = Label(fig[3, 2:3], "Observations")

########################################################################
# Assimilation
axa = Axis(fig[1, 2], title="(a) Net CO₂ assimilation Aₙ")
index = 1:length(df.PPFD)
index = ifelse.(df.PPFD .> 1400, index, 9999.0)
filter!(x -> x != 9999.0, index)
indexPE = 1:length(df.PPFD)
indexPE = ifelse.((df.PPFD .> 1400) .& (df.PEfailed .== 0.0), indexPE, 9999.0)
filter!(x -> x != 9999.0, indexPE)
RMSPB = round(RMSE(df.A[index], df.AsimPB[index]), digits=3)
RMSLG = round(RMSE(df.A[index], df.AsimLG[index]), digits=3)
RMSPE = round(RMSE(df.A[indexPE], df.AsimPE[indexPE]), digits=3)
EF(df.A[index], df.AsimPB[index])
EF(df.A[index], df.AsimLG[index])
EF(df.A[indexPE], df.AsimPE[indexPE])
Makie.lines!(axa, df.AsimPB[index], df.AsimPB[index], color=(:grey, 0.4), linewidth=4)
LG = Makie.scatter!(axa, df.A[index], df.AsimLG[index], color=(:black, 0.5), markersize=8, marker='□')
PE = Makie.scatter!(axa, df.A[indexPE], df.AsimPE[indexPE], color=(:blue, 0.5), markersize=8, marker='∆')
PB = Makie.scatter!(axa, df.A[index], df.AsimPB[index], color=(:red, 0.8), markersize=8, marker='o')
axislegend(axa, [PB, PE, LG], ["PlantBiophysics.jl (RMSE=" * string(RMSPB) * ")", "plantecophys (RMSE=" * string(RMSPE) * ")", "LeafGasExchange.jl (RMSE=" * string(RMSLG) * ")"], "", position=:rb,
    orientation=:vertical, labelsize=13)

########################################################################
# Transpiration
index = 1:length(df.PPFD)
index = ifelse.(df.Cₐ .> 140.0, index, 9999.0)
filter!(x -> x != 9999.0, index)
indexPE = 1:length(df.PPFD)
indexPE = ifelse.((df.Cₐ .> 140.0) .& (df.PEfailed .== 0.0), indexPE, 9999.0)
filter!(x -> x != 9999.0, indexPE)

axb = Axis(fig[1, 3], title="(b) Transpiration rate E")

RMSPB = round(RMSE(df.Trmmol[index], df.EsimPB[index]), digits=3)
RMSLG = round(RMSE(df.Trmmol[index], df.EsimLG[index]), digits=3)
RMSPE = round(RMSE(df.Trmmol[indexPE], df.EsimPE[indexPE]), digits=3)
EF(df.Trmmol[index], df.EsimPB[index])
EF(df.Trmmol[index], df.EsimLG[index])
EF(df.Trmmol[indexPE], df.EsimPE[indexPE])
Makie.lines!(axb, df.EsimPB[index], df.EsimPB[index], linewidth=4, color=(:grey, 0.4))
LG = Makie.scatter!(axb, df.Trmmol[index], df.EsimLG[index], color=(:black, 0.5), markersize=8, marker='□')
PE = Makie.scatter!(axb, df.Trmmol[indexPE], df.EsimPE[indexPE], color=(:blue, 0.5), markersize=8, marker='∆')
PB = Makie.scatter!(axb, df.Trmmol[index], df.EsimPB[index], color=(:red, 0.8), markersize=8, marker='o')
axislegend(axb, [PB, PE, LG], ["PlantBiophysics.jl (RMSE=" * string(RMSPB) * ")", "plantecophys (RMSE=" * string(RMSPE) * ")", "LeafGasExchange.jl (RMSE=" * string(RMSLG) * ")"], "", position=:rb,
    orientation=:vertical, labelsize=13)

########################################################################
# Leaf temperature
indexPE = 1:length(df.PPFD)
indexPE = ifelse.((df.PEfailed .== 0.0), indexPE, 9999.0)
filter!(x -> x != 9999.0, indexPE)

axd = Axis(fig[2, 3], title="(d) Leaf temperature Tₗ")

RMSPB = round(RMSE(df.Tₗ, df.TlsimPB), digits=3)
RMSLG = round(RMSE(df.Tₗ, df.TlsimLG), digits=3)
RMSPE = round(RMSE(df.Tₗ[indexPE], df.TlsimPE[indexPE]), digits=3)
EF(df.Tₗ[index], df.TlsimPB[index])
EF(df.Tₗ[index], df.TlsimLG[index])
EF(df.Tₗ[indexPE], df.TlsimPE[indexPE])
Makie.lines!(axd, df.TlsimPB, df.TlsimPB, linewidth=4, color=(:grey, 0.4))
LG = Makie.scatter!(axd, df.Tₗ, df.TlsimLG, color=(:black, 0.5), markersize=8, marker='□')
PE = Makie.scatter!(axd, df.Tₗ[indexPE], df.TlsimPE[indexPE], color=(:blue, 0.5), markersize=8, marker='∆')
PB = Makie.scatter!(axd, df.Tₗ, df.TlsimPB, color=(:red, 0.8), markersize=8, marker='o')
axislegend(axd, [PB, PE, LG], ["PlantBiophysics.jl (RMSE=" * string(RMSPB) * ")", "plantecophys (RMSE=" * string(RMSPE) * ")", "LeafGasExchange.jl (RMSE=" * string(RMSLG) * ")"], "", position=:rb,
    orientation=:vertical, labelsize=13)

########################################################################
# Stomatal conductance
index = 1:length(df.PPFD)
index = ifelse.(df.Cₐ .> 140.0, index, 9999.0)
filter!(x -> x != 9999.0, index)
indexPE = 1:length(df.PPFD)
indexPE = ifelse.((df.Cₐ .> 140.0) .& (df.PEfailed .== 0.0), indexPE, 9999.0)
filter!(x -> x != 9999.0, indexPE)
axc = Axis(fig[2, 2], title="(c) CO₂ stomatal conductance Gₛ")


RMSLG = round(RMSE(df.Gₛ[index], df.GssimLG[index]), digits=3)
RMSPB = round(RMSE(df.Gₛ[index], df.GssimPB[index]), digits=3)
RMSPE = round(RMSE(df.Gₛ[indexPE], df.GssimPE[indexPE]), digits=3)
EF(df.Gₛ[index], df.GssimPB[index])
EF(df.Gₛ[index], df.GssimLG[index])
EF(df.Gₛ[indexPE], df.GssimPE[indexPE])
Makie.lines!(axc, df.GssimPB, df.GssimPB, linewidth=4, color=(:grey, 0.4))
LG = Makie.scatter!(axc, df.Gₛ[index], df.GssimLG[index], color=(:black, 0.5), markersize=8, marker='□')
PE = Makie.scatter!(axc, df.Gₛ[indexPE], df.GssimPE[indexPE], color=(:blue, 0.5), markersize=8, marker='∆')
PB = Makie.scatter!(axc, df.Gₛ[index], df.GssimPB[index], color=(:red, 0.8), markersize=8, marker='o')
Makie.ylims!(axc, (0.0, 0.4))
Makie.xlims!(axc, (0.0, 0.4))
axislegend(axc, [PB, PE, LG], ["PlantBiophysics.jl (RMSE=" * string(RMSPB) * ")", "plantecophys (RMSE=" * string(RMSPE) * ")", "LeafGasExchange.jl (RMSE=" * string(RMSLG) * ")"], "", position=:rb,
    orientation=:vertical, labelsize=13)
fig

if saving_figure
    save(path * "../tutorials/figure_global_simulation.png", fig, px_per_unit=2) # output size = 1600 x 1200 pixels
end
