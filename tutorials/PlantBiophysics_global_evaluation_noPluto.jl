########################################################################
# Global simulation of ACis curves from Medlyn data
# Simon Treillou, 2022
########################################################################

# using Plots
using Revise
using CSV, Statistics, DataFrames, Downloads, Dates
using AlgebraOfGraphics, CairoMakie, Colors
using PlantBiophysics, RCall, LeafGasExchange, Cropbox
using MonteCarloMeasurements

constants = Constants()
unsafe_comparisons(true)

includet("functions.jl")

# Do you want to save the simulations ?
saving_simulations = false

# Do you want to save the figure ?
saving_figure = false

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
df = read_licor6400(Downloads.download("https://figshare.com/ndownloader/files/3402635"))

transform!(
    df,
    [:Date, :Time] => ((x, y) -> Date.(x, dateformat"Y/m/d") .+ y) => :Date
)

df.VcMaxRef .= df.JMaxRef .= df.RdRef .= df.TPURef .= df.g0 .= df.g1 .= df.Tᵣ .= 0.0
df.VcMaxRefPE .= df.JMaxRefPE .= df.RdRefPE .= df.TPURefPE .= df.g0PE .= df.g1PE .= 0.0
for i in unique(df.Curve)
    dfi = filter(x -> x.Curve == i, df)
    sort!(dfi, :Cᵢ)

    g0, g1 = PlantBiophysics.fit(Medlyn, dfi)
    df.g0[df.Curve.==i] .= g0
    df.g1[df.Curve.==i] .= g1

    filter!(x -> x.PPFD > 1400.0, dfi)

    VcMaxRef, JMaxRef, RdRef, TPURef, Tᵣ = PlantBiophysics.fit(Fvcb, dfi)
    df.VcMaxRef[df.Curve.==i] .= VcMaxRef
    df.JMaxRef[df.Curve.==i] .= JMaxRef
    df.RdRef[df.Curve.==i] .= RdRef
    df.TPURef[df.Curve.==i] .= TPURef
    df.Tᵣ[df.Curve.==i] .= Tᵣ

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
        fitmethod = "default",
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

########################################################################
# Simulation
########################################################################
# Define characteristic length (here we have the Licor6400 area, we deduce a
# characteristic length in m
d = sqrt(df.Area[1]) / 100#sqrt(df.Area[1]/pi) / 100
Wind = 20.0
@rput Wind
@rput d
Leaf_abs = 0.86 # default from plantecophys
emissivity = 0.95 # default from plantecophys


########################################################################
# PlantBiophysics.jl

df.AsimPB .= df.EsimPB .= df.TlsimPB .= df.GssimPB .= 0.0
for i in unique(df.Curve)
    dfi = filter(x -> x.Curve == i, df)

    cols = fieldnames(Atmosphere)
    dfiMeteo = select(dfi, names(dfi, x -> Symbol(x) in cols))
    dfiMeteo.Wind .= Wind
    # Note that as we only use A-Ci curves, there is no NIR in the Licor6400
    dfiMeteo.Ri_SW_f .= dfi.PPFD .* Leaf_abs ./ (4.57)
    meteo = Weather(dfiMeteo)

    leaf = ModelList(
        energy_balance=Monteith(
            aₛₕ=2,
            aₛᵥ=1,
            ε=emissivity, # Matching the value in plantecophys (https://github.com/RemkoDuursma/plantecophys/blob/c9749828041f10ca47c6691436678e0a5632cfb8/R/LeafEnergyBalance.R#L112)
            maxiter=100
        ),
        photosynthesis=Fvcb(
            Tᵣ=dfi.Tᵣ[1],
            VcMaxRef=dfi.VcMaxRef[1],
            JMaxRef=dfi.JMaxRef[1],
            RdRef=dfi.RdRef[1],
            TPURef=dfi.TPURef[1]
        ),
        stomatal_conductance=Medlyn(dfi.g0[1], dfi.g1[1]),
        status=(
            Rₛ=meteo[:Ri_SW_f],
            sky_fraction=1.0,
            PPFD=dfi.PPFD,
            d=d
        )
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
    configs = []
    for i in 1:length(dfi.T)
        config = :Weather => (
            PFD=dfi.PPFD[i],
            CO2=dfi.Cₐ[i],
            RH=dfi.Rh[i] * 100,
            T_air=dfi.T[i],
            wind=Wind,
            P_air=dfi.P[i],
            g0=dfi.g0[i],
            g1=1.57 * dfi.g1[i],
            Vc25=dfi.VcMaxRef[i],
            Jm25=dfi.JMaxRef[i],
            Rd25=dfi.RdRef[i],
            Tp25=dfi.TPURef[i],
            Ear=46.39,
            Haj=29.68,
            w=d,
            d=d,
            EaVc=58.55,
            ϵ=emissivity,
            Dh=21.5,
            α_s=1 - Leaf_abs
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
    VcMaxRef = dfi$VcMaxRef[1]
    JMaxRef = dfi$JMaxRef[1]
    TPURef = dfi$TPURef[1]
    RdRef = dfi$RdRef[1]
    g0 = dfi$g0[1]
    g1 = dfi$g1[1]
    res = PhotosynEB(
        Tair = dfi$T,
        Wind = Wind,
        VPD=dfi$VPD,
        Wleaf = d,
        Ca = Ca,
        # StomatalRatio = 1, LeafAbs = 0.86, # default values
        gsmodel = "BBOpti",
        g0 = g0, g1 = g1,
        EaV = 58550.0,EdVC = 2e+05, delsC = 629.26,
        EaJ = 29680.0,EdVJ = 2e+05,delsJ = 631.88,
        alpha = 0.24,theta = 0.7, Jmax = JMaxRef,
        Vcmax = VcMaxRef, TPU = TPURef,Rd = RdRef,
        RH = dfi$Rh*100,
        PPFD=dfi$PPFD,
        Patm = dfi$P,gk=0.,
        Tcorrect = FALSE
    )
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
    CSV.write("out/Medlyn_ACis_simulations.csv", df)
end

########################################################################
# Stacking the results
########################################################################

meas =
    stack(
        select(
            df,
            [:Date => :Date, :Cₐ => :Cₐ, :A => :A, :Trmmol => :E, :Tₗ => :Tl, :Gₛ => :Gs]
        ),
        [:A, :E, :Tl, :Gs],
        [:Date, :Cₐ],
        value_name=:measured
    )

sim_PB =
    stack(
        select(
            df,
            [:Date => :Date, :AsimPB => :A, :EsimPB => :E, :TlsimPB => :Tl, :GssimPB => :Gs]
        ),
        [:A, :E, :Tl, :Gs],
        :Date,
        value_name=:simulated
    )
sim_PB.origin .= "PlantBiophysics.jl"

sim_LG =
    stack(
        select(
            df,
            [:Date => :Date, :AsimLG => :A, :EsimLG => :E, :TlsimLG => :Tl, :GssimLG => :Gs]
        ),
        [:A, :E, :Tl, :Gs],
        :Date,
        value_name=:simulated
    )
sim_LG.origin .= "LeafGasExchange.jl"

sim_PE = stack(
    select(
        df,
        [:Date => :Date, :AsimPE => :A, :EsimPE => :E, :TlsimPE => :Tl, :GssimPE => :Gs]
    ),
    [:A, :E, :Tl, :Gs],
    :Date,
    value_name=:simulated
)
sim_PE.origin .= "plantecophys"

df_all = vcat(sim_PB, sim_LG, sim_PE)
df_res = leftjoin(df_all, meas, on=[:Date, :variable])

########################################################################
# Filtering
########################################################################
# We need to filter-out some points because they are not reliable for the stomatal conductance

filter!(x -> x.Cₐ > 150, df_res)
filter!(x -> x.Cₐ > 150, df)

########################################################################
# Statistics
########################################################################

stats =
    combine(
        groupby(df_res, [:variable, :origin], sort=true),
        [:measured, :simulated] => ((x, y) -> RMSE(x, y)) => :RMSE,
        [:measured, :simulated] => ((x, y) -> nRMSE(x, y)) => :nRMSE,
        [:measured, :simulated] => ((x, y) -> Bias(x, y)) => :Bias,
        [:measured, :simulated] => ((x, y) -> nBias(x, y)) => :nBias,
        [:measured, :simulated] => ((x, y) -> EF(x, y)) => :EF
    )

########################################################################
# Plotting with Makie.jl
########################################################################

function rgb(r, g, b)
    return RGB(r / 255, g / 255, b / 255)
end

function rgb(r, g, b, a)
    return RGBA(r / 255, g / 255, b / 255, a)
end

begin
    transparency_col = 1.0
    transparency_fill = 0.5
    color_pb = rgb(148, 180, 159, transparency_col)
    color_lg = rgb(223, 120, 97, transparency_col)
    color_pe = rgb(118, 84, 154, transparency_col)
    fill_pb = rgb(148, 180, 159, transparency_fill)
    fill_lg = rgb(223, 120, 97, transparency_fill)
    fill_pe = rgb(118, 84, 154, transparency_fill)
    stw = 1.5 # strokewidth
    ms = 7 # markersize
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    legend_lab_size = 10
    size_inches = (10, 10)
    size_pt = 72 .* size_inches
    fig = Figure(
        font=noto_sans,
        resolution=size_pt,
        fontsize=12,
        xminorgridstyle=true
    )

    sideinfo1 = Label(fig[1:2, 1], "Simulations", rotation=pi / 2, textsize=12)
    sideinfo2 = Label(fig[3, 2:3], "Observations", textsize=12)

    # Assimilation
    axa = Axis(fig[1, 2], title="(a) Aₙ: Net CO₂ assimilation", aspect=1)
    xlims!(-10.0, 50.0)
    ylims!(-10.0, 50.0)

    abline!(axa, 0, 1, color=(:grey, 0.4), linewidth=4)

    LG = scatter!(
        axa, df.A, df.AsimLG,
        color=fill_lg,
        markersize=ms,
        strokecolor=color_lg,
        strokewidth=stw
    )
    PE = scatter!(axa, df.A, df.AsimPE,
        color=fill_pe,
        markersize=ms,
        strokecolor=color_pe,
        strokewidth=stw
    )
    PB = scatter!(axa, df.A, df.AsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )
    axislegend(
        axa,
        [PB, PE, LG],
        [
            "nRMSE: " * string(filter(x -> x.variable == "A" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            "nRMSE: " * string(filter(x -> x.variable == "A" && x.origin == "plantecophys", stats).nRMSE[1]),
            "nRMSE: " * string(filter(x -> x.variable == "A" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
        ],
        "",
        position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )

    # Transpiration
    axb = Axis(fig[1, 3], title="(b) E: Transpiration rate", aspect=1)
    xlims!(-0.5, 10.0)
    ylims!(-0.5, 10.0)

    abline!(axb, 0, 1, color=(:grey, 0.4), linewidth=4)

    LG = scatter!(axb, df.Trmmol, df.EsimLG,
        color=fill_lg,
        markersize=ms,
        strokecolor=color_lg,
        strokewidth=stw
    )
    PE = scatter!(axb, df.Trmmol, df.EsimPE,
        color=fill_pe,
        markersize=ms,
        strokecolor=color_pe,
        strokewidth=stw
    )
    PB = scatter!(axb, df.Trmmol, df.EsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )
    axislegend(
        axb,
        [PB, PE, LG],
        [
            "nRMSE: " * string(filter(x -> x.variable == "E" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            "nRMSE: " * string(filter(x -> x.variable == "E" && x.origin == "plantecophys", stats).nRMSE[1]),
            "nRMSE: " * string(filter(x -> x.variable == "E" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
        ],
        "", position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )

    # Stomatal conductance
    axc = Axis(
        fig[2, 2],
        title="(c) Gₛ: CO₂ stomatal conductance",
        aspect=1
    )
    xlims!(-0.05, 0.85)
    ylims!(-0.05, 0.85)

    abline!(axc, 0, 1, color=(:grey, 0.4), linewidth=4)

    LG = scatter!(axc, df.Gₛ, df.GssimLG,
        color=fill_lg,
        markersize=ms,
        strokecolor=color_lg,
        strokewidth=stw,
        label="LeafGasExchange.jl"
    )
    PE = scatter!(axc, df.Gₛ, df.GssimPE,
        color=fill_pe,
        markersize=ms,
        strokecolor=color_pe,
        label="plantecophys",
        strokewidth=stw,
    )
    PB = scatter!(axc, df.Gₛ, df.GssimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw,
        label="PlantBiophysics.jl"
    )

    axislegend(
        axc,
        [PB, PE, LG],
        [
            "nRMSE: " * string(filter(x -> x.variable == "Gs" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            "nRMSE: " * string(filter(x -> x.variable == "Gs" && x.origin == "plantecophys", stats).nRMSE[1]),
            "nRMSE: " * string(filter(x -> x.variable == "Gs" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
        ],
        "",
        position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )


    # Leaf temperature
    axd = Axis(fig[2, 3], title="(d) Tₗ: Leaf temperature", aspect=1)
    abline!(axd, 0, 1, color=(:grey, 0.4), linewidth=4)

    LG = scatter!(axd, df.Tₗ, df.TlsimLG,
        color=fill_lg,
        markersize=ms,
        strokecolor=color_lg,
        strokewidth=stw
    )
    PE = scatter!(axd, df.Tₗ, df.TlsimPE,
        color=fill_pe,
        markersize=ms,
        strokecolor=color_pe,
        strokewidth=stw
    )
    PB = scatter!(axd, df.Tₗ, df.TlsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )

    xlims!(10.0, 36.0)
    ylims!(10.0, 36.0)

    axislegend(
        axd, [PB, PE, LG],
        [
            "nRMSE: " * string(filter(x -> x.variable == "Tl" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            "nRMSE: " * string(filter(x -> x.variable == "Tl" && x.origin == "plantecophys", stats).nRMSE[1]),
            "nRMSE: " * string(filter(x -> x.variable == "Tl" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
        ],
        "", position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )

    Legend(
        fig[4, 1:end],
        axc,
        orientation=:horizontal,
        framevisible=false,
        padding=0.0
    )
    fig
end

save("out/figure_global_simulation.png", fig, px_per_unit=3)
