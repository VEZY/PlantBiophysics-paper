########################################################################
# Global simulation of ACis curves from Medlyn data
# Simon Treillou, 2022
########################################################################

# using Plots
# using Revise  # Comment out Revise if not available
using CSV, Statistics, DataFrames, Downloads, Dates
using CairoMakie, Colors
using PlantBiophysics, PlantSimEngine, PlantMeteo, RCall, Cropbox, LeafGasExchange
using MonteCarloMeasurements

constants = Constants()
unsafe_comparisons(true)
include(joinpath(@__DIR__, "functions.jl"))

# Do you want to save the simulations ?
saving_simulations = false

# Do you want to save the figure ?
saving_figure = false

Leaf_abs = 0.85 # Von Caemmerer et al. (2009)
emissivity = 0.95 # default from plantecophys

# Loading R packages
R"""
if (!require('plantecophys')) install.packages('plantecophys', repos = "https://cloud.r-project.org"); library('plantecophys')
if (!require('dplyr')) install.packages('dplyr', repos = "https://cloud.r-project.org"); library('dplyr')
if (!require('readr')) install.packages('readr', repos = "https://cloud.r-project.org"); library('readr')
if (!require('microbenchmark')) install.packages('microbenchmark', repos = "https://cloud.r-project.org"); library('microbenchmark')
"""

########################################################################
# Parameters : fit them or download already fitted parameters
########################################################################

df = load_medlyn_data(abs=Leaf_abs)
rename!(df, :PARi => :PPFD) # Rename Ttop to Tmin for consistency

transform!(df, [:Date, :Time] => ((x, y) -> Date.(x, dateformat"Y/m/d") .+ y) => :Date)

# Initialize fitting and simulation result columns if they don't exist
required_columns = [
    :VcMaxRef, :JMaxRef, :RdRef, :TPURef, :g0, :g1, :Tᵣ,
    :VcMaxRefPE, :JMaxRefPE, :RdRefPE, :TPURefPE, :g0PE, :g1PE,
    :AsimPB, :EsimPB, :TlsimPB, :GssimPB,
    :AsimLG, :EsimLG, :TlsimLG, :GssimLG,
    :AsimPE, :EsimPE, :TlsimPE, :GssimPE, :PEfailed
]

for col in required_columns
    if !(col in names(df))
        df[!, col] .= 0.0
    end
end

for i in unique(df.Curve)
    dfi = filter(x -> x.Curve == i, df)
    sort!(dfi, :Cᵢ)

    g0, g1 = PlantSimEngine.fit(Medlyn, dfi)
    df.g0[df.Curve.==i] .= g0
    df.g1[df.Curve.==i] .= g1

    filter!(x -> x.PPFD > 1400.0, dfi)

    VcMaxRef, JMaxRef, RdRef, TPURef, Tᵣ = collect(PlantSimEngine.fit(Fvcb, dfi, α=0.425, θ=0.7))
    df.VcMaxRef[df.Curve.==i] .= VcMaxRef
    df.JMaxRef[df.Curve.==i] .= JMaxRef
    df.RdRef[df.Curve.==i] .= RdRef
    df.TPURef[df.Curve.==i] .= TPURef
    df.Tᵣ[df.Curve.==i] .= Tᵣ

    dfiPE = dfi[:, 3:end]
    P = mean(dfi.P)
    # rename!(dfiPE, :Tₗ => :Tleaf, :Cᵢ => :Ci, :Cₐ => :Ca, :Gₛ => :gs)
    @rput P
    @rput dfiPE
    R"""
    fit = fitaci(
        dfiPE,
        varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "aPPFD"),
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
        PPFD = dfiPE$aPPFD, # Note: in the code PPFD is used as aPPFD, there is no correction by the leaf absorptance, so it is indeed aPPFD
        Tleaf = dfiPE$Tleaf,
        alpha = 0.425,
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
    # rename!(dfiPE, :Tₗ => :Tleaf, :Cᵢ => :Ci, :Cₐ => :Ca, :Gₛ => :gs, :Dₗ => :VpdL)
    rename!(dfiPE, :Cₐ => :Ca, :Gₛ => :gs)
    transform!(dfiPE, :gs => (x -> PlantBiophysics.gsc_to_gsw.(x)) => :gs) # Convert Gs to conductance to water vapor
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

########################################################################
# PlantBiophysics.jl

atm_cols = keys(Atmosphere(T=25.0, Rh=0.5, Wind=10.0))
df.AsimPB .= df.EsimPB .= df.TlsimPB .= df.GssimPB .= 0.0
for i in unique(df.Curve) # i = 1
    dfi = filter(x -> x.Curve == i, df)
    dfiMeteo = select(dfi, names(dfi, x -> Symbol(x) in atm_cols))
    dfiMeteo.Wind .= Wind
    # There is no NIR in the Licor6400 light, only PAR, so the shortwave radiation is only PAR
    # incident PAR to the leaf (we correct by the leaf absorptance because it is aPPFD in input)
    Ra_SW_f = dfi.aPPFD ./ 4.57
    dfiMeteo.check .= false
    meteo = Weather(dfiMeteo)

    leaf = ModelList(
        energy_balance=Monteith(
            aₛₕ=2,
            aₛᵥ=1,
            ε=emissivity, # Matching the value in plantecophys (https://github.com/RemkoDuursma/plantecophys/blob/c9749828041f10ca47c6691436678e0a5632cfb8/R/LeafEnergyBalance.R#L112)
            maxiter=100,
        ),
        photosynthesis=Fvcb(
            Tᵣ=dfi.Tᵣ[1],
            VcMaxRef=dfi.VcMaxRef[1],
            JMaxRef=dfi.JMaxRef[1],
            RdRef=dfi.RdRef[1],
            TPURef=dfi.TPURef[1],
            α=0.425,
            θ=0.7,
        ),
        stomatal_conductance=Medlyn(dfi.g0[1], dfi.g1[1]),
        status=(Ra_SW_f=Ra_SW_f, sky_fraction=1.0, aPPFD=dfi.aPPFD, d=d),
    )

    sim = run!(leaf, meteo)
    df.AsimPB[df.Curve.==i, :] = sim.A
    df.EsimPB[df.Curve.==i, :] = sim.λE ./ (meteo[:λ] * constants.Mₕ₂ₒ) * 1000.0
    df.TlsimPB[df.Curve.==i, :] = sim.Tₗ
    df.GssimPB[df.Curve.==i, :] = sim.Gₛ
end

########################################################################
# LeafGasExchange.jl

df.AsimLG .= df.EsimLG .= df.TlsimLG .= df.GssimLG .= 0.0
for i in unique(df.Curve)
    dfi = filter(x -> x.Curve == i, df)
    configs = []
    for i = 1:length(dfi.T)
        config =
            :Weather => (
                PFD=dfi.PPFD[i], # Note that we use PPFD here, not aPPFD because we pass the absorptance to the model
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
                δ=1 - Leaf_abs,
                θ=0.7,
                f=0.15 # corresponds to the α value used above, because α = (1-f)/2
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
        StomatalRatio = 1, LeafAbs = $Leaf_abs,
        gsmodel = "BBOpti",
        g0 = g0, g1 = g1,
        EaV = 58550.0,EdVC = 2e+05, delsC = 629.26,
        EaJ = 29680.0,EdVJ = 2e+05,delsJ = 631.88,
        alpha = 0.425,theta = 0.7, Jmax = JMaxRef,
        Vcmax = VcMaxRef, TPU = TPURef,Rd = RdRef,
        RH = dfi$Rh*100,
        PPFD=dfi$PPFD, # Note that we use PPFD here, not aPPFD because we pass the absorptance to the model
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
    CSV.write("Medlyn_ACis_simulations.csv", df)
end

########################################################################
# Stacking the results
########################################################################

meas = stack(
    select(
        df,
        [:Date => :Date, :Cₐ => :Cₐ, :A => :A, :Trmmol => :E, :Tₗ => :Tl, :Gₛ => :Gs],
    ),
    [:A, :E, :Tl, :Gs],
    [:Date, :Cₐ],
    value_name=:measured,
)

sim_PB = stack(
    select(
        df,
        [:Date => :Date, :AsimPB => :A, :EsimPB => :E, :TlsimPB => :Tl, :GssimPB => :Gs],
    ),
    [:A, :E, :Tl, :Gs],
    :Date,
    value_name=:simulated,
)
sim_PB.origin .= "PlantBiophysics.jl"

sim_LG = stack(
    select(
        df,
        [:Date => :Date, :AsimLG => :A, :EsimLG => :E, :TlsimLG => :Tl, :GssimLG => :Gs],
    ),
    [:A, :E, :Tl, :Gs],
    :Date,
    value_name=:simulated,
)
sim_LG.origin .= "LeafGasExchange.jl"

sim_PE = stack(
    select(
        df,
        [:Date => :Date, :AsimPE => :A, :EsimPE => :E, :TlsimPE => :Tl, :GssimPE => :Gs],
    ),
    [:A, :E, :Tl, :Gs],
    :Date,
    value_name=:simulated,
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

if saving_simulations
    CSV.write("df_res.csv", df_res)
end

########################################################################
# Statistics
########################################################################

stats = combine(
    groupby(df_res, [:variable, :origin], sort=true),
    [:measured, :simulated] => ((x, y) -> RMSE(x, y)) => :RMSE,
    [:measured, :simulated] => ((x, y) -> nRMSE(x, y)) => :nRMSE,
    [:measured, :simulated] => ((x, y) -> Bias(x, y)) => :Bias,
    [:measured, :simulated] => ((x, y) -> nBias(x, y)) => :nBias,
    [:measured, :simulated] => ((x, y) -> EF(x, y)) => :EF,
)

if saving_simulations
    CSV.write("statistics.csv", stats)
end

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
    col_pb = (244, 124, 124)
    col_lg = (93, 174, 139)
    col_pe = (112, 161, 215)
    transparency_col = 0.6
    transparency_fill = 0.4

    color_pb = rgb(col_pb..., transparency_col)
    color_lg = rgb(col_lg..., transparency_col)
    color_pe = rgb(col_pe..., transparency_col)
    fill_pb = rgb(col_pb..., transparency_fill)
    fill_lg = rgb(col_lg..., transparency_fill)
    fill_pe = rgb(col_pe..., transparency_fill)

    stw = 1.5 # strokewidth
    ms = 7 # markersize
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    legend_lab_size = 10
    size_inches = (8.5, 8.5)
    size_pt = 72 .* size_inches
    fig = Figure(
        font=noto_sans,
        size=size_pt,
        fontsize=12,
        xminorgridstyle=true,
        backgroundcolor=:transparent,
    )

    sideinfo1 = Label(fig[1:2, 1], "Simulations", rotation=pi / 2, fontsize=18)
    sideinfo2 = Label(fig[3, 2:3], "Observations", fontsize=18)

    # Assimilation
    axa = Axis(
        fig[1, 2],
        title="a) Net CO₂ assimilation (Aₙ)",
        titlealign=:left,
        backgroundcolor=:transparent,
    )

    xlims!(-10.0, 50.0)
    ylims!(-10.0, 50.0)

    ablines!(axa, 0, 1, color=(:grey, 0.4), linewidth=4)

    LG = scatter!(
        axa,
        df.A,
        df.AsimLG,
        color=fill_lg,
        markersize=ms,
        strokecolor=color_lg,
        strokewidth=stw,
        label=string(
            filter(
                x -> x.variable == "A" && x.origin == "LeafGasExchange.jl",
                stats,
            ).nRMSE[1],
        )
    )
    PE = scatter!(
        axa,
        df.A,
        df.AsimPE,
        color=fill_pe,
        markersize=ms,
        strokecolor=color_pe,
        strokewidth=stw,
        label=string(
            filter(x -> x.variable == "A" && x.origin == "plantecophys", stats).nRMSE[1],
        )
    )
    PB = scatter!(
        axa,
        df.A,
        df.AsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw,
        label=string(
            filter(
                x -> x.variable == "A" && x.origin == "PlantBiophysics.jl",
                stats,
            ).nRMSE[1],
        )
    )
    axislegend(axa, "nRMSE:", position=:rb, orientation=:vertical, labelsize=legend_lab_size, padding=0.0, framevisible=false,)

    # Transpiration
    axb = Axis(
        fig[1, 3],
        title="b) Transpiration rate (E)",
        titlealign=:left,
        backgroundcolor=:transparent,
    )
    xlims!(-0.5, 10.0)
    ylims!(-0.5, 10.0)

    ablines!(axb, 0, 1, color=(:grey, 0.4), linewidth=4)

    LG = scatter!(
        axb,
        df.Trmmol,
        df.EsimLG,
        color=fill_lg,
        markersize=ms,
        strokecolor=color_lg,
        strokewidth=stw,
        label=string(
            filter(
                x -> x.variable == "E" && x.origin == "LeafGasExchange.jl",
                stats,
            ).nRMSE[1],
        )
    )
    PE = scatter!(
        axb,
        df.Trmmol,
        df.EsimPE,
        color=fill_pe,
        markersize=ms,
        strokecolor=color_pe,
        strokewidth=stw,
        label=string(
            filter(x -> x.variable == "E" && x.origin == "plantecophys", stats).nRMSE[1],
        )
    )
    PB = scatter!(
        axb,
        df.Trmmol,
        df.EsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw,
        label=string(
            filter(
                x -> x.variable == "E" && x.origin == "PlantBiophysics.jl",
                stats,
            ).nRMSE[1],
        )
    )
    axislegend(axb, "nRMSE:", position=:rb, orientation=:vertical, labelsize=legend_lab_size, padding=0.0, framevisible=false,)

    # Stomatal conductance
    axc = Axis(
        fig[2, 2],
        title="c) CO₂ stomatal conductance (Gₛ)",
        titlealign=:left,
        backgroundcolor=:transparent,
    )
    xlims!(-0.05, 0.85)
    ylims!(-0.05, 0.85)

    ablines!(axc, 0, 1, color=(:grey, 0.4), linewidth=4)

    LG = scatter!(
        axc,
        df.Gₛ,
        df.GssimLG,
        color=fill_lg,
        markersize=ms,
        strokecolor=color_lg,
        strokewidth=stw,
        label=string(
            filter(
                x -> x.variable == "Gs" && x.origin == "LeafGasExchange.jl",
                stats,
            ).nRMSE[1],
        )
    )
    PE = scatter!(
        axc,
        df.Gₛ,
        df.GssimPE,
        color=fill_pe,
        markersize=ms,
        strokecolor=color_pe,
        strokewidth=stw,
        label=string(
            filter(x -> x.variable == "Gs" && x.origin == "plantecophys", stats).nRMSE[1],
        ),
    )
    PB = scatter!(
        axc,
        df.Gₛ,
        df.GssimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw,
        label=string(
            filter(
                x -> x.variable == "Gs" && x.origin == "PlantBiophysics.jl",
                stats,
            ).nRMSE[1],
        )
    )

    axislegend(axc, "nRMSE:", position=:rb, orientation=:vertical, labelsize=legend_lab_size, padding=0.0, framevisible=false,)

    # Leaf temperature
    axd = Axis(
        fig[2, 3],
        title="d) Leaf temperature (Tₗ)",
        titlealign=:left,
        backgroundcolor=:transparent,
    )
    ablines!(axd, 0, 1, color=(:grey, 0.4), linewidth=4)

    LG = scatter!(
        axd,
        df.Tₗ,
        df.TlsimLG,
        color=fill_lg,
        markersize=ms,
        strokecolor=color_lg,
        strokewidth=stw,
        label="LeafGasExchange.jl",
    )
    PE = scatter!(
        axd,
        df.Tₗ,
        df.TlsimPE,
        color=fill_pe,
        markersize=ms,
        strokecolor=color_pe,
        strokewidth=stw,
        label="plantecophys",
    )
    PB = scatter!(
        axd,
        df.Tₗ,
        df.TlsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw,
        label="PlantBiophysics.jl",
    )

    xlims!(10.0, 36.0)
    ylims!(10.0, 36.0)

    axislegend(
        axd,
        [LG, PE, PB,],
        [
            string(
                filter(
                    x -> x.variable == "Tl" && x.origin == "LeafGasExchange.jl",
                    stats,
                ).nRMSE[1],
            ),
            string(
                filter(x -> x.variable == "Tl" && x.origin == "plantecophys", stats).nRMSE[1],
            ),
            string(
                filter(
                    x -> x.variable == "Tl" && x.origin == "PlantBiophysics.jl",
                    stats,
                ).nRMSE[1],
            ),],
        "nRMSE:",
        position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false,
    )

    Legend(fig[4, 1:end], axd, orientation=:horizontal, framevisible=false, padding=0.0)

    # colgap!(fig.layout, 5.0)
    # rowgap!(fig.layout, 0.0)

    # colsize!(fig.layout, 1, Aspect(1, 0.005))
    rowsize!(fig.layout, 3, Aspect(1, 0.005))
    rowsize!(fig.layout, 4, Aspect(1, 0.005))
    colgap!(fig.layout, 1, Relative(0.005))
    rowgap!(fig.layout, 1, Relative(0.02))
    rowgap!(fig.layout, 2, Relative(0.03))
    # colsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)
    fig
end

save("figure_global_simulation.png", fig, px_per_unit=6)
