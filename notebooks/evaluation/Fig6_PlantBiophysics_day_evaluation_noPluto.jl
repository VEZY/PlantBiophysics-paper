using PlantBiophysics, PlantSimEngine, PlantMeteo
import PlantBiophysics.gs_closure
using DataFrames, CSV, Downloads
using Statistics
using MonteCarloMeasurements
import MonteCarloMeasurements: ±
using CairoMakie, Colors
using Dates
unsafe_comparisons(true)
constants = Constants();

include(joinpath(@__DIR__, "functions.jl"))

df_Medlyn = let
    df_ = load_medlyn_data(type="spot"; abs=0.85)
    df_.Asim .= df_.Esim .= df_.Gssim .= df_.Dlsim .= df_.Tlsim .= 0.0 ± 0.0
    transform!(
        df_,
        [:Date, :Time] => ((x, y) -> Date.(x, dateformat"d/m/Y") .+ y) => :date
    )
    df_
end;

df_curves = let
    df_ = load_medlyn_data(; abs=0.85)
    df_.Asim .= df_.Esim .= df_.Gssim .= df_.Dlsim .= df_.Tlsim .= 0.0 ± 0.0
    transform!(
        df_,
        :Date => (x -> Dates.format.(Date.(x, dateformat"Y/m/d"), dateformat"d/m/Y")) => :Date,
        [:Date, :Time] => ((x, y) -> Date.(x, dateformat"Y/m/d") .+ y) => :date
    )
    df_
end;

dates_tree = Dict(tree => unique(filter(x -> x.Tree == tree, df_Medlyn).Date) for tree in unique(df_Medlyn.Tree))
tree = 3
date = dates_tree[tree][1]

df = filter(
    row -> row.Date == date && row.Tree == tree,# && row["Age Class"] == 1,
    df_Medlyn
)

df_curve_leaf = filter(
    row -> row.Date == date && row["Leaf Age"] == 1,
    df_curves
)

meteo = let
    meteo_df = select(
        df,
        :date,
        :T => (x -> x ± 0.1) => :T,
        :T => (x -> 40.0 ± 10.0) => :Wind,
        :P => (x -> x ± (0.001 * x)) => :P,
        :Rh => (x -> x ± 0.01) => :Rh,
        :Cₐ => (x -> x ± 10.0) => :Cₐ
    )
    Weather(meteo_df, (site="Tumbarumba",))
end;

# Fitting parameters

VcMaxRef, JMaxRef, RdRef, TPURef = PlantSimEngine.fit(Fvcb, df_curve_leaf)

# `g0` and `g1` are not fitted because the dataset does not present a Gₛ~VPD curve, and the snap measurements are already used for validation. See the script for an evaluation with a forcing of Gₛ.

# Simulation

function plot_var(ax, var_meas, var_sim)
    line_color = Colors.RGB(([67, 101, 139] ./ 255)...)
    error_color = line_color
    point_color = Colors.RGB(253 / 255, 100 / 255, 103 / 255)
    point_fill = Colors.RGBA(253 / 255, 100 / 255, 103 / 255, 0.5)

    scatter!(
        ax,
        var_meas,
        color=point_fill,
        markersize=12,
        label="Measurement",
        strokecolor=point_color,
        strokewidth=3
    )

    lines!(
        ax,
        pmean.(var_sim),
        color=line_color,
        linewidth=2.5,
        label="Simulation ± 95% confidence interval"
    )

    errorbars!(
        ax,
        eachindex(var_sim),
        pmean.(var_sim),
        pstd.(var_sim),
        pstd.(var_sim),
        color=error_color,
        whiskerwidth=5,
        linewidth=2
    )
end

function nRMSE(obs, sim; digits=2)
    return round(sqrt(sum((obs .- sim) .^ 2) / length(obs)) / (findmax(obs)[1] - findmin(obs)[1]), digits=digits)
end

function EF(obs, sim, digits=2)
    SSres = sum((obs - sim) .^ 2)
    SStot = sum((obs .- mean(obs)) .^ 2)
    return round(1 - SSres / SStot, digits=digits)
end

function Bias(obs, sim, digits=4)
    return round(mean(sim .- obs), digits=digits)
end

function nBias(obs, sim; digits=2)
    return round(mean((sim .- obs)) / (findmax(obs)[1] - findmin(obs)[1]), digits=digits)
end

"""
    ForcedGs()

A stomatal conductance model that forces the stomatal conductance to the value of `Gₛ` in the status.

It usually is used to force the stomatal conductance to the value measured in the chamber.

"""
struct ForcedGs <: PlantBiophysics.AbstractStomatal_ConductanceModel
    g0
end

ForcedGs() = ForcedGs(0.0)

# We implement a method for gs_closure as it is used in the photosynthesis model:
function PlantBiophysics.gs_closure(::ForcedGs, models, status, meteo=missing, constants=nothing, extra=nothing)
    # first iteration, we take measured Gs as a proxy:
    if status.A < 1e-9
        status.Gₛ
    else
        # Then we compute it using A from the previous iteration:
        status.A / (status.Gₛ - models.stomatal_conductance.g0)
    end
end;

# We implement the model as a method for run!:
function PlantSimEngine.run!(
    ::ForcedGs,
    models,
    status,
    meteo::M,
    constants=Constants(),
    extra=nothing,
) where {M<:PlantMeteo.AbstractAtmosphere}
    status.Gₛ
end

# We also implement a method for run! with a gs_closure (this is called from FvCB):
function PlantSimEngine.run!(::ForcedGs, models, status, gs_closure, extra=nothing)
    status.Gₛ
end

# Now we declared the inputs:
function PlantSimEngine.inputs_(::ForcedGs)
    (Gₛ=-Inf,)
end

# And the outputs:
function PlantSimEngine.outputs_(::ForcedGs)
    (Gₛ=-Inf,)
end

# ╔═╡ a1be516a-3b10-4c4a-a219-7e34d3e8edd7
df_sim = let

    leaf = ModelList(
        energy_balance=Monteith(maxiter=100),
        photosynthesis=Fvcb(VcMaxRef=VcMaxRef, JMaxRef=JMaxRef, RdRef=RdRef, TPURef=TPURef),
        stomatal_conductance=Medlyn(0.01, 3.42),
        status=(
            Ra_SW_f=(df.aPPFD ± (0.1 * df.aPPFD)) / 4.57, # not / 0.48 because it is in the chamber, the source is only PAR
            sky_fraction=1.0,
            aPPFD=df.aPPFD ± (0.1 * df.aPPFD),
            d=Particles(Uniform(0.01, 0.10))
        ),
        type_promotion=Dict(Float64 => Particles{Float64,2000}),
        variables_check=false
    )

    # Make the simulation:
    sim = run!(leaf, meteo)

    # Extract the outputs:
    df_sim =
        select(
            DataFrame(sim),
            :A => :Asim,
            :λE => (x -> x ./ (meteo[:λ] .* constants.Mₕ₂ₒ) .* 1000.0) => :Esim,
            :Gₛ => :Gssim,
            :Dₗ => :Dlsim,
            :Tₗ => :Tlsim,
            :aPPFD
        )

    df_sim
end

let
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    x_labs = (eachindex(df.Tlsim), Dates.format.(df.Time, dateformat"HH:MM"))
    # size_inches = (6.7, 5)
    size_inches = (10, 7)
    size_pt = 72 .* size_inches
    fig = Figure(
        font=noto_sans,
        size=size_pt,
        fontsize=12,
        xminorgridstyle=true
    )

    ga = fig[1, 1] = GridLayout()

    # axDl = Axis(ga[1, 1], ylabel="Dₗ (kPa)")
    # plot_var(axDl, df.Dₗ, df_sim.Dlsim)
    axDl = Axis(ga[1, 1], ylabel="Gₛ (mol m⁻² s⁻¹)")
    plot_var(axDl, df.Gₛ, df_sim.Gssim)

    axTl = Axis(ga[1, 2], ylabel="Tₗ (°C)")
    plot_var(axTl, df.Tₗ, df_sim.Tlsim)

    hidexdecorations!(axDl, grid=false)
    hidexdecorations!(axTl, grid=false)

    axA = Axis(ga[2, 1], xlabel="Time (HH:MM)", ylabel="A (μmol m⁻² s⁻¹)")
    axA.xticks = deepcopy(x_labs)
    plot_var(axA, df.A, df_sim.Asim)

    axE = Axis(ga[2, 2], xlabel="Time (HH:MM)", ylabel="Tr (mol m⁻² s⁻¹)")
    axE.xticks = deepcopy(x_labs)
    plot_var(axE, df.Trmmol, df_sim.Esim)

    rowgap!(ga, 10)
    Legend(fig[2, 1], axDl, orientation=:horizontal, framevisible=false, padding=0.0)
    fig
end

# ╔═╡ 612c00e7-6abd-4023-a27a-c75c95a73e38
df_sim_forcedGs = let
    leaf = ModelList(
        Monteith(maxiter=100),
        Fvcb(VcMaxRef=VcMaxRef, JMaxRef=JMaxRef, RdRef=RdRef, TPURef=TPURef),
        ForcedGs(),
        status=(
            Ra_SW_f=(df.aPPFD ± (0.1 * df.aPPFD)) / 4.57,
            sky_fraction=1.0,
            aPPFD=df.aPPFD ± (0.1 * df.aPPFD),
            d=Particles(Uniform(0.01, 0.10)),
            Gₛ=df.Gₛ ± 0.0,
        ),
        type_promotion=Dict(Float64 => Particles{Float64,2000}),
        variables_check=false,
    )

    # Make the simulation:
    sim = run!(leaf, meteo)

    # Extract the outputs:
    df_ = select(
        DataFrame(sim),
        :A => :Asim,
        :λE => (x -> x ./ (meteo[:λ] .* constants.Mₕ₂ₒ) .* 1000.0) => :Esim,
        :Gₛ => :Gssim,
        :Dₗ => :Dlsim,
        :Tₗ => :Tlsim,
        :aPPFD,
    )
    df_
end

# ╔═╡ 108e2a1e-d078-43c6-a490-b4e93d5f1da1
fig = let
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    x_labs = (eachindex(df.Tlsim), Dates.format.(df.Time, dateformat"HH:MM"))
    # size_inches = (6.7, 5)
    size_inches = (10, 7)
    size_pt = 72 .* size_inches
    fig = Figure(font=noto_sans, size=size_pt, fontsize=12, xminorgridstyle=true)

    ga = fig[1, 1] = GridLayout()

    axDl = Axis(ga[1, 1], ylabel="Dₗ (kPa)")
    plot_var(axDl, df.Dₗ, df_sim_forcedGs.Dlsim)
    # axDl = Axis(ga[1, 1], ylabel="Gₛ (mol m⁻² s⁻¹)")
    # plot_var(axDl, df.Gₛ, df_sim.Gssim)

    axTl = Axis(ga[1, 2], ylabel="Tₗ (°C)")
    plot_var(axTl, df.Tₗ, df_sim_forcedGs.Tlsim)

    hidexdecorations!(axDl, grid=false)
    hidexdecorations!(axTl, grid=false)

    axA = Axis(ga[2, 1], xlabel="Time (HH:MM)", ylabel="A (μmol m⁻² s⁻¹)")
    axA.xticks = deepcopy(x_labs)
    plot_var(axA, df.A, df_sim_forcedGs.Asim)

    axE = Axis(ga[2, 2], xlabel="Time (HH:MM)", ylabel="Tr (mol m⁻² s⁻¹)")
    axE.xticks = deepcopy(x_labs)
    plot_var(axE, df.Trmmol, df_sim_forcedGs.Esim)

    rowgap!(ga, 10)
    Legend(fig[2, 1], axDl, orientation=:horizontal, framevisible=false, padding=0.0)
    fig
end

# ╔═╡ 518deea7-f7a5-4c41-87b1-2ecade7eadb4
save("figure_day.png", fig, px_per_unit=3);

# ╔═╡ a46a3ae6-cf66-4071-b546-931ba1c8abb7
let
    df_vec = []
    vars = [:Dₗ, :Tₗ, :A, :Tr]
    for var in [df.Dₗ => df_sim_forcedGs.Dlsim, df.Tₗ => df_sim_forcedGs.Tlsim, df.A => df_sim_forcedGs.Asim, df.Trmmol => df_sim_forcedGs.Esim]
        var_ = popfirst!(vars)
        var_vec = Any[:variable=>var_]
        for fn in [RMSE, nRMSE, EF, Bias, nBias]
            push!(var_vec, Symbol(fn) => fn(var.first, var.second))
        end
        push!(df_vec, (; var_vec...))
    end

    DataFrame(df_vec)
end