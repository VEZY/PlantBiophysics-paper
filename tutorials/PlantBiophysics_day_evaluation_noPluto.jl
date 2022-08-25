using PlantBiophysics
using DataFrames
using CSV
using Statistics
using MonteCarloMeasurements
using CairoMakie, Colors
using Dates
unsafe_comparisons(true)
constants = Constants()
CairoMakie.activate!()

# Defining data path and choosing day and tree displayed
datapath = "../data/"
date = "14/11/2001"
tree = 3

# Reading file
file = datapath * "/TumbarumbaGasex_Spot_Medlyn.csv"
df = read_licor6400(file)
filter!(x -> x.Date == date && x.Tree == tree, df)

# Fitting paramaters (photosynthesis and stomatal conductance)
df.Asim .= df.Esim .= df.Gssim .= df.Dlsim .= df.Tlsim .= 0.0 ± 0.0
VcMaxRef, JMaxRef, RdRef, TPURef = PlantBiophysics.fit(Fvcb, df)
g0, g1 = PlantBiophysics.fit(Medlyn, df)

# Adding dummy uncertainty to the meteorological data:
df_uncertain = select(
    df,
    :T => (x -> x ± 0.1) => :T,
    :T => (x -> 40.0 ± 10.0) => :Wind,
    :P => (x -> x ± (0.001 * x)) => :P,
    :Rh => (x -> x ± 0.01) => :Rh,
    :Cₐ => (x -> x ± 10.0) => :Cₐ,
)

# Instantiate the weather data:
meteo = Weather(df_uncertain)
# Instantiate the models:
leaf = ModelList(
    energy_balance=Monteith(maxiter=100),
    photosynthesis=Fvcb(VcMaxRef=VcMaxRef, JMaxRef=JMaxRef, RdRef=RdRef, TPURef=TPURef, θ=0.9),
    stomatal_conductance=Medlyn(g0, g1),
    status=(
        Rₛ=(df.PPFD ± (0.1 * df.PPFD)) / 4.57, # not / 0.48 because it is in the chamber, the source is only PAR
        sky_fraction=1.0,
        PPFD=df.PPFD ± (0.1 * df.PPFD),
        d=Particles(Uniform(0.01, 0.10))
    ),
    type_promotion=Dict(Float64 => Particles{Float64,2000}),
    variables_check=false
)

# Make the simulation:
energy_balance!(leaf, meteo)

# Extract the outputs:
df.Asim .= leaf[:A]
@. df.Esim = leaf[:λE] / (meteo[:λ] * constants.Mₕ₂ₒ) * 1000.0
df.Gssim .= leaf[:Gₛ]
df.Dlsim .= leaf[:Dₗ]
df.Tlsim .= leaf[:Tₗ]

# Plotting

begin
    # line_color = Colors.RGB(36 / 255, 36 / 255, 64 / 255)
    line_color = Colors.RGB(100 / 255, 130 / 255, 150 / 255)
    error_color = Colors.RGBA(100 / 255, 130 / 255, 150 / 255, 0.5)
    point_color = Colors.RGB(253 / 255, 100 / 255, 103 / 255)
    point_fill = Colors.RGBA(253 / 255, 100 / 255, 103 / 255, 0.5)

    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    size_inches = (6.7, 5)
    size_pt = 72 .* size_inches
    fig = Figure(
        font=noto_sans,
        # dpi=300,
        resolution=size_pt,
        fontsize=10,
        xminorgridstyle=true
    )

    ax = Axis(
        fig[1, 1],
        # title=string("Tree N°", tree, ", ", date),
        xlabel="Time (HH:MM)",
        ylabel="Leaf temperature (°C)",
        titlesize=10,
    )

    ax.xticks = (eachindex(df.Tlsim), Dates.format.(df.Time, dateformat"HH:MM"))
    dat = scatter!(
        ax,
        df.Tₗ,
        color=point_fill,
        markersize=12,
        label="Measurement",
        strokecolor=point_color,
        strokewidth=3
    )

    sim = lines!(
        ax,
        eachindex(df.Tlsim),
        pmean.(df.Tlsim),
        color=line_color,
        linewidth=3,
        label="Simulation"
    )

    err = errorbars!(
        ax,
        eachindex(df.Tlsim),
        pmean.(df.Tlsim),
        pstd.(df.Tlsim),
        pstd.(df.Tlsim),
        color=error_color,
        whiskerwidth=8,
        linewidth=2.5,
        label="95% conf. interval"
    )

    Legend(fig[2, 1], ax, orientation=:horizontal, framevisible=false, padding=0.0)
    fig
end

save("figure_day.png", fig, px_per_unit=3)
