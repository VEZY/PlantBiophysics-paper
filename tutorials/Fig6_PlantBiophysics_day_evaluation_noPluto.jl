using PlantBiophysics, PlantSimEngine, PlantMeteo
using DataFrames, CSV
using Statistics
using MonteCarloMeasurements
using CairoMakie, Colors
using Dates
using Revise
using Downloads
includet("./functions.jl")

unsafe_comparisons(true)
constants = Constants()
CairoMakie.activate!()

# Defining data path and choosing day and tree displayed
date = "14/11/2001"
tree = 3

# Reading spot measurements:
df = let
    df_ = read_licor6400(Downloads.download("https://figshare.com/ndownloader/files/3402638"))
    df_.Asim .= df_.Esim .= df_.Gssim .= df_.Dlsim .= df_.Tlsim .= 0.0 ± 0.0
    transform!(df_, [:Date, :Time] => ((x, y) -> Date.(x, dateformat"d/m/Y") .+ y) => :date)

    filter!(
        row -> row.Date == date && row.Tree == tree,# && row["Age Class"] == 1,
        df_,
    )
    df_
end;

# Reading A-Ci curves:
df_curve_leaf = let
    df_ = read_licor6400(Downloads.download("https://figshare.com/ndownloader/files/3402635"))
    df_.Asim .= df_.Esim .= df_.Gssim .= df_.Dlsim .= df_.Tlsim .= 0.0 ± 0.0
    transform!(
        df_,
        :Date =>
            (x -> Dates.format.(Date.(x, dateformat"Y/m/d"), dateformat"d/m/Y")) =>
                :Date,
        [:Date, :Time] => ((x, y) -> Date.(x, dateformat"Y/m/d") .+ y) => :date,
    )

    filter!(row -> row.Date == date && row["Leaf Age"] == 1, df_)
    df_
end;

# Select the columns related to the meteorology and add some uncertainty to the measurement:
meteo = let
    meteo_df = select(
        df,
        :date,
        :T => (x -> x ± 0.1) => :T,
        :T => (x -> 40.0 ± 10.0) => :Wind,
        :P => (x -> x ± (0.001 * x)) => :P,
        :Rh => (x -> x ± 0.01) => :Rh,
        :Cₐ => (x -> x ± 10.0) => :Cₐ,
    )
    Weather(meteo_df, (site="Tumbarumba",))
end

# Fitting parameters (photosynthesis and stomatal conductance)
VcMaxRef, JMaxRef, RdRef, TPURef = PlantBiophysics.fit(Fvcb, df_curve_leaf)
g0, g1 = PlantBiophysics.fit(Medlyn, df)

df_sim = let
    leaf = ModelList(
        energy_balance=Monteith(maxiter=100),
        photosynthesis=Fvcb(VcMaxRef=VcMaxRef, JMaxRef=JMaxRef, RdRef=RdRef, TPURef=TPURef),
        # stomatal_conductance=Medlyn(g0, g1),
        stomatal_conductance=ForcedGs(0.00001, 0.0),
        status=(
            Rₛ=(df.PPFD ± (0.1 * df.PPFD)) / 4.57, # not / 0.48 because it is in the chamber, the source is only PAR
            sky_fraction=1.0,
            PPFD=df.PPFD ± (0.1 * df.PPFD),
            d=Particles(Uniform(0.01, 0.10)),
            Gₛ=df.Gₛ ± 0.0,
        ),
        type_promotion=Dict(Float64 => Particles{Float64,2000}),
        variables_check=false,
    )

    # Make the simulation:
    energy_balance!(leaf, meteo)

    # Extract the outputs:
    df_sim = select(
        DataFrame(leaf),
        :A => :Asim,
        :λE => (x -> x ./ (meteo[:λ] .* constants.Mₕ₂ₒ) .* 1000.0) => :Esim,
        :Gₛ => :Gssim,
        :Dₗ => :Dlsim,
        :Tₗ => :Tlsim,
        :PPFD,
    )

    df_sim
end

# Plotting
begin
    noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
    x_labs = (eachindex(df.Tlsim), Dates.format.(df.Time, dateformat"HH:MM"))
    # size_inches = (6.7, 5)
    size_inches = (10, 7)
    size_pt = 72 .* size_inches
    fig = Figure(font=noto_sans, resolution=size_pt, fontsize=12, xminorgridstyle=true)

    ga = fig[1, 1] = GridLayout()

    axDl = Axis(ga[1, 1], ylabel="Dₗ (kPa)")
    plot_var(axDl, df.Dₗ, df_sim.Dlsim)
    # axDl = Axis(ga[1, 1], ylabel="Gₛ (mol m⁻² s⁻¹)")
    # plot_var(axDl, df.Gₛ, df_sim.Gssim)

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

save("out/figure_day.png", fig, px_per_unit=3)
