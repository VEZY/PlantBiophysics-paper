### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 74ae23d5-71bc-4ad8-acdf-0a38a1b129bd
begin
	using CSV, Statistics, DataFrames, Downloads, Dates
	using AlgebraOfGraphics, CairoMakie, Colors
	using PlantBiophysics
end

# ╔═╡ 1100f4ba-d45e-11ec-25f5-b79c9d3b9851
md"""
# Global evaluation of PlantBiophysics.jl: observations vs simulations

This Pluto notebook presents the computation of Fig. 3 from the scientific article. The notebook does not compute anything because it would imply a dependency on R (plantecophys) and Python (LeafGasExchange.jl), which is possible but not performant in Pluto an Github actions. Instead, we display the code and the results only. if you want to reproduce the results, execute the code provided here, or the code from the script provided [here](https://github.com/VEZY/PlantBiophysics-paper/blob/main/tutorials/Fig3_PlantBiophysics_global_evaluation_noPluto.jl).

## Importing the dependencies:

!!! note 
	Make sure to have R installed on your computer first.
"""

# ╔═╡ 1609ce3c-4066-45ca-acce-8dbf081bfb89
md"""
Loading the Julia packages:
```julia
using CSV, Statistics, DataFrames, Downloads, Dates
using AlgebraOfGraphics, CairoMakie, Colors
using PlantBiophysics, RCall, LeafGasExchange, Cropbox
using MonteCarloMeasurements
```
"""

# ╔═╡ aff7314d-7a8a-4a47-8a77-252e96377e9c
md"""
Loading the R package (note the use of the `R\"\"\"` macro here that sends the code to R):
```julia
R\"\"\"
library(plantecophys)
\"\"\"
```
"""

# ╔═╡ c2e2c3c8-288c-491a-a8c7-8adf4517103c
md"""
## Reading the data
"""

# ╔═╡ ab97a0a5-d155-4758-a1d4-70e3a0b8c796
md"""
The data comes from Medlyn et al. (2015), see [here](https://figshare.com/articles/dataset/Tumbarumba_Gas_Exchange/1538079?file=3402641) for more details.

"""

# ╔═╡ 3c5bcd7e-87fb-4aef-86d8-326921b4423c
df = let 
	df_ = read_licor6400(Downloads.download("https://figshare.com/ndownloader/files/3402635"))

	# Computing DateTime from Date and Time:
	transform!(
		df_, 
		[:Date, :Time] => ((x, y) -> Date.(x, dateformat"Y/m/d") .+ y) => :Date
	)

	# Initializing the columns:
	df_.VcMaxRef .= df_.JMaxRef .= df_.RdRef .= df_.TPURef .= df_.g0 .= df_.g1 .= df_.Tᵣ .= 0.0
	df_.AsimPB .= df_.EsimPB .= df_.TlsimPB .= df_.GssimPB .= 0.0
	df_.VcMaxRefPE .= df_.JMaxRefPE .= df_.RdRefPE .= df_.TPURefPE .= df_.g0PE .= df_.g1PE .= 0.0
	df_
end

# ╔═╡ e8961d19-85d4-484b-a029-359a9dbea800
md"""
## Fitting

Photosynthesis and stomatal conductance parameters are fitted using PlantBiophysics.jl here, but it is also possible to use plantecophys instead (see Julia the script for more info). Comparison results are close with both.
"""

# ╔═╡ f9be64a3-b788-48d8-a925-5954e532effb
df_fit = let df_ = df
	for i in unique(df.Curve)
	    dfi = filter(x -> x.Curve == i, df_)
	    sort!(dfi, :Cᵢ)
	
	    g0, g1 = PlantBiophysics.fit(Medlyn, dfi)
	    df_.g0[df.Curve.==i] .= g0
	    df_.g1[df.Curve.==i] .= g1
	
	    filter!(x -> x.PPFD > 1400.0, dfi)
	
	    VcMaxRef, JMaxRef, RdRef, TPURef, Tᵣ = PlantBiophysics.fit(Fvcb, dfi)
	    df_.VcMaxRef[df_.Curve.==i] .= VcMaxRef
	    df_.JMaxRef[df_.Curve.==i] .= JMaxRef
	    df_.RdRef[df_.Curve.==i] .= RdRef
	    df_.TPURef[df_.Curve.==i] .= TPURef
	    df_.Tᵣ[df_.Curve.==i] .= Tᵣ
	end
end

# ╔═╡ 3592d0a1-d658-4137-8ea8-e8ea9ff3dc91
md"""
## Simulation

We set the wind velocity to 20 m/s considering that the Licor-6400 chamber is well ventilated. The characteristic length is set to the square root of the chamber area in $m^2$. The leaf absorptance and emissivity are set to the default value from plantecophys.
"""

# ╔═╡ a65d0175-f6cf-4bc9-8a2f-fbde4d2692b6
begin
	d = sqrt(df.Area[1]) / 100 # Characteristic dimension
	Wind = 20.0 # Wind, in m/s
	Leaf_abs = 0.86 # default from plantecophys
	emissivity = 0.95 # default from plantecophys
end

# ╔═╡ 44db3882-e5b3-4427-b781-e5185e2d70af
df_PB = let df = df  
	constants = Constants()
	
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
	df
end

# ╔═╡ 92deadb1-537b-4c3d-b0d2-3fd3cba171e4
md""" 
### With plantecophys and LeafGasExchange.jl

Below is presented the code for running the exact same simulation with the two other packages.

!!! note 
	The code shown here for plantecophys and LeafGasExchange.jl is not executed, just displayed. If you want to reproduce the computation, please use the script instead.
"""

# ╔═╡ 6e3a7997-17dd-453e-879e-8f8d6b46ca2d
md"""
#### plantecophys
"""

# ╔═╡ 8fe5544c-8f61-4306-b424-01d62e25fbcc
md"""
```julia
begin
	if !simulated
				
		R\"\"\"
		A_sim = c()
		E_sim = c()
		Tl_sim = c()
		Gs_sim = c()
		failed =c()
		\"\"\"

		df.AsimPE .= df.EsimPE .= df.TlsimPE .= df.GssimPE .= df.PEfailed .= 0.
		for i in unique(df.Curve)
		    dfi = filter(x->x.Curve == i,df)
		    dfi = dfi[:,3:end]
		    @rput dfi
		    Ca = dfi.Cₐ
		    @rput Ca
		    R\"\"\"
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
		    \"\"\"
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
	end
end
```
"""

# ╔═╡ 96af7e3e-38ce-4a8a-b503-163e88ddfb8f
md"""
#### LeafGasExchange.jl
"""

# ╔═╡ 232de02a-9440-4df5-8fbc-1e2eb6dae24a
md"""
```julia
begin
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
end
```
"""

# ╔═╡ 56705c8b-58ac-4c2d-acc5-dcbb03f25bf0
md"""
## Stacking the results

Here we stack the results in a long-format DataFrame for computing the statistics:
"""

# ╔═╡ db357873-44b5-4110-a80b-98d3adf785e6
begin
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
	            df_PB,
	            [:Date => :Date, :AsimPB => :A, :EsimPB => :E, :TlsimPB => :Tl, :GssimPB => :Gs]
	        ),
	        [:A, :E, :Tl, :Gs],
	        :Date,
	        value_name=:simulated
	    )
	sim_PB.origin .= "PlantBiophysics.jl"
	sim_PB
end

# ╔═╡ 5ac831c1-cc9e-47e9-a11b-bd191a1d55a2
md"""
Same for the other packages, though the data is not available in the notebook:
"""

# ╔═╡ 7a09574b-b06e-4c36-9fea-7757c20fd763
md"""
```julia
begin 
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
end
```
"""

# ╔═╡ 3d4375ea-2020-4829-ba70-5406c921c943
md"""
And finally we can stack the results together:
"""

# ╔═╡ 53a0fa15-d042-4e00-84ea-f9caa5fc3c3e
md"""
```julia
begin 
	df_all = vcat(sim_PB, sim_LG, sim_PE)
	df_res = leftjoin(df_all, meas, on=[:Date, :variable])
end
```
"""

# ╔═╡ 9f4d46ab-cc75-47e1-aef1-b7ada0e01e29
md"""
Note that the real code executed here is the following, because neither `sim_LG` or `sim_PE` are available:
"""

# ╔═╡ 182787c9-5075-451f-bdc7-5a8add68d141
begin
df_all = sim_PB
df_res = leftjoin(df_all, meas, on=[:Date, :variable])
# Filtering out the results for very low Cₐ:
filter!(x -> x.Cₐ > 150, df_res)
filter!(x -> x.Cₐ > 150, df)
df_res
end

# ╔═╡ 2fabe562-5eb9-4ef1-897e-f6e380d1518c
md"""
## Statistics

Computing the statistics for each model:
"""

# ╔═╡ e8171ee6-36fa-4227-8f78-2b0bd8675cdf
md"""
## Plotting

 - Comparison of simulations and observations for $A_n$, $E$, $G_s$ and $T_l$.
"""

# ╔═╡ f180ca3d-490c-4a64-8608-4f2f5e57c41f
md"""
### Observations vs simulations (Fig. 3)
"""

# ╔═╡ dd192660-2c6d-46a2-a06a-f0c1d01a9598
md"""
*Figure 1. Measured (x) and simulated (y) net carbon assimilation (a), net transpiration rate (b), stomatal conductance for CO2 (c) and leaf temperature (d). Data comes from Eucalyptus delegatensis (Medlyn et al., 2015). All simulations were performed using a photosynthesis-stomatal conductance-energy balance coupled model with PlantBiophysics.jl. Grey line represents x=y. All simulations were done with Ca > 150 ppm*
"""

# ╔═╡ 87b9a664-1746-4dba-b61f-31ed626b402b
md"""
The full validation that includes all three packages would give the following plot:

![Full validation plot](https://raw.githubusercontent.com/VEZY/PlantBiophysics-paper/main/tutorials/out/figure_global_simulation.png)
"""

# ╔═╡ 0ac70149-83d3-4056-9e1b-fc91dd7d44ca
md"""
## References
"""

# ╔═╡ 7dcd8225-c3eb-48d4-9f3c-7a2db6c7a679
"""
	rgb(r, g, b, a)

Like `Colors.RGBA` but accepts colors in the 0-255 range.
"""
function rgb(r, g, b, a)
    return RGBA(r / 255, g / 255, b / 255, a)
end

# ╔═╡ c432a0c9-9254-4cf8-af04-2756f80f748e
"""
    RMSE(obs,sim)

Returns the Root Mean Squared Error between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function RMSE(obs, sim, digits=2)
    return round(sqrt(sum((obs .- sim) .^ 2) / length(obs)), digits=digits)
end

# ╔═╡ bf9a93fa-bc77-49f0-b91b-c94718af289e
"""
    nRMSE(obs,sim)

Returns the normalized Root Mean Squared Error between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function nRMSE(obs, sim; digits=2)
    return round(sqrt(sum((obs .- sim) .^ 2) / length(obs)) / (findmax(obs)[1] - findmin(obs)[1]), digits=digits)
end

# ╔═╡ d47fefa3-64c4-4ffd-8303-02f99ef98d5a
"""
    EF(obs,sim)

Returns the Efficiency Factor between observations `obs` and simulations `sim` using NSE (Nash-Sutcliffe efficiency) model.
More information can be found at https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient.
The closer to 1 the better.
"""
function EF(obs, sim, digits=2)
    SSres = sum((obs - sim) .^ 2)
    SStot = sum((obs .- mean(obs)) .^ 2)
    return round(1 - SSres / SStot, digits=digits)
end

# ╔═╡ 3d63422e-d8ba-4649-b97d-0704280a5646
"""
	    Bias(obs,sim)

	Returns the bias between observations `obs` and simulations `sim`.
	The closer to 0 the better.
	"""
function Bias(obs, sim, digits=4)
    return round(mean(sim .- obs), digits=digits)
end

# ╔═╡ a2cbcc7d-891d-4c1f-8d4c-bdd1e55c76fc
"""
	nBias(obs,sim; digits = 2)

Returns the normalised bias (%) between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function nBias(obs, sim; digits=2)
    return round(mean((sim .- obs)) / (findmax(obs)[1] - findmin(obs)[1]), digits=digits)
end

# ╔═╡ 4195aa80-f8e4-4e05-be83-ca97b0f6041e
stats =
    combine(
        groupby(df_res, [:variable, :origin], sort=true),
        [:measured, :simulated] => ((x, y) -> RMSE(x, y)) => :RMSE,
        [:measured, :simulated] => ((x, y) -> nRMSE(x, y)) => :nRMSE,
        [:measured, :simulated] => ((x, y) -> Bias(x, y)) => :Bias,
        [:measured, :simulated] => ((x, y) -> nBias(x, y)) => :nBias,
        [:measured, :simulated] => ((x, y) -> EF(x, y)) => :EF
	)

# ╔═╡ b040102a-bf17-47d7-acfd-3a80b843fc4b
begin
    transparency_col = 0.8
    transparency_fill = 0.3
    color_pb = rgb(0, 0, 0, transparency_col)
    color_lg = rgb(223, 120, 97, transparency_col)
    color_pe = rgb(118, 84, 154, transparency_col)
    fill_pb = rgb(0, 0, 0, transparency_fill)
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
    axa = Axis(fig[1, 2], title="a) Net CO₂ assimilation (Aₙ)", aspect=1, titlealign=:left)
    xlims!(-10.0, 50.0)
    ylims!(-10.0, 50.0)

    ablines!(axa, 0, 1, color=(:grey, 0.4), linewidth=4)

    # LG = scatter!(
    #     axa, df.A, df.AsimLG,
    #     color=fill_lg,
    #     markersize=ms,
    #     strokecolor=color_lg,
    #     strokewidth=stw
    # )
    # PE = scatter!(axa, df.A, df.AsimPE,
    #     color=fill_pe,
    #     markersize=ms,
    #     strokecolor=color_pe,
    #     strokewidth=stw
    # )
    PB = scatter!(axa, df.A, df.AsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )
    axislegend(
        axa,
        [PB],
        [
            "nRMSE: " * string(filter(x -> x.variable == "A" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
        ],
        "",
        position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )

    # Transpiration
    axb = Axis(fig[1, 3], title="b) Transpiration rate (E)", aspect=1, titlealign=:left)
    xlims!(-0.5, 10.0)
    ylims!(-0.5, 10.0)

    ablines!(axb, 0, 1, color=(:grey, 0.4), linewidth=4)

    # LG = scatter!(axb, df.Trmmol, df.EsimLG,
    #     color=fill_lg,
    #     markersize=ms,
    #     strokecolor=color_lg,
    #     strokewidth=stw
    # )
    # PE = scatter!(axb, df.Trmmol, df.EsimPE,
    #     color=fill_pe,
    #     markersize=ms,
    #     strokecolor=color_pe,
    #     strokewidth=stw
    # )
    PB = scatter!(axb, df.Trmmol, df.EsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )
    axislegend(
        axb,
        [PB],
        [
            "nRMSE: " * string(filter(x -> x.variable == "E" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -> x.variable == "E" && x.origin == "plantecophys", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -> x.variable == "E" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
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
        title="(c) CO₂ stomatal conductance (Gₛ)",
        aspect=1, 
		titlealign=:left
    )
    xlims!(-0.05, 0.85)
    ylims!(-0.05, 0.85)

    ablines!(axc, 0, 1, color=(:grey, 0.4), linewidth=4)

    # LG = scatter!(axc, df.Gₛ, df.GssimLG,
    #     color=fill_lg,
    #     markersize=ms,
    #     strokecolor=color_lg,
    #     strokewidth=stw,
    #     label="LeafGasExchange.jl"
    # )
    # PE = scatter!(axc, df.Gₛ, df.GssimPE,
    #     color=fill_pe,
    #     markersize=ms,
    #     strokecolor=color_pe,
    #     label="plantecophys",
    #     strokewidth=stw,
    # )
    PB = scatter!(axc, df.Gₛ, df.GssimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw,
        label="PlantBiophysics.jl"
    )

    axislegend(
        axc,
        [PB],
        [
            "nRMSE: " * string(filter(x -> x.variable == "Gs" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -> x.variable == "Gs" && x.origin == "plantecophys", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -> x.variable == "Gs" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
        ],
        "",
        position=:rb,
        orientation=:vertical,
        labelsize=legend_lab_size,
        padding=0.0,
        framevisible=false
    )


    # Leaf temperature
    axd = Axis(fig[2, 3], title="(d) Leaf temperature (Tₗ)", aspect=1, titlealign=:left)
    ablines!(axd, 0, 1, color=(:grey, 0.4), linewidth=4)

    # LG = scatter!(axd, df.Tₗ, df.TlsimLG,
    #     color=fill_lg,
    #     markersize=ms,
    #     strokecolor=color_lg,
    #     strokewidth=stw
    # )
    # PE = scatter!(axd, df.Tₗ, df.TlsimPE,
    #     color=fill_pe,
    #     markersize=ms,
    #     strokecolor=color_pe,
    #     strokewidth=stw
    # )
    PB = scatter!(axd, df.Tₗ, df.TlsimPB,
        color=fill_pb,
        markersize=ms,
        strokecolor=color_pb,
        strokewidth=stw
    )

    xlims!(10.0, 36.0)
    ylims!(10.0, 36.0)

    axislegend(
        axd, [PB],
        [
            "nRMSE: " * string(filter(x -> x.variable == "Tl" && x.origin == "PlantBiophysics.jl", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -> x.variable == "Tl" && x.origin == "plantecophys", stats).nRMSE[1]),
            # "nRMSE: " * string(filter(x -> x.variable == "Tl" && x.origin == "LeafGasExchange.jl", stats).nRMSE[1]),
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
PlantBiophysics = "7ae8fcfa-76ad-4ec6-9ea7-5f8f5e2d6ec9"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
AlgebraOfGraphics = "~0.6.11"
CSV = "~0.10.4"
CairoMakie = "~0.8.13"
Colors = "~0.12.8"
DataFrames = "~1.3.4"
PlantBiophysics = "~0.4.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0"
manifest_format = "2.0"
project_hash = "3eebe7ece32a2e2e184e2632b0ae040dc2fdd1c2"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "5c0b629df8a5566a06f5fef5100b53ea56e465a0"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.2"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "1b99dfa6ccc27abc03ac2aee0433eb2e3c2d4d4d"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.11"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "40debc9f72d0511e12d817c7ca06a721b6423ba3"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.17"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA"]
git-tree-sha1 = "387e0102f240244102814cf73fe9fbbad82b9e9e"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.8.13"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "78bee250c6826e1cf805a88b7f1e86025275d208"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.46.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "96dc5c5c8994be519ee3420953c931c55657a3f2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.24"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "992a23afdb109d0d2f8802a30cf5ae4b1fe7ea68"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "334a5896c1534bb1aa7aa2a642d30ba7707357ef"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.68"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "0fa3b52a04a4e210aeb1626def9c90df3ae65268"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.1.0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "94f5101b96d2d968ace56f7f2db19d0a5f592e28"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.15.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "316daa94fad0b7a008ebd573e002efd6609d85ac"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.19"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "3399bbad4c9e9a2fd372a54d7b67b3c7121b6402"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.3"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "5a2cff9b6b77b33b89f3d97a4d367747adce647e"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.15.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "187198a4ed8ccd7b5d99c41b69c679269ea2b2d4"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.32"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "b5c7fe9cea653443736d264b85466bad8c574f4a"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.9"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "039118892476c2bf045a43b88fcb75ed566000ff"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.8.0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "a7a97895780dab1085a97769316aa348830dc991"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.3"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "a6d30bdc378d340912f48abf01281aab68c0dec8"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.2"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "53c7e69a6ffeb26bd594f5a1421b889e7219eeaa"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "d19f9edd8c34760dca2de2b503f969d8700ed288"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "64f138f9453a018c8f3562e7bae54edc059af249"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.4"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "076bb0da51a8c8d1229936a1af7bdfacd65037e1"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.2"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "46efcea75c890e5d820e670516dc156689851722"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.4"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "91aa1442e63a77f101aff01dec5a821a17f43922"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.12.1"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "41d162ae9c868218b1f3fe78cba878aa348c2d26"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "UnicodeFun"]
git-tree-sha1 = "b0323393a7190c9bf5b03af442fc115756df8e59"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.17.13"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "fbf705d2bdea8fc93f1ae8ca2965d8e03d4ca98c"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.4.0"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "114ef48a73aea632b8aebcb84f796afcc510ac7c"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.4.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MetaGraphsNext]]
deps = ["Graphs", "JLD2"]
git-tree-sha1 = "6f3371197acf57728e3d27e0b3f113824f808c95"
uuid = "fa8bd995-216d-47f1-8a91-f3b68fbeb377"
version = "0.3.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MultiScaleTreeGraph]]
deps = ["AbstractTrees", "DataFrames", "Dates", "DelimitedFiles", "Graphs", "MetaGraphsNext", "MutableNamedTuples", "OrderedCollections", "Printf", "SHA", "XLSX"]
git-tree-sha1 = "af07d47dc843fe3c207ee065a84d71d1df84fa06"
uuid = "dd4a991b-8a45-4075-bede-262ee62d5583"
version = "0.7.0"

[[deps.MutableNamedTuples]]
git-tree-sha1 = "44964b7961bdb59c88f42ee681061dd38c384780"
uuid = "af6c499f-54b4-48cc-bbd2-094bba7533c7"
version = "0.1.2"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "dfd8d34871bc3ad08cd16026c1828e271d554db9"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.1"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "1ea784113a6aa054c5ebd95945fa5e52c2f378e7"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.7"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "e925a64b8585aa9f4e3047b8d2cdc3f0e79fd4e4"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.16"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "1155f6f937fa2b94104162f01fa400e192e4272f"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.2"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a121dfbba67c94a5bec9dde613c3d0cbcf3a12b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.3+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlantBiophysics]]
deps = ["CSV", "DataFrames", "Dates", "LsqFit", "MultiScaleTreeGraph", "MutableNamedTuples", "OrderedCollections", "RecipesBase", "Statistics", "Test", "YAML"]
git-tree-sha1 = "0588e820b66b3383144dad70bf9f8c2bd4b797b8"
uuid = "7ae8fcfa-76ad-4ec6-9ea7-5f8f5e2d6ec9"
version = "0.4.4"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "3c009334f45dfd546a16a57960a821a1a023d241"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.5.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
git-tree-sha1 = "7dbc15af7ed5f751a82bf3ed37757adf76c32402"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.1"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "dfec37b90740e3b9aa5dc2613892a3fc155c3b42"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.6"

[[deps.StaticArraysCore]]
git-tree-sha1 = "ec2bd695e905a3c755b33026954b119ea17f2d22"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.3.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "f8ba54b202c77622a713e25e7616d618308b34d3"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.31"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "50ccd5ddb00d19392577902f0079267a72c5ab04"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.5"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "8c6ac65ec9ab781af05b08ff305ddc727c25f680"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.12"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "70e6d2da9210371c927176cb7a56d41ef1260db7"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XLSX]]
deps = ["Artifacts", "Dates", "EzXML", "Printf", "Tables", "ZipFile"]
git-tree-sha1 = "187f7c392881b1dd3453ac0cab36c3f70d1bb239"
uuid = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"
version = "0.8.3"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "3c6e8b9f5cdaaa21340f841653942e1a6b6561e5"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.7"

[[deps.ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "ef4f23ffde3ee95114b461dc667ea4e6906874b2"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.10.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╟─1100f4ba-d45e-11ec-25f5-b79c9d3b9851
# ╟─1609ce3c-4066-45ca-acce-8dbf081bfb89
# ╟─74ae23d5-71bc-4ad8-acdf-0a38a1b129bd
# ╟─aff7314d-7a8a-4a47-8a77-252e96377e9c
# ╟─c2e2c3c8-288c-491a-a8c7-8adf4517103c
# ╟─ab97a0a5-d155-4758-a1d4-70e3a0b8c796
# ╠═3c5bcd7e-87fb-4aef-86d8-326921b4423c
# ╟─e8961d19-85d4-484b-a029-359a9dbea800
# ╠═f9be64a3-b788-48d8-a925-5954e532effb
# ╟─3592d0a1-d658-4137-8ea8-e8ea9ff3dc91
# ╠═a65d0175-f6cf-4bc9-8a2f-fbde4d2692b6
# ╠═44db3882-e5b3-4427-b781-e5185e2d70af
# ╟─92deadb1-537b-4c3d-b0d2-3fd3cba171e4
# ╟─6e3a7997-17dd-453e-879e-8f8d6b46ca2d
# ╟─8fe5544c-8f61-4306-b424-01d62e25fbcc
# ╟─96af7e3e-38ce-4a8a-b503-163e88ddfb8f
# ╟─232de02a-9440-4df5-8fbc-1e2eb6dae24a
# ╟─56705c8b-58ac-4c2d-acc5-dcbb03f25bf0
# ╠═db357873-44b5-4110-a80b-98d3adf785e6
# ╟─5ac831c1-cc9e-47e9-a11b-bd191a1d55a2
# ╟─7a09574b-b06e-4c36-9fea-7757c20fd763
# ╟─3d4375ea-2020-4829-ba70-5406c921c943
# ╠═53a0fa15-d042-4e00-84ea-f9caa5fc3c3e
# ╟─9f4d46ab-cc75-47e1-aef1-b7ada0e01e29
# ╠═182787c9-5075-451f-bdc7-5a8add68d141
# ╟─2fabe562-5eb9-4ef1-897e-f6e380d1518c
# ╠═4195aa80-f8e4-4e05-be83-ca97b0f6041e
# ╟─e8171ee6-36fa-4227-8f78-2b0bd8675cdf
# ╟─f180ca3d-490c-4a64-8608-4f2f5e57c41f
# ╟─b040102a-bf17-47d7-acfd-3a80b843fc4b
# ╟─dd192660-2c6d-46a2-a06a-f0c1d01a9598
# ╟─87b9a664-1746-4dba-b61f-31ed626b402b
# ╟─0ac70149-83d3-4056-9e1b-fc91dd7d44ca
# ╟─7dcd8225-c3eb-48d4-9f3c-7a2db6c7a679
# ╟─c432a0c9-9254-4cf8-af04-2756f80f748e
# ╟─bf9a93fa-bc77-49f0-b91b-c94718af289e
# ╟─d47fefa3-64c4-4ffd-8303-02f99ef98d5a
# ╟─3d63422e-d8ba-4649-b97d-0704280a5646
# ╟─a2cbcc7d-891d-4c1f-8d4c-bdd1e55c76fc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
