########################################################################
# Taylor diagram for global simulation of ACis curves from Medlyn data
# Simon Treillou, 2023
########################################################################

using TaylorDiag
using Plots
using CSV
using DataFrames
using Colors
using Measures
plotly()

# Loading previous results
df_res = DataFrame(CSV.File("out/df_res.csv"))

########################################################################
# Defining accurate colors for plotting
########################################################################

function rgb(r, g, b, a)
    return RGBA(r / 255, g / 255, b / 255, a)
end

col_pb = (244, 124, 124)    # PlantBiophysics.jl
col_lg = (93, 174, 139)     # LeafGasExchange.jl
col_pe = (112, 161, 215)    # plantecophys
transparency_col = 0.6      # Transparency for edges
transparency_fill = 0.4     # Transparency for fill

color_pb = rgb(col_pb..., transparency_col)
color_pe = rgb(col_lg..., transparency_col)
color_lg = rgb(col_pe..., transparency_col)
fill_pb = rgb(col_pb..., transparency_fill)
fill_pe = rgb(col_lg..., transparency_fill)
fill_lg = rgb(col_pe..., transparency_fill)

cols = [:black,fill_pb,fill_lg,fill_pe] 
strkcols = [:black,color_pb,color_lg,color_pe] 

# Parameters for fitting
msize = 7                   # Marker size
stw = 1.5                   
legend_lab_size=10
xleg = 1.6
yleg = 0.4

########################################################################
# Computing and plotting standard deviations and correlations for carbon assimilation
########################################################################

var = "A"
obs   = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"PlantBiophysics.jl")==0),df_res).measured
modPB = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"PlantBiophysics.jl")==0),df_res).simulated
modPE = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"plantecophys")==0),df_res).simulated
modLG = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"LeafGasExchange.jl")==0),df_res).simulated

S = [STD(obs),STD(modPB), STD(modPE), STD(modLG)]
S = [S[i]/S[1] for i in eachindex(S)]
C = [COR(obs,obs),COR(obs,modPB),COR(obs,modPE),COR(obs,modLG)]

# Plotting initial Taylor diagram
nms = ["","","",""]
fig =taylordiagram([S[1]],[C[1]],[nms[1]],normalize=true,ang=pi/2)
nms = ["Data","PlantBiophysics.jl","plantecophys","LeafGasExchange.jl"]

# Defining polar coordinates
rho   = S
theta = to_polar(C)
limSTD   = findmax(S)[1]*2

    # Plotting reference and model points
    for i in 1:length(theta)
        Plots.scatter!([cos.(theta[i]).*rho[i]], [sin.(theta[i]).*rho[i]],markerstrokecolor=strkcols[i],markercolor=cols[i],markershape=:circle,markersize=msize,markerstrokewidth=stw,label=nms[i])
    end
plot!(fontfamily="NotoSans-Regular.ttf",fontsize=12)
plot!(size = (700,700))
scatter!([xleg],[yleg], markerstrokecolor=:grey,markercolor=:grey,markershape=:circle,markersize=msize,markerstrokewidth=stw, markerstrokealpha=0.6,markeralpha=0.4,label="")
annotate!(xleg+0.05, yleg-0.005, text(string("A"), :left, legend_lab_size))


########################################################################
# Computing and plotting standard deviations and correlations for transpiration
########################################################################

var = "E"
obs   = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"PlantBiophysics.jl")==0),df_res).measured
modPB = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"PlantBiophysics.jl")==0),df_res).simulated
modPE = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"plantecophys")==0),df_res).simulated
modLG = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"LeafGasExchange.jl")==0),df_res).simulated

S = [STD(obs),STD(modPB), STD(modPE), STD(modLG)]
S = [S[i]/S[1] for i in eachindex(S)]
C = [COR(obs,obs),COR(obs,modPB),COR(obs,modPE),COR(obs,modLG)]


# Defining polar coordinates
rho   = S
theta = to_polar(C)
limSTD   = findmax(S)[1]*2

    # Plotting reference and model points
    for i in 2:length(theta)
        Plots.scatter!([cos.(theta[i]).*rho[i]], [sin.(theta[i]).*rho[i]],markerstrokecolor=strkcols[i],color=cols[i],markershape=:rect,markersize=msize,markerstrokewidth=stw,label="")
    end
Plots.scatter!()
scatter!([xleg],[yleg-0.1], markerstrokecolor=:grey,markercolor=:grey,markershape=:rect,markersize=msize,markerstrokewidth=stw, markerstrokealpha=0.6,markeralpha=0.4,label="")
annotate!(xleg+0.05, yleg-0.1-0.005, text(string("E"), :left, legend_lab_size))

########################################################################
# Computing and plotting standard deviations and correlations for leaf temperature
########################################################################

var = "Tl"
obs   = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"PlantBiophysics.jl")==0),df_res).measured
modPB = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"PlantBiophysics.jl")==0),df_res).simulated
modPE = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"plantecophys")==0),df_res).simulated
modLG = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"LeafGasExchange.jl")==0),df_res).simulated

S = [STD(obs),STD(modPB), STD(modPE), STD(modLG)]
S = [S[i]/S[1] for i in eachindex(S)]
C = [COR(obs,obs),COR(obs,modPB),COR(obs,modPE),COR(obs,modLG)]

# Defining polar coordinates
rho   = S
theta = to_polar(C)
limSTD   = findmax(S)[1]*2

    # Plotting reference and model points
    for i in 2:length(theta)
        Plots.scatter!([cos.(theta[i]).*rho[i]], [sin.(theta[i]).*rho[i]],markerstrokecolor=strkcols[i],color=cols[i],markershape=:diamond,markersize=msize,markerstrokewidth=stw,label="")
    end
Plots.scatter!()
scatter!([xleg],[yleg-0.2], markerstrokecolor=:grey,markercolor=:grey,markershape=:diamond,markersize=msize,markerstrokewidth=stw, markerstrokealpha=0.6,markeralpha=0.4,label="")
annotate!(xleg+0.05, yleg-0.2-0.005, text(string("Tₗ"), :left, legend_lab_size))

########################################################################
# Computing and plotting standard deviations and correlations for stomatal conductance
########################################################################

var = "Gs"
obs   = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"PlantBiophysics.jl")==0),df_res).measured
modPB = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"PlantBiophysics.jl")==0),df_res).simulated
modPE = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"plantecophys")==0),df_res).simulated
modLG = filter(x->(cmp.(x.variable,var)==0)&(cmp.(x.origin,"LeafGasExchange.jl")==0),df_res).simulated

S = [STD(obs),STD(modPB), STD(modPE), STD(modLG)]
S = [S[i]/S[1] for i in eachindex(S)]
C = [COR(obs,obs),COR(obs,modPB),COR(obs,modPE),COR(obs,modLG)]

# Defining polar coordinates
rho   = S
theta = to_polar(C)
limSTD   = findmax(S)[1]*2

    # Plotting reference and model points
    for i in 2:length(theta)
        Plots.scatter!([cos.(theta[i]).*rho[i]], [sin.(theta[i]).*rho[i]],markerstrokecolor=strkcols[i],color=cols[i],markershape=:utriangle,markersize=msize,markerstrokewidth=stw,label="")
    end
Plots.scatter!()
scatter!([xleg],[yleg-0.3], markerstrokecolor=:grey,markercolor=:grey,markershape=:utriangle,markersize=msize,markerstrokewidth=stw, markerstrokealpha=0.6,markeralpha=0.4,label="")
annotate!(xleg+0.05, yleg-0.3-0.005, text(string("Gₛ"), :left, legend_lab_size))

plot!(legend=:topright, foreground_color_legend = nothing)
