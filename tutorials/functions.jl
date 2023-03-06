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
        strokewidth=3,
    )

    lines!(
        ax,
        pmean.(var_sim),
        color=line_color,
        linewidth=2.5,
        label="Simulation ± 95% confidence interval",
    )

    errorbars!(
        ax,
        eachindex(var_sim),
        pmean.(var_sim),
        pstd.(var_sim),
        pstd.(var_sim),
        color=error_color,
        whiskerwidth=5,
        linewidth=2,
    )
end


struct ForcedGs <: PlantBiophysics.AbstractStomatal_ConductanceModel
    g0::Float64
    gs_min::Float64
end
ForcedGs() = ForcedGs(0.0, 0.0)

function PlantBiophysics.gs_closure(::ForcedGs, models, status, meteo=missing)
    # (status.Gₛ - models.stomatal_conductance.g0) / status.A
    status.Gₛ
end

function PlantSimEngine.run!(::ForcedGs, models, status, gs_closure)
    status.Gₛ
end

function PlantSimEngine.run!(
    ::ForcedGs,
    models,
    status,
    meteo::M,
    constants=Constants(),
) where {M<:PlantMeteo.AbstractAtmosphere}
    status.Gₛ
end

function PlantSimEngine.inputs_(::ForcedGs)
    NamedTuple()
end

function PlantSimEngine.outputs_(::ForcedGs)
    (Gₛ=-Inf,)
end



"""
    RMSE(obs,sim)

Returns the Root Mean Squared Error between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function RMSE(obs, sim, digits=2)
    return round(sqrt(sum((obs .- sim) .^ 2) / length(obs)), digits=digits)
end

"""
    nRMSE(obs,sim)

Returns the normalized Root Mean Squared Error between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function nRMSE(obs, sim; digits=2)
    return round(
        sqrt(sum((obs .- sim) .^ 2) / length(obs)) / (findmax(obs)[1] - findmin(obs)[1]),
        digits=digits,
    )
end

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

"""
	    Bias(obs,sim)

	Returns the bias between observations `obs` and simulations `sim`.
	The closer to 0 the better.
	"""
function Bias(obs, sim, digits=4)
    return round(mean(sim .- obs), digits=digits)
end

"""
	nBias(obs,sim; digits = 2)

Returns the normalised bias (%) between observations `obs` and simulations `sim`.
The closer to 0 the better.
"""
function nBias(obs, sim; digits=2)
    return round(mean((sim .- obs)) / (findmax(obs)[1] - findmin(obs)[1]), digits=digits)
end
