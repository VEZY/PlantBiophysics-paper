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


struct ForcedGs <: PlantBiophysics.AbstractGsModel
    g0::Float64
    gs_min::Float64
end
ForcedGs() = ForcedGs(0.0, 0.0)

function PlantBiophysics.gs_closure(::ForcedGs, models, status, meteo=missing)
    # (status.Gₛ - models.stomatal_conductance.g0) / status.A
    status.Gₛ
end

function PlantBiophysics.stomatal_conductance!_(::ForcedGs, models, status, gs_closure)
    status.Gₛ
end

function PlantBiophysics.stomatal_conductance!_(::ForcedGs, models, status, meteo::M, constants=Constants()) where {M<:PlantBiophysics.AbstractAtmosphere}
    status.Gₛ
end

function PlantBiophysics.inputs_(::ForcedGs)
    NamedTuple()
end

function PlantBiophysics.outputs_(::ForcedGs)
    (Gₛ=-999.99,)
end
