module PlotFields
export displacement_field, flux_fields

using DelimitedFiles
using PyCall
using PyPlot
matplotlib = pyimport("matplotlib")
matplotlib.use("Agg")

function displacement_field(coordinate, element, u, filepath=nothing)
    fig = figure()
    ax = fig.add_subplot(111, projection="3d")
    for j in axes(element, 1)
        ax.plot_trisurf(coordinate[element[j, :], 1], coordinate[element[j, :], 2], ones(3) .* u[j], linewidth=1.5)
    end
    ax.set_title("Displacement field")
    if !isnothing(filepath)
        fig.savefig(filepath * "displacement_field.png")
    end
end

function flux_fields(coordinate, element, p, filepath=nothing)
    fig = figure()
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")
    for j in axes(element, 1)

        ax1.plot_trisurf(coordinate[element[j, :], 1], coordinate[element[j, :], 2], p[3*(j-1).+[1, 2, 3], 1])
        ax2.plot_trisurf(coordinate[element[j, :], 1], coordinate[element[j, :], 2], p[3*(j-1).+[1, 2, 3], 2])
    end
    ax1.set_title("Flux field in x-direction")
    ax2.set_title("Flux field in y-direction")

    if !isnothing(filepath)
        fig.savefig(filepath * "flux_fields.png")
    end
end

end