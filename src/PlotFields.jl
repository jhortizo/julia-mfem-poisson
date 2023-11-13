module PlotFields
export displacement_field, flux_fields

using DelimitedFiles
using PyCall
using PyPlot
using Statistics
matplotlib = pyimport("matplotlib")
matplotlib.use("Agg")

function estimate_values_colors(z, colormap="viridis")
    colormap = get_cmap(colormap)
    normed_values = (z .- minimum(z)) ./ (maximum(z) - minimum(z))
    colors = colormap(normed_values)
    return colors
end

function displacement_field(coordinate, element, u, filepath=nothing)
    
    colors = estimate_values_colors(u)
    fig = figure()
    ax = fig.add_subplot(111, projection="3d")
    for j in axes(element, 1)
        x = coordinate[element[j, :], 1]
        y = coordinate[element[j, :], 2]
        z = ones(3) .* u[j]
        color = colors[j, :]
        println(color)
        ax.plot_trisurf(x, y, z, linewidth=1.5, color=color)
    end
    ax.set_title("Displacement field")
    ax.view_init(elev=30, azim=200)
    if !isnothing(filepath)
        fig.savefig(filepath * "displacement_field.png")
    end
end

function flux_fields(coordinate, element, p, filepath=nothing)
    # estimate mean flux for each element
    mean_xs = zeros(size(element, 1))
    mean_ys = zeros(size(element, 1))
    for j in axes(element, 1)
        mean_xs[j] = mean(p[3*(j-1).+[1, 2, 3], 1])
        mean_ys[j] = mean(p[3*(j-1).+[1, 2, 3], 2])
    end

    colors_x = estimate_values_colors(mean_xs, "plasma")
    colors_y = estimate_values_colors(mean_ys, "plasma")

    fig = figure()
    ax1 = fig.add_subplot(121, projection="3d")
    ax2 = fig.add_subplot(122, projection="3d")
    for j in axes(element, 1)
        x = coordinate[element[j, :], 1]
        y = coordinate[element[j, :], 2]
        flux_x = p[3*(j-1).+[1, 2, 3], 1]
        flux_y = p[3*(j-1).+[1, 2, 3], 2]

        color_x = colors_x[j, :]
        color_y = colors_y[j, :]
        
        ax1.plot_trisurf(x, y, flux_x, linewidth=1.5, color=color_x)
        ax2.plot_trisurf(x, y, flux_y, linewidth=1.5, color=color_y)
    end
    ax1.set_title("Flux field in x-direction")
    ax1.view_init(elev=30, azim=200)

    ax2.set_title("Flux field in y-direction")
    ax2.view_init(elev=30, azim=200)

    if !isnothing(filepath)
        fig.savefig(filepath * "flux_fields.png")
    end
end

end