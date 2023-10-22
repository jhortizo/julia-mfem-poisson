module PlotFields
export displacement_field

using DelimitedFiles
using PyPlot
pygui(true)

function displacement_field(coordinate, element, u)
    for j = axes(element, 1)
        plot_trisurf(coordinate[element[j, :], 1], coordinate[element[j, :], 2], ones(3) .* u[j], linewidth=1.5)
    end
    show()
end

end