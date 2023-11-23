module InputPrep
export load_from_examplename, build_aux_structures, create_triangulation, plot_triangulation

include("UserFunctions.jl")

using DelimitedFiles
using SparseArrays
using LinearAlgebra
using Triangulate
using PyCall
using PyPlot
matplotlib = pyimport("matplotlib")
matplotlib.use("TkAgg")

using .UserFunctions

function read_file(filepath, filename, is_mandatory=true, vartype=Int)

    if isfile(filepath * filename)
        return readdlm(filepath * filename, ' ', vartype, '\n')
    else
        if is_mandatory
            error("File $filename not found in $filepath")
        else
            return []
        end
    end
end

function load_from_examplename(example::String)

    if example == "9.1"
        filepath = "data/section9.1/"
        f = fcn_zeros
        u_D = fcn_zeros
        g = g_91
    elseif example == "3.1"
        filepath = "data/section3.1/"
        f = fcn_ones
        u_D = fcn_zeros
        g = fcn_zeros
    end

    coordinate = read_file(filepath, "coordinate.dat", true, Float64)
    element = read_file(filepath, "element.dat")
    dirichlet = read_file(filepath, "dirichlet.dat", false)
    neumann = read_file(filepath, "neumann.dat", false)

    return coordinate, element, dirichlet, neumann, f, u_D, g, filepath

end


function create_triangulation(points)
    triin = TriangulateIO()
    triin.pointlist = points
    (triout, _) = Triangulate.triangulate("cQ", triin)

    coordinate = (triout.pointlist)'
    element = (triout.trianglelist)'

    return coordinate, element
end

function plot_triangulation(coordinate, element)
    fig = figure()
    ax = fig.add_subplot(111)
    for j in axes(element, 1)
        x = coordinate[element[j, :], 1]
        y = coordinate[element[j, :], 2]
        ax.plot(x, y, ".k", markersize=10)
        ax.plot([x[1], x[2]], [y[1], y[2]], "-b")
        ax.plot([x[2], x[3]], [y[2], y[3]], "-b")
        ax.plot([x[3], x[1]], [y[3], y[1]], "-b")
    end
    ax.set_title("Triangulation")
    ax.set_aspect("equal")

    show()
end


end