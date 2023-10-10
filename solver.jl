include("UserFunctions.jl")

using DelimitedFiles
using LinearAlgebra
using SparseArrays
using .UserFunctions

filepath = "data/section3.1/"

f = fcn_ones
u_D = fcn_zeros


function stimaB(coord)
    N = coord[:] * ones(1, 3) - repeat(coord, 3, 1)
    C = diagm([norm(N[[5, 6], 2]), norm(N[[1, 2], 3]), norm(N[[1, 2], 2])])
    M = spdiagm(-4 => ones(2), -2 => ones(4), 0 => 2 * ones(6), 2 => ones(4), 4 => ones(2))
    C * N' * M * N * C / (24 * det([[1 1 1]; coord]))
end

function read_file(filepath, filename, is_mandatory=True)

    if isfile(filepath * filename)
        return readdlm(filepath * filename, ' ', Int, '\n')
    else
        if is_mandatory
            error("File $filename not found in $filepath")
        else
            return []
        end
    end
end

coordinate = read_file(filepath, "coordinate.dat")
element = read_file(filepath, "element.dat")
dirichlet = read_file(filepath, "dirichlet.dat", false)
neumann = read_file(filepath, "neumann.dat", false)

nodes2element = spzeros(size(coordinate, 1), size(coordinate, 1))
for j = axes(element, 1)
    nodes2element[element[j, :], element[j, [2, 3, 1]]] .+= j .* Matrix(1I, 3, 3)
end

B = nodes2element + nodes2element';
indices = findall(!iszero, triu(B))
rows, cols = [i[1] for i in indices], [i[2] for i in indices]
nodes2edge = sparse(rows, cols, 1:size(rows, 1), size(coordinate, 1), size(coordinate, 1))
nodes2edge = nodes2edge + nodes2edge'

noedges = size(rows, 1)
edge2element = zeros(noedges, 4)
for m = axes(element, 1)
    for k = 1:3
        initial_node = element[m, k]
        end_node = element[m, k%3+1]
        p = nodes2edge[initial_node, end_node]
        if edge2element[p, 1] == 0
            edge2element[p, :] = [initial_node end_node nodes2element[initial_node, end_node] nodes2element[end_node, initial_node]]
        end
    end
end


B = spzeros(noedges, noedges)
C = spzeros(noedges, size(element, 1))
for j = axes(element, 1)
    local coord = coordinate[element[j, :], :]'
    dummy = nodes2edge[element[j, [2 3 1]], element[j, [3 1 2]]]
    dummy = dropdims(dummy, dims=Dims(findall(size(dummy) .== 1)))
    rows = diag(dummy)
    signum = ones(1, 3)
    signum[findall(j .== edge2element[rows, 4])] .= -1
    n = coord[:, [3 1 2]][:, 1, :] - coord[:, [2 3 1]][:, 1, :]
    B[rows, rows] .+= diagm(signum[:]) * stimaB(coord) * diagm(signum[:])
    C[rows, j] = diagm(signum[:]) * [norm(n[:, 1]) norm(n[:, 2]) norm(n[:, 3])]'
end

A = spzeros(noedges + size(element, 1), noedges + size(element, 1))
A = [B C; C' spzeros(size(element, 1), size(element, 1))]


b = zeros(noedges + size(element, 1), 1)
for l = axes(element, 1)
    coord = coordinate[element[l, :], :]
    b[noedges+l] = -det([[1 1 1]; coord']) * f(sum(coord) / 3)[1] / 6
end

for k = axes(dirichlet, 1)
    this_diri = dirichlet[k, :]
    this_edge = nodes2edge[this_diri[1], this_diri[2]]
    b[this_edge] = norm(coordinate[this_diri[1], :] - coordinate[this_diri[2], :]) * u_D(sum(coordinate[this_diri, :]) / 2)[1]
end

# Neumann Condition


x = A \ b

writedlm(filepath * "solution.dat", x, ' ')