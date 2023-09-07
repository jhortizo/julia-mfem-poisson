using DelimitedFiles
using LinearAlgebra
using SparseArrays

filepath = "data/"

coordinate = readdlm(filepath * "coordinate.dat", ' ', Float64, '\n');
element = readdlm(filepath * "element.dat", ' ', Int, '\n');
dirichlet = readdlm(filepath * "dirichlet.dat", ' ', Int, '\n');

nodes2element = spzeros(size(coords, 1), size(coords, 1))
for j = axes(element, 1)
    nodes2element[element[j, :], element[j, [2, 3, 1]]] .+= j .* Matrix(1I, 3, 3)
end

B = nodes2element + nodes2element';
indices = findall(!iszero, triu(B))
rows, cols = [i[1] for i in indices], [i[2] for i in indices]
nodes2edge = sparse(rows, cols, 1:size(rows, 1), size(coords, 1), size(coords, 1))
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

function stimaB(coord)
    N = coord[:] * ones(1, 3) - repmat(coord, 3, 1)
    C = diag([norm(N[5,6], 2), norm(N[1, 2], 3), norm(N[1, 2], 2)])
    M = spdiags([ones(6,1), ones(6,1), 2 * ones(6,1), ones(6,1)], [-4,-2,0, 2, 4], 6, 6)
    C * N' * M * N * C / (24 * det([1, 1, 1; coord]))
end 

B = sparse(noedges, noedges)
C = sparse(noedges, size(element, 1))
for j=axes(element, 1)
    coord = coordinate(element[j, :], :)'
    rows = diag(nodes2edge[element[j, [2 3 1]], element[j, [3 1 2]]])
    signum = ones(1, 3)
    signum[findall(j .== edge2element[rows, 4])] = -1
    n = coord[:, [3 1 2]] - coord[:, [2 3 1]]
    B[rows, rows] .+= diag(signum) * stimaB(coord) * diag(signum)
    C[rows, j] = diag(signum) * [norm(n[:, 1]) norm(n[:, 2]) norm(n[:, 3])]'
end

A = sparse(noedges + size(element, 1), noedges + size(element, 1))
A = [B C; C' sparse(size(element, 1), size(element, 1))]
