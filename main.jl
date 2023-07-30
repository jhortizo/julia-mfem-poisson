using DelimitedFiles
using LinearAlgebra
using SparseArrays

filepath = "data/"

coords = readdlm(filepath * "coordinates.dat", ' ', Float64, '\n');
elements = readdlm(filepath * "element.dat", ' ', Int, '\n');
dirichlet = readdlm(filepath * "dirichlet.dat", ' ', Int, '\n');

nodes2element = spzeros(size(coords, 1), size(coords, 1))
for j = axes(elements, 1)
    nodes2element[elements[j, :], elements[j, [2, 3, 1]]] .= j .* Matrix(1I, 3, 3)
end

