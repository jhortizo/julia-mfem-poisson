using DelimitedFiles
using PyPlot
pygui(true)

filepath = "data/section3.1/"

coordinate = readdlm(filepath * "coordinate.dat", ' ', Float64, '\n');
element = readdlm(filepath * "element.dat", ' ', Int, '\n');

x = readdlm(filepath * "solution.dat", ' ', Float64, '\n');
u = x[end-size(element,1) + 1:end, 1];
for j=axes(element, 1)
    plot_trisurf(coordinate[element[j, :], 1], coordinate[element[j, :], 2], ones(3) .* u[j], linewidth=1.5)
end
show()  