using Triangulate
using PyCall
using PyPlot
matplotlib = pyimport("matplotlib")
matplotlib.use("TkAgg")

# Create a grid of points in a square
side_length = 1
num_points = 5
x = LinRange(0, side_length, num_points)
y = LinRange(0, side_length, num_points)
points = hcat([[i, j] for i in x for j in y]...)

# Create Delaunay triangulation
triin = TriangulateIO()
triin.pointlist = points
(triout, vorout) = Triangulate.triangulate("cV", triin)

coordinate = (triout.pointlist)'
element = (triout.trianglelist)'

# Plot the triangulation
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
