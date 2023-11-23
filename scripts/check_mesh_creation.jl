import MfemPoisson

# Create a grid of points in a square
side_length = 1
num_points = 5
x = LinRange(0, side_length, num_points)
y = LinRange(0, side_length, num_points)
points = hcat([[i, j] for i in x for j in y]...)

coordinate, element, _ = MfemPoisson.create_triangulation(points)
MfemPoisson.plot_triangulation(coordinate, element)