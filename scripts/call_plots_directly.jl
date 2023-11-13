using MfemPoisson
using DelimitedFiles


example="3.1"

coordinate, element, dirichlet, neumann, f, u_D, g, filepath = load_from_examplename(example)
_, nodes2edge, edge2element, noedges = build_aux_structures(coordinate, element)

x = readdlm(filepath * "solution.dat", ' ', Float64, '\n')
u, p = calculate_fields(coordinate, element, nodes2edge, edge2element, x)

plot_and_save(coordinate, element, x, u, p, )