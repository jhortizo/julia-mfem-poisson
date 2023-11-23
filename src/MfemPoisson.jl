module MfemPoisson
export run_exact_paper_example, load_from_examplename, build_aux_structures, plot_and_save, calculate_fields


include("PlotFields.jl")
include("InputPrep.jl")
include("MfemSolver.jl")


using LinearAlgebra
using SparseArrays
using .PlotFields
using .InputPrep
using .MfemSolver


function run_exact_paper_example(example::String)
    # TODO: Handle case where no border condition is defined

    coordinate, element, dirichlet, neumann, f, u_D, g, filepath = load_from_examplename(example)

    x, u, p = EBMfemSolver(coordinate, element, dirichlet, neumann, f, u_D, g)

    plot_and_save(coordinate, element, x, u, p, filepath)
end

end