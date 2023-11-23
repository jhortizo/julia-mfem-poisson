module MfemSolver

export EBMfemSolver, build_aux_structures

using LinearAlgebra
using SparseArrays

function build_aux_structures(coordinate, element)
    # calculate matrices nodes2element, nodes2edge and edge2element and number of edges (noedges)
    nodes2element = spzeros(size(coordinate, 1), size(coordinate, 1))
    for j in axes(element, 1)
        nodes2element[element[j, :], element[j, [2, 3, 1]]] .+= j .* Matrix(1I, 3, 3)
    end

    B = nodes2element + nodes2element'
    indices = findall(!iszero, triu(B))
    rows, cols = [i[1] for i in indices], [i[2] for i in indices]
    nodes2edge = sparse(rows, cols, 1:size(rows, 1), size(coordinate, 1), size(coordinate, 1))
    nodes2edge = nodes2edge + nodes2edge'

    noedges = size(rows, 1)
    edge2element = zeros(noedges, 4)
    for m in axes(element, 1)
        for k in 1:3
            initial_node = element[m, k]
            end_node = element[m, k%3+1]
            p = nodes2edge[initial_node, end_node]
            if edge2element[p, 1] == 0
                edge2element[p, :] = [initial_node end_node nodes2element[initial_node, end_node] nodes2element[end_node, initial_node]]
            end
        end
    end

    return nodes2edge, edge2element, noedges
end

function calculate_local_stiffness_matrix(coord)
    N = coord[:] * ones(1, 3) - repeat(coord, 3, 1)
    C = diagm([norm(N[[5, 6], 2]), norm(N[[1, 2], 3]), norm(N[[1, 2], 2])])
    M = spdiagm(-4 => ones(2), -2 => ones(4), 0 => 2 * ones(6), 2 => ones(4), 4 => ones(2))
    C * N' * M * N * C / (24 * det([[1 1 1]; coord]))
end

function assemble_global_stiffness_matrix(coordinate, element, nodes2edge, edge2element, noedges)
    B = spzeros(noedges, noedges)
    C = spzeros(noedges, size(element, 1))
    for j in axes(element, 1)
        local coord = coordinate[element[j, :], :]'
        dummy = nodes2edge[element[j, [2 3 1]], element[j, [3 1 2]]]
        dummy = dropdims(dummy, dims=Dims(findall(size(dummy) .== 1)))
        local rows = diag(dummy)
        signum = ones(1, 3)
        signum[findall(j .== edge2element[rows, 4])] .= -1
        n = coord[:, [3 1 2]][:, 1, :] - coord[:, [2 3 1]][:, 1, :]
        B_local = calculate_local_stiffness_matrix(coord)
        B[rows, rows] .+= diagm(signum[:]) * B_local * diagm(signum[:])
        C[rows, j] = diagm(signum[:]) * [norm(n[:, 1]) norm(n[:, 2]) norm(n[:, 3])]'
    end

    A = spzeros(noedges + size(element, 1), noedges + size(element, 1))
    A = [B C; C' spzeros(size(element, 1), size(element, 1))]

    return A

end

function calculate_constants_vector(coordinate, element, nodes2edge, noedges, f, u_D, dirichlet)

    # Add volume forces
    b = zeros(noedges + size(element, 1), 1)
    for l in axes(element, 1)
        coord = coordinate[element[l, :], :]
        b[noedges+l] = -det([[1 1 1]; coord']) * f(sum(coord) / 3)[1] / 6
    end

    # Add Dirichlet condition
    # TODO: Handle case where no Dirichlet condition
    for k in axes(dirichlet, 1)
        this_diri = dirichlet[k, :]
        this_edge = nodes2edge[this_diri[1], this_diri[2]]
        b[this_edge] = norm(coordinate[this_diri[1], :] - coordinate[this_diri[2], :]) * u_D(sum(coordinate[this_diri, :]) / 2)[1]
    end

    return b
end

function solve_system(A, b, neumann, coordinate, element, g, noedges, nodes2edge)

    if !isempty(neumann)
        # Add Neumann conditions
        tmp = zeros(noedges + size(element, 1), 1)
        tmp[diag(nodes2edge[neumann[:, 1], neumann[:, 2]])] = ones(size(diag(nodes2edge[neumann[:, 1], neumann[:, 2]]), 1), 1)
        FreeEdge = findall(iszero, tmp)
        FreeRows = [i[1] for i in FreeEdge]
        x = zeros(noedges + size(element, 1), 1)
        CN = coordinate[neumann[:, 2], :] - coordinate[neumann[:, 1], :]
        for j in axes(neumann, 1)
            x[nodes2edge[neumann[j, 1], neumann[j, 2]]] = g(sum(coordinate[neumann[j, :], :], dims=1) / 2, CN[j, :]' * [0 -1; 1 0] / norm(CN[j, :]))
        end
        b = b - A * x
        # Solve for unrestricted edges
        x[FreeRows] = A[FreeRows, FreeRows] \ b[FreeRows]
    else
        # Solve in case no Neumann condition is defined
        x = A \ b
    end
    return x
end


function calculate_fields(coordinate, element, nodes2edge, edge2element, x)
    # Calculate displacement field u and flux field p based on solution vector x
    u = x[end-size(element, 1)+1:end, 1] # elementwise displacement field

    p = zeros(3 * size(element, 1), 2)
    for j in axes(element, 1)
        signum = ones(1, 3)
        dummy = diag(nodes2edge[element[j, [2 3 1]][:], element[j, [3 1 2]][:]])
        condition = j .== edge2element[dummy, 4]
        signum[findall(!iszero, condition)] .= -1
        c = coordinate[element[j, [2 3 1]][:], :] - coordinate[element[j, [3 1 2]][:], :]
        n = [norm(c[1, :]), norm(c[2, :]), norm(c[3, :])]
        coord = coordinate[element[j, :], :]'
        N = coord[:] * ones(1, 3) - repeat(coord, 3, 1)
        pc = reshape(N * diagm(vec(signum)) * diagm(n) * x[dummy]det([1 1 1; coordinate[element[j, :], :]']), 2, 3)
        p[3*(j-1).+[1, 2, 3], :] = [pc[1, :] pc[2, :]]
    end

    return u, p
end

function EBMfemSolver(coordinate, element, dirichlet, neumann, f, u_D, g)

    nodes2edge, edge2element, noedges = build_aux_structures(coordinate, element)
    A = assemble_global_stiffness_matrix(coordinate, element, nodes2edge, edge2element, noedges)
    b = calculate_constants_vector(coordinate, element, nodes2edge, noedges, f, u_D, dirichlet)

    x = solve_system(A, b, neumann, coordinate, element, g, noedges, nodes2edge)
    u, p = calculate_fields(coordinate, element, nodes2edge, edge2element, x)
    return x, u, p

end

end