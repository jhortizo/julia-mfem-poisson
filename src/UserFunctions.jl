module UserFunctions
export fcn_zeros, fcn_ones, g_91

using LinearAlgebra


function fcn_zeros(x)
    zeros(size(x, 1), 1)
end

function fcn_ones(x)
    ones(size(x, 1), 1)
end

function g_91(x, n)
    a = angle.((x[:, 1] + x[:, 2] * im) * (-1 - im) / sqrt(2)) .+ pi * 3 / 4 
    r = sqrt.(x[:, 1] .^ 2 + x[:, 2] .^ 2)
    dot((2 / 3 * r .^ (-1 / 3) .* [-sin.(a / 3) cos.(a / 3)]), n) # TODO: solve dot product 
end

end