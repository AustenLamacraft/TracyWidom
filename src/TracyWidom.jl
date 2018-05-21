module TracyWidom

using SpecialFunctions, FastGaussQuadrature

export F2

function F2(s::Real, N::Integer)
    nodes, weights = gausslegendre(N)
    sqrt_weights = sqrt.(weights)
    weights_matrix = kron(transpose(sqrt_weights),sqrt_weights)
    K_matrix = [K2tilde(ξ,η,s) for ξ in nodes, η in nodes]
    det(eye(N) - weights_matrix .* K_matrix)
end

function _airy_kernel(x, y)
    if x==y
        return (airyaiprime(x))^2 - x * (airyai(x))^2
    else
        return (airyai(x) * airyaiprime(y) - airyai(y) * airyaiprime(x)) / (x - y)
    end
end

_ϕ(ξ, s) =  s + 10*tan(π*(ξ+1)/4)
_ϕprime(ξ) = (5π/2)*(sec(π*(ξ+1)/4))^2
_K2tilde(ξ,η,s) = sqrt(ϕprime(ξ) * ϕprime(η)) * airy_kernel(ϕ(ξ,s), ϕ(η,s));

end # module
