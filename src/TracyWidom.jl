module TracyWidom

using SpecialFunctions, FastGaussQuadrature

export F2, F1

function F2(s::Real; num_points::Integer=25)
    quad = gausslegendre(num_points)
    _F2(s, quad)
end

function F2(s_vals::AbstractArray{T}; num_points::Integer=25) where {T<:Real}
    quad = gausslegendre(num_points)
    [_F2(s, quad) for s in s_vals]
end

function F1(s::Real; num_points::Integer=25)
    quad = gausslegendre(num_points)
    _F1(s, quad)
end

function F1(s_vals::AbstractArray{T}; num_points::Integer=25) where {T<:Real}
    quad = gausslegendre(num_points)
    [_F1(s, quad) for s in s_vals]
end

function _F2(s::Real, quad::Tuple{Array{T,1},Array{T,1}}) where {T<:Real}
    kernel = ((ξ,η) -> _K2tilde(ξ,η,s))
    _fredholm_det(kernel, quad)
end

function _F1(s::Real, quad::Tuple{Array{T,1},Array{T,1}}) where {T<:Real}
    kernel = ((ξ,η) -> _K1tilde(ξ,η,s))
    _fredholm_det(kernel, quad)
end

function _fredholm_det(kernel::Function, quad::Tuple{Array{T,1},Array{T,1}}) where {T<:Real}
    nodes, weights = quad
    N = length(nodes)
    sqrt_weights = sqrt.(weights)
    weights_matrix = kron(transpose(sqrt_weights),sqrt_weights)
    K_matrix = [kernel(ξ,η) for ξ in nodes, η in nodes]
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
_K2tilde(ξ,η,s) = sqrt(_ϕprime(ξ) * _ϕprime(η)) * _airy_kernel(_ϕ(ξ,s), _ϕ(η,s))

# For the GOE Tracy-Widom
_A_kernel(x,y) = airyai((x+y)/2) / 2
_K1tilde(ξ,η,s) = sqrt(_ϕprime(ξ) * _ϕprime(η)) * _A_kernel(_ϕ(ξ,s), _ϕ(η,s))

end # module
