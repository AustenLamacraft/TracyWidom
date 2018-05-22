using TracyWidom
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end


# See https://arxiv.org/abs/0904.1581
@test F1(0) ≈ 0.83190806620295 atol=1e-14
@test F2(0) ≈ 0.96937282835526 atol=1e-14
