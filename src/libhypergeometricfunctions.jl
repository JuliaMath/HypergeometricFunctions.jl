if get(ENV, "HF_BUILD_FROM_SOURCE", "false") == "true"
    using Libdl
    const libhypergeometricfunctions = find_library("libhypergeometricfunctions", [joinpath(dirname(@__DIR__), "deps")])
    if libhypergeometricfunctions ≡ nothing || length(libhypergeometricfunctions) == 0
        error("HypergeometricFunctions is not properly installed. Please run Pkg.build(\"HypergeometricFunctions\") ",
              "and restart Julia.")
    end
else
    #using HypergeometricFunctions_jll
end

function drummondpFq(::Tuple{}, ::Tuple{}, z::Float64; kmax::Int=10_000)
    return ccall((:hf_drummond0f0, libhypergeometricfunctions), Cdouble, (Cdouble, Cint), z, kmax)
end

function drummondpFq(::Tuple{}, ::Tuple{}, z::BigFloat; kmax::Int=10_000)
    ret = BigFloat()
    ccall((:hf_drummond0f0_mpfr, libhypergeometricfunctions), Cvoid, (Ref{BigFloat}, Ref{BigFloat}, Cint, Clong, Int32), ret, z, kmax, precision(BigFloat), Base.MPFR.ROUNDING_MODE[])
    return ret
end

function wenigerpFq(::Tuple{}, ::Tuple{}, z::Float64; kmax::Int=10_000)
    return ccall((:hf_weniger0f0, libhypergeometricfunctions), Cdouble, (Cdouble, Cint), z, kmax)
end

function wenigerpFq(α::Tuple{Float64}, ::Tuple{}, z::Float64; kmax::Int=10_000)
    return ccall((:hf_weniger1f0, libhypergeometricfunctions), Cdouble, (Cdouble, Cdouble, Cint), α[1], z, kmax)
end

function wenigerpFq(α::Tuple{Float64, Float64}, ::Tuple{}, z::Float64; kmax::Int=10_000)
    return ccall((:hf_weniger2f0, libhypergeometricfunctions), Cdouble, (Cdouble, Cdouble, Cdouble, Cint), α[1], α[2], z, kmax)
end

function wenigerpFq(α::Tuple{Float64, Float64}, β::Tuple{Float64}, z::Float64; kmax::Int=10_000)
    return ccall((:hf_weniger2f1, libhypergeometricfunctions), Cdouble, (Cdouble, Cdouble, Cdouble, Cdouble, Cint), α[1], α[2], β[1], z, kmax)
end
