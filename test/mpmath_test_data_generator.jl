using PyCall
@pyimport mpmath
@pyimport numpy

const scale = 3.0

function printer(M, N, functionstring, numbertype, splat="")
for m in M, n in N
  i = 0
  while i < 10
    a = numbertype.((rand(m) .- 0.5) * 2 * scale)
    b = numbertype.((rand(n) .- 0.5) * 2 * scale)
    c = numbertype.((rand() - 0.5) * 2 * scale)
    result = 0.0
    try
      result = numpy.float(mpmath.hyper(a, b, c))
      println("(a, b, c, result) = NumberType.($a), NumberType.($b), NumberType($(Float64(c))), $result")
      println("@test $functionstring(a$splat, b$splat, c) ≈ result atol=eps() rtol=rtol")
      i += 1
    catch
      continue
    end
  end
end
end

println("@testset \"mFn vs mpmath\" begin")
printer(1:3, 1:3, "mFn", Float64)
println("end")

println("@testset \"_₂F₁ vs mpmath\" begin")
printer(2, 1, "_₂F₁", Float64, "...")
println("end")

println("@testset \"_₃F₂ vs mpmath\" begin")
printer(3, 2, "_₃F₂", Float64, "...")
println("end")
