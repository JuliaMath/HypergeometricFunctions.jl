using PyCall
@pyimport mpmath
@pyimport numpy

const scale = 3.0
const N = 10

for m in 1:3, n in 1:3
  @show m, n
  for i in 1:N
    a = (rand(m) .- 0.5) * 2 * scale
    b = (rand(n) .- 0.5) * 2 * scale
    c = (rand() - 0.5) * 2 * scale
    result = 0.0
    try
      result = numpy.float(mpmath.hyper(a, b, c))
      @show a, b, c, result
      println("@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol")
    catch
      continue
    end
  end
end

