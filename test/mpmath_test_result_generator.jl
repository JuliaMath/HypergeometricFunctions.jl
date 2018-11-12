using PyCall
@pyimport mpmath
@pyimport numpy

const scale = 3.0

function printer(M, N)
for m in M, n in N
  @show m, n
  i = 0
  while i < 10
    a = BigFloat.((rand(m) .- 0.5) * 2 * scale)
    b = BigFloat.((rand(n) .- 0.5) * 2 * scale)
    c = BigFloat.((rand() - 0.5) * 2 * scale)
    result = 0.0
    try
      result = numpy.float(mpmath.hyper(a, b, c))
      println("(a, b, c, result) = $a, $b, BigFloat($(Float64(c))), $result")
      println("@test mFn(a, b, c) ≈ result atol=eps() rtol=rtol")
      i += 1
    catch
      continue
    end
  end
end
end

@show "mFn"
printer(1:3, 1:3)

@show "2F1"
printer(2, 1)

@show "3F2"
printer(3, 2)
