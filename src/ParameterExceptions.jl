getindex(exceptions::ParameterExceptions, i::Integer) = exceptions.exceptions[i]

eltype(::Type{ParameterExceptions}) = UPolyFrac

length(exceptions::ParameterExceptions) = length(exceptions.exceptions)

iterate(exceptions::ParameterExceptions, state::Integer=1) =
  state > length(exceptions) ? nothing : (exceptions[state], state + 1)

@doc raw"""
    is_rational_number(x::UPolyFrac)

Return if `x` represents a rational number.
"""
is_rational_number(x::UPolyFrac) = is_constant(denominator(x)) && is_constant(numerator(x))

@doc raw"""
    is_integer_quotient(x::UPolyFrac, y::UPolyFrac)
    is_integer_quotient(p::UPoly, q::UPoly)

Return if `x//y` represents an integer.

Forming `x//y` would cancel the result via polynomial gcds, which dominates the
runtime of e.g. `scalar_product`. Instead note that `x//y` is an integer `n`
exactly if `x == n*y`, so `n` can be read off a single coefficient and only
needs to be verified.
"""
function is_integer_quotient(x::UPolyFrac, y::UPolyFrac)
  # The uncanonicalised numerators and denominators differ from the canonical
  # ones by a unit each, and those cancel in the quotient.
  nx, dx = numerator(x, false), denominator(x, false)
  ny, dy = numerator(y, false), denominator(y, false)
  dx == dy && return is_integer_quotient(nx, ny)
  return is_integer_quotient(nx * dy, dx * ny)
end

function is_integer_quotient(p::UPoly, q::UPoly)
  iszero(q) && throw(DivideError())
  iszero(p) && return true
  length(p) == length(q) || return false
  n = leading_coefficient(p)//leading_coefficient(q)
  isone(denominator(n)) || return false
  return p == n * q
end

@doc raw"""
    add_exception!(a::ParameterExceptions, exception::UPolyFrac)

Include `exception` into `a`. This also removes all now redundant exceptions from `a`.
"""
function add_exception!(a::ParameterExceptions, exception::UPolyFrac)
  @req !is_integer_quotient(numerator(exception), denominator(exception)) "The computation has failed, possibly due to illegal parameter combinations"
  is_integer_quotient(denominator(exception), numerator(exception)) && return nothing
  is_rational_number(exception) && return nothing
  new_exception = sign(leading_coefficient(numerator(exception))) * exception
  for old_exception in a.exceptions
    if is_integer_quotient(old_exception, new_exception)
      return nothing
    end
  end
  filter!(e -> !is_integer_quotient(new_exception, e), a.exceptions)
  push!(a.exceptions, new_exception)
end

@doc raw"""
    remove_exceptions!(a::ParameterExceptions, b::ParameterExceptions)

Remove the exceptions in `b` from `a`.
"""
function remove_exceptions!(a::ParameterExceptions, b::ParameterExceptions)
  setdiff!(a.exceptions, b.exceptions)
end

@doc raw"""
    merge(x::ParameterExceptions, y::ParameterExceptions)

Return a new collection of parameter exceptions composed of `x` and `y` where all redundant exceptions are omitted.
"""
function merge(x::ParameterExceptions, y::ParameterExceptions)
  is_restriction(x) || return y
  is_restriction(y) || return x
  exceptions = deepcopy(x)
  for exception in y.exceptions
    add_exception!(exceptions, exception)
  end
  return exceptions
end

@doc raw"""
    is_restriction(x::ParameterExceptions)

Return if `x` actually restricts something.
"""
is_restriction(x::ParameterExceptions) = !isempty(x)

function show(io::IO, x::ParameterExceptions)
  for (i, exception) in enumerate(x)
    if isone(denominator(exception))
      print(io, "$(numerator(exception)) ∈ ℤ")
    else
      print(io, "$(numerator(exception)) ∈ ($(denominator(exception)))ℤ")
    end
    if i < length(x)
      print(io, "\n")
    end
  end
end

# evaluate

function evaluate(x::ParameterExceptions, vars::Vector{Int64}, vals::Vector{<:RingElement})
  exceptions = parameter_exceptions()
  for exception in x
    add_exception!(exceptions, evaluate(exception, vars, vals))
  end
  return exceptions
end
