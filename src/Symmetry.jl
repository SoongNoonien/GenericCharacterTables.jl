function exceptions_apply(exceptions, vars, vals)
  for exception in exceptions
    if is_integer(evaluate(exception, vars, vals))
      return true
    end
  end
  return false
end

function parameter_matrix(char::GenericCharacter, val::Int64)
  q = parameters(parent(char))[1]
  params = parameters(char)
  exceptions = params.exceptions
  ranges = UnitRange{Int64}[]
  vars = Int64[1]
  for param in params
    push!(vars, var_index(param.var))
    upper = evaluate(param.modulus, [q], [val])
    push!(ranges, 1:numerator(constant_coefficient(upper)))
  end
  mat = Array{Union{Vector{GenericCyclo}, Nothing}}(undef, length.(ranges)...)
  for vals_with_index in Iterators.product(enumerate.(ranges)...)
    indices = first.(vals_with_index)
    vals = last.(vals_with_index)
    vals = collect((val, vals...))
    if exceptions_apply(exceptions, vars, vals)
      mat[indices...] = nothing
    else
      ordinary_char = evaluate.(char.values, Ref(vars), Ref(vals))
      mat[indices...] =  ordinary_char
    end
  end
  return mat
end

function reduce_parameter_matrix(mat::Array{Union{Vector{GenericCyclo}, Nothing}})
  reduced_mat = zeros(Int64, size(mat)...)
  reduced_val = 1
  found = true
  while found
    found = false
    val = nothing
    for i in CartesianIndices(reduced_mat)
      if iszero(reduced_mat[i])
        val = mat[i]
        if val !== nothing
          found = true
          break
        end
      end
    end
    if found
      for i in CartesianIndices(mat)
        if mat[i] == val
          reduced_mat[i] = reduced_val
        end
      end
      reduced_val += 1
    end
  end
  return reduced_mat
end

reduced_parameter_matrix(char::GenericCharacter, val::Int64) = reduce_parameter_matrix(parameter_matrix(char, val))
