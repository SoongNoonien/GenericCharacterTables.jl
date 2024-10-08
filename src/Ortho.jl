import Oscar: norm, scalar_product, class_multiplication_coefficient

@doc raw"""
    class_multiplication_coefficient(t::CharTable, class1::Int64, class2::Int64, class3::Int64)

Return the (generic) class multiplication coefficient of the class types `class1`, `class2` and `class3` of the table `t`.

# Examples
```jldoctest
julia> g=generic_character_table("SL2.0");

julia> class_multiplication_coefficient(g,2,2,4)
q + 1
With exceptions:
  a3 ∈ (q + 1)ℤ
```
"""
function class_multiplication_coefficient(t::CharTable, class1::Int64, class2::Int64, class3::Int64)
	sum=0
	for char in t
		val1=shift_class_parameters(t, char[class1], 1)
		val2=shift_class_parameters(t, char[class2], 2)
		val3=shift_class_parameters(t, char[class3], 3)
		sum+=order(t)*charsum(char, val1*val2*conj(val3))//degree(char)  # TODO move order(t)* to the end of the function
	end
	return shrink((t.classlength[class1]*t.classlength[class2])*sum//(order(t)^2))
end

@doc raw"""
    class_multiplication_coefficient(t::SimpleCharTable{T}, class1::Int64, class2::Int64, class3::Int64) where T <: NfPoly

Return the (generic) class multiplication coefficient of the class types `class1`, `class2` and `class3` of the table `t`.

# Examples
```jldoctest
julia> g=green_function_table("GL3");

julia> class_multiplication_coefficient(g,1,2,3)
(q + 3)//(q^5 - 2*q^4 + q^3)
```
"""
function class_multiplication_coefficient(t::SimpleCharTable{T}, class1::Int64, class2::Int64, class3::Int64) where T <: NfPoly  # TODO is correct?
	sum=0
	for char in t
		sum1=char[class1]*char[class2]*char[class3]
		sum+=order(t)*sum1//degree(char)  # TODO move order(t)* to the end of the function
	end
	return (t.classlength[class1]*t.classlength[class2])*sum//(order(t)^2)
end

@doc raw"""
    norm(t::CharTable, class::Int64)

Return the (generic) norm of the class type `class`.

# Examples
```jldoctest
julia> g=generic_character_table("GL2");

julia> norm(g,2)
1
```
"""
function norm(t::CharTable, class::Int64)
	sum=0
	for char in t
		val=char[class]
		sum+=charsum(char, val*conj(val))
	end
	return shrink(t.classlength[class]*sum//order(t))
end

@doc raw"""
    norm(t::SimpleCharTable{T}, class::Int64) where T <: NfPoly

Return the (generic) norm of the class type `class`.

# Examples
```jldoctest
julia> g=green_function_table("GL3");

julia> norm(g,2)
(5*q^2 + 2*q + 3)//(q^5 - 2*q^4 + q^3)
```
"""
function norm(t::SimpleCharTable{T}, class::Int64) where T <: NfPoly  # TODO is correct?
	sum=0
	for char in t
		sum+=char[class]^2
	end
	return t.classlength[class]*sum//order(t)
end

@doc raw"""
    scalar_product(t::CharTable, class1::Int64, class2::Int64)

Return the (generic) scalar product between the class types `class1` and `class2`.

# Examples
```jldoctest
julia> g=generic_character_table("GL2");

julia> scalar_product(g,3,2)
0
With exceptions:
  i1 + j1 - 2*i2 ∈ (q - 1)ℤ
  i1 - i2 ∈ (q - 1)ℤ
  j1 - i2 ∈ (q - 1)ℤ
```
"""
function scalar_product(t::CharTable, class1::Int64, class2::Int64)
	sum=0
	for char in t
		val1=shift_class_parameters(t, char[class1], 1)
		val2=shift_class_parameters(t, char[class2], 2)
		sum+=charsum(char, val1*conj(val2))
	end
	return shrink(t.classlength[class1]*sum//order(t))
end

@doc raw"""
    scalar_product(t::SimpleCharTable{T}, class1::Int64, class2::Int64) where T <: NfPoly

Return the (generic) scalar product between the class types `class1` and `class2`.

# Examples
```jldoctest
julia> g=green_function_table("GL3");

julia> scalar_product(g,2,2)
(5*q^2 + 2*q + 3)//(q^5 - 2*q^4 + q^3)
```
"""
function scalar_product(t::SimpleCharTable{T}, class1::Int64, class2::Int64) where T <: NfPoly  # TODO is correct?
	sum=0
	for char in t
		sum+=char[class1]*char[class2]
	end
	return t.classlength[class1]*sum//order(t)
end
