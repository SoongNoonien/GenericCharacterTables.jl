```@meta
CurrentModule = GenericCharacterTables
DocTestSetup = :(using GenericCharacterTables, Oscar)
```

# An extended example

The following is based on an example in [OSCAR](@cite),
which in turn is based on Section 5.3 of [GHLMP96](@cite).
To fully understand what is going on we strongly recommend to read the latter
source concurrently with the present text.

To start we load the generic character table for $\mathrm{SL}_3(q)$ with
$q\not\equiv 1\pmod 3$.
To learn about the origin of the table, we could enter `printinfotab(T)`.
```jldoctest book
julia> T = genchartab("SL3.n1")
Generic character table SL3.n1
  of order q^8 - q^6 - q^5 + q^3
  with 8 irreducible character types
  with 8 class types
  with parameters (a, b, m, n)
```

As we can see this table has four parameters in addition to $q$.
The entries of the table are “generalized cyclotomics”, that is,
linear combinations over $\mathbb{Q}(q)$ of symbolic “roots of unity”
depending on the parameters listed above.
```jldoctest book
julia> T[4,4]
(q + 1)*exp(2π𝑖((a*n)//(q - 1))) + exp(2π𝑖((-2*a*n)//(q - 1)))
```

Denoting row 4 of the table by $\chi_4$, we note that is not a single
character, but a *character type*, describing a whole family of
characters.
```jldoctest book
julia> χ₄ = T[4]
Generic character of SL3.n1
  with parameters
    n ∈ {1,…, q - 1} except n ∈ (q - 1)ℤ
  of degree q^2 + q + 1
  with values
    q^2 + q + 1
    q + 1
    1
    (q + 1)*exp(2π𝑖((a*n)//(q - 1))) + exp(2π𝑖((-2*a*n)//(q - 1)))
    exp(2π𝑖((a*n)//(q - 1))) + exp(2π𝑖((-2*a*n)//(q - 1)))
    exp(2π𝑖((a*n)//(q - 1))) + exp(2π𝑖((b*n)//(q - 1))) + exp(2π𝑖((-a*n - b*n)//(q - 1)))
    exp(2π𝑖((a*n)//(q - 1)))
    0
```

With these generic characters, we can compute norms, scalar products, and more.
For this demonstration we tensor the second character type $\chi_2$ with itself.
```jldoctest book
julia> h = T[2] * T[2]
Generic character of SL3.n1
  of degree q^4 + 2*q^3 + q^2
  with values
    q^4 + 2*q^3 + q^2
    q^2
    0
    q^2 + 2*q + 1
    1
    4
    0
    1
```

We may now attempt to decompose this character type `h` by computing its scalar product with the irreducible character types.
This returns a “generic” scalar product, plus a (possibly empty)
list of parameter exceptions for which the general result may not hold.
For example:
```jldoctest book
julia> scalar(T[4], h)
0
With exceptions:
  2*n1 ∈ (q - 1)ℤ
```
This scalar product is $0$
except when $q-1$ divides $2n_1$, where $n_1$ indicates the value of the parameter $n$ for the first factor in the scalar product (i.e. the character type $\chi_4$).
For a “generic decomposition” we need to compute all the scalar products and
exceptions.
```jldoctest book
julia> for i in 1:8 println("<$i, h> = ", scalar(T[i], h)) end
<1, h> = 1
<2, h> = 2
<3, h> = 2
<4, h> = 0
With exceptions:
  2*n1 ∈ (q - 1)ℤ
<5, h> = 0
With exceptions:
  2*n1 ∈ (q - 1)ℤ
<6, h> = 0
With exceptions:
  2*m1 - n1 ∈ (q - 1)ℤ
  m1 - 2*n1 ∈ (q - 1)ℤ
  m1 + n1 ∈ (q - 1)ℤ
  m1 ∈ (q - 1)ℤ
  m1 - n1 ∈ (q - 1)ℤ
  n1 ∈ (q - 1)ℤ
<7, h> = 0
With exceptions:
  n1 ∈ (q - 1)ℤ
<8, h> = 0
With exceptions:
  q*n1 ∈ (q^2 + q + 1)ℤ
  n1 ∈ (q^2 + q + 1)ℤ
  q*n1 + n1 ∈ (q^2 + q + 1)ℤ
```

This suggest a decomposition of $\chi_2\otimes\chi_2$ into $\chi_1+2\chi_2+2\chi_3$ “in general”.
But comparing the respective degrees, we notice a discrepancy:
```jldoctest book
julia> degree(lincomb([1,2,2],[T[1],T[2],T[3]]))
2*q^3 + 2*q^2 + 2*q + 1

julia> degree(h)
q^4 + 2*q^3 + q^2
```

To resolve this, we need to work through the exceptions.
Recall that earlier we saw that $\langle\chi_4,\chi_2^2\rangle=0$
except when $q-1$ divides $2n_1$, where $n_1$ was the value of the parameter $n$
used in $\chi_4$. Further restrictions apply to $n$:
```jldoctest book
julia> printcharparam(T,4)
4	n ∈ {1,…, q - 1} except n ∈ (q - 1)ℤ
```

So $n$ may take on any value between $1$ and $q-1$ not divisible by $q-1$.
Hence the only possible exception is $n=(q-1)/2$ which can only
occur if $q$ is odd.
It thus makes sense to consider $q$ odd and $q$ even separately.

We demonstrate this for $q$ even. Then $\langle\chi_4,\chi_2^2\rangle=0$
and with a similar argument $\langle\chi_5,\chi_2^2\rangle=0$.
We now construct a copy of the table but with the congruence equation $q\equiv 0\pmod 2$ applied:
```jldoctest book
julia> T2 = setcongruence(T, (0,2));

```
Inspecting the list of exceptions for $\langle\chi_6,\chi_2^2\rangle$, the first
occurs when $q-1$ divides $m+n$. To study this case we specialize $m$:
```jldoctest book
julia> (q, (a, b, m, n)) = params(T2);

julia> x = param(T2, "x");  # create an additional "free" variable

julia> speccharparam!(T2, 6, m, -n + (q-1)*x)  # force m = -n (mod q-1)
```

Recomputing the scalar product gives a new result
```jldoctest book
julia> h = T2[2] * T2[2];

julia> scalar(T2[6], h)
1
With exceptions:
  2*n1 ∈ (q - 1)ℤ
  3*n1 ∈ (q - 1)ℤ
```

The exceptions both cannot occur as $q$ is even and the table we are considering
is only defined for $q\not\equiv 1\pmod 3$.
By working through the other possible exceptions and irreducible character types, and
handling duplicates, one finally obtains
```math
\chi_2^2 = \chi_1+2\chi_2+2\chi_3
    +\frac12\sum_{n=1}^{q-2} \chi_6(n,q-1-n)
    +\frac12\sum_{n=1}^{q} \chi_7(n(q-1)).
```
Where $\chi_6(n,q-1-n)$ shows that the $6$th character in the table `T2`
is a family on two parameters: $n$ and $q-1-n$, while $\chi_7$ depends on only
one, namely $n(q-1)$.
A similar result can be obtained for odd $q$ albeit with a few more cases that need
to be dealt with, but all in essentially the same manner.
