using ..GenericCharacterTables
import ..GenericCharacterTables: Parameters, Parameter, CharTable, GenericCyclo, GenericCycloFrac
using Oscar
R = universal_polynomial_ring(QQ; cached=false)
q = gen(R, "q")
S = generic_cyclotomic_ring(R; congruence=QQ.((-1,3)))
a,b,m,n, _...=gens(R, ["a", "b", "m", "n", "a1", "b1", "m1", "n1", "a2", "b2", "m2", "n2", "a3", "b3", "m3", "n3", "at1", "bt1", "mt1", "nt1", "at2", "bt2", "mt2", "nt2"])

order = q^3*(q+1)^2*(q-1)*(q^2-q+1)*1//3
table = GenericCyclo[[
	S(1),
	S(1),
	S(1),
	S(1),
	S(1),
	S(1),
	S(1),
	S(1),
	S(1),
	S(1),
	S(1)] [
	(q^2-q)*S(1),
	(-q)*S(1),
	S(0),
	S(0),
	S(0),
	(1-q)*S(1),
	S(1),
	(2)*S(1),
	(2)*S(1),
	S(0),
	-S(1)] [
	(q^3)*S(1),
	S(0),
	S(0),
	S(0),
	S(0),
	(q)*S(1),
	S(0),
	-S(1),
	-S(1),
	S(1),
	-S(1)] [
	((q^2-q+1))*S(1),
	((1-q))*S(1),
	S(1),
	S(1),
	S(1),
	(1-q)*S(1, exponent=(3*n*a)*1//R((q+1)))+S(1, exponent=(-6*n*a)*1//R((q+1))),
	S(1, exponent=(3*n*a)*1//R((q+1)))+S(1, exponent=(-6*n*a)*1//R((q+1))),
	(3)*S(1),
	S(1, exponent=(3*n*a)*1//R((q+1)))+S(1, exponent=(3*n*b)*1//R((q+1)))+S(1, exponent=(3*n*(-a-b))*1//R((q+1))),
	S(1, exponent=(3*n*a)*1//R((q+1))),
	S(0)] [
	(q*(q^2-q+1))*S(1),
	(q)*S(1),
	S(0),
	S(0),
	S(0),
	(q-1)*S(1, exponent=(3*n*a)*1//R((q+1)))+q*S(1, exponent=(-6*n*a)*1//R((q+1))),
	-S(1, exponent=(3*n*a)*1//R((q+1))),
	(-3)*S(1),
	-(S(1, exponent=(3*n*a)*1//R((q+1)))+S(1, exponent=(3*n*b)*1//R((q+1)))+S(1, exponent=(3*n*(-a-b))*1//R((q+1)))),
	S(1, exponent=(3*n*a)*1//R((q+1))),
	S(0)] [
	((q+1)*(q^2-q+1))*S(1),
	S(1),
	S(1),
	S(1),
	S(1),
	(q+1)*S(1, exponent=(3*n*a)*1//R((q+1))),
	S(1, exponent=(3*n*a)*1//R((q+1))),
	S(0),
	S(0),
	S(1, exponent=(3*n*a)*1//R((q^2-1)))+S(1, exponent=(-3*n*a*q)*1//R((q^2-1))),
	S(0)] [
	((q-1)*(q^2-q+1)*1//3)*S(1),
	((2*q-1)*1//3)*S(1),
	(q-(q+1)*1//3)*S(1),
	(-(q+1)*1//3)*S(1),
	(-(q+1)*1//3)*S(1),
	(q-1)*S(1),
	-S(1),
	-(S(1, exponent=((q+1)*1//3)*1//R((3)))+S(1, exponent=(2*(q+1)*1//3)*1//R((3)))),
	-(S(1, exponent=(a-b)*1//R((3)))+S(1, exponent=(b-a)*1//R((3)))),
	S(0),
	S(0)] [
	((q-1)*(q^2-q+1)*1//3)*S(1),
	((2*q-1)*1//3)*S(1),
	(-(q+1)*1//3)*S(1),
	(q-(q+1)*1//3)*S(1),
	(-(q+1)*1//3)*S(1),
	(q-1)*S(1),
	-S(1),
	-(S(1, exponent=((q+1)*1//3)*1//R((3)))+S(1, exponent=(2*(q+1)*1//3)*1//R((3)))),
	-(S(1, exponent=(a-b)*1//R((3)))+S(1, exponent=(b-a)*1//R((3)))),
	S(0),
	S(0)] [
	((q-1)*(q^2-q+1)*1//3)*S(1),
	((2*q-1)*1//3)*S(1),
	(-(q+1)*1//3)*S(1),
	(-(q+1)*1//3)*S(1),
	(q-(q+1)*1//3)*S(1),
	(q-1)*S(1),
	-S(1),
	-(S(1, exponent=((q+1)*1//3)*1//R((3)))+S(1, exponent=(2*(q+1)*1//3)*1//R((3)))),
	-(S(1, exponent=(a-b)*1//R((3)))+S(1, exponent=(b-a)*1//R((3)))),
	S(0),
	S(0)] [
	((q-1)*(q^2-q+1))*S(1),
	((2*q-1))*S(1),
	-S(1),
	-S(1),
	-S(1),
	(q-1)*(S(1, exponent=(3*n*a-3*m*a)*1//R((q+1)))+S(1, exponent=(3*m*a-6*n*a)*1//R((q+1)))+S(1, exponent=(3*n*a)*1//R((q+1)))),
	-(S(1, exponent=(3*n*a-3*m*a)*1//R((q+1)))+S(1, exponent=(3*m*a-6*n*a)*1//R((q+1)))+S(1, exponent=(3*n*a)*1//R((q+1)))),
	-3*(S(1, exponent=(m)*1//R((3)))+S(1, exponent=(2*m)*1//R((3)))),
	-(S(1, exponent=(m*a+(3*n-m)*b)*1//R((q+1)))+S(1, exponent=(m*b+(3*n-m)*a)*1//R((q+1)))+S(1, exponent=(m*(-a-b)+(3*n-m)*a)*1//R((q+1)))+S(1, exponent=(m*a+(3*n-m)*(-a-b))*1//R((q+1)))+S(1, exponent=(m*(-a-b)+(3*n-m)*b)*1//R((q+1)))+S(1, exponent=(m*b+(3*n-m)*(-a-b))*1//R((q+1)))),
	S(0),
	S(0)] [
	((q+1)^2*(q-1))*S(1),
	(-(q+1))*S(1),
	-S(1),
	-S(1),
	-S(1),
	S(0),
	S(0),
	S(0),
	S(0),
	S(0),
	-(S(1, exponent=(3*n*a)*1//R((q^2-q+1)))+S(1, exponent=(-3*n*a*q)*1//R((q^2-q+1)))+S(1, exponent=(3*n*a*q^2)*1//R((q^2-q+1))))]]
classinfo = Vector{Any}[
	[[1,0],["^2A_2",[1,1,1]]],
	[[1,1],["^2A_2",[2,1]]],
	[[1,2],["^2A_2",[2,1]]],
	[[1,3],["^2A_2",[3]]],
	[[1,4],["^2A_2",[3]]],
	[[2,0],["A_1",[1,1]]],
	[[2,1],["A_1",[2]]],
	[[6,0],["A_0:3"[1]]],
	[[3,0],["A_0",[1]]],
	[[4,0],["A_0",[1]]],
	[[5,0],["A_0",[1]]]]
classlength = R.([
	1,
	(q+1)*(q-1)*(q^2-q+1),
	q*(q+1)^2*(q-1)*(q^2-q+1)*1//3,
	q*(q+1)^2*(q-1)*(q^2-q+1)*1//3,
	q*(q+1)^2*(q-1)*(q^2-q+1)*1//3,
	q^2*(q^2-q+1),
	q^2*(q+1)*(q-1)*(q^2-q+1),
	q^3*(q-1)*(q^2-q+1)*1//3,
	q^3*(q-1)*(q^2-q+1),
	q^3*(q+1)*(q^2-q+1),
	q^3*(q+1)^2*(q-1)])
charinfo = Vector{Any}[
	[[1,0],["^2A_2",[3]]],
	[[1,1],["^2A_2",[2,1]]],
	[[1,2],["^2A_2",[1,1,1]]],
	[[2,0],["A_1",[1,1]]],
	[[2,1],["A_1",[1,1]]],
	[[4,0],["A_0",[1]]],
	[[3,0,1],["A_0:3",[1]]],
	[[3,0,2],["A_0:3","\\omega"]],
	[[3,0,3],["A_0:3","\\omega^2"]],
	[[5,0],["A_0",[1]]],
	[[8,0],["A_0",[1]]]]
chardegree = R.([
	1,
	q^2-q,
	q^3,
	q^2-q+1,
	q*(q^2-q+1),
	(q+1)*(q^2-q+1),
	(q-1)*(q^2-q+1)*1//3,
	(q-1)*(q^2-q+1)*1//3,
	(q-1)*(q^2-q+1)*1//3,
	(q-1)*(q^2-q+1),
	(q+1)^2*(q-1)])

classsums=[
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	ss6=nesum(tt, a, 1, q)
	tt1=eesubs(tt, [a], [(q+1)*1//3])
	tt2=eesubs(tt, [a], [2*(q+1)*1//3])
	1//3*ss6-1//3*tt1-1//3*tt2
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	ss6=nesum(tt, a, 1, q)
	tt1=eesubs(tt, [a], [(q+1)*1//3])
	tt2=eesubs(tt, [a], [2*(q+1)*1//3])
	1//3*ss6-1//3*tt1-1//3*tt2
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	ss5=nesum(tt, a, 0, q)
	ss8=nesum(ss5, b, 0, q)
	tt1=eesubs(tt, [b], [a])
	ss9=nesum(tt1, a, 0, q)
	tt2=eesubs(tt1, [a], [(q+1)*1//3*a])
	ss10=nesum(tt2, a, 0, 2)
	tt5=eesubs(tt, [b], [a+(q+1)*1//3])
	tt4=eesubs(tt5, [a], [(q+1)*1//3*a])
	ss11=nesum(tt4, a, 0, 2)
	1//18*ss8-1//6*ss9+1//9*ss10-1//9*ss11
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	ss7=nesum(tt, a, 0, q^2-2)
	tt1=eesubs(tt, [a], [(q-1)*a])
	ss8=nesum(tt1, a, 0, q)
	1//6*ss7-1//6*ss8
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	ss10=nesum(tt, a, 0, q^2-q)
	tt1=eesubs(tt, [a], [(q^2-q+1)*1//3*a])
	ss11=nesum(tt1, a, 0, 2)
	1//9*ss10-1//9*ss11
end
]

charsums=[
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	nesum(tt, n, 1, (q-2)*1//3)
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	nesum(tt, n, 1, (q-2)*1//3)
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	ss7=nesum(tt, n, 0, (q^2-4)*1//3)
	tt1=eesubs(tt, [n], [(q-1)*n])
	ss8=nesum(tt1, n, 0, (q-2)*1//3)
	1//2*ss7-1//2*ss8
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	tt
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	ss5a=nesum(tt, m, 1, q)
	ss5=nesum(ss5a, n, 1, (q+1)*1//3)
	tt1=eesubs(tt, [n], [2*n])
	tt2=eesubs(tt1, [m], [3*n])
	ss6=nesum(tt2, n, 1, (q-2)*1//3)
	tt3=eesubs(tt, [m], [3*n])
	ss8=nesum(tt3, n, 1, (q-2)*1//3)
	ss7=eesubs(tt, [m,n], [(q+1)*1//3,(q+1)*1//3])
	1//6*ss5-1//6*ss6-1//3*ss7-1//6*ss8
end,
function (tt::Union{GenericCyclo, GenericCycloFrac})
	ss8=nesum(tt, n, 1, (q^2-q-2)*1//3)
	1//3*ss8
end
]

classparams=[
Parameters(Parameter[]),
Parameters(Parameter[]),
Parameters(Parameter[]),
Parameters(Parameter[]),
Parameters(Parameter[]),
Parameters([Parameter(a, q+1)], [((a)*1//((q+1)*1//3))]),
Parameters([Parameter(a, q+1)], [((a)*1//((q+1)*1//3))]),
Parameters(Parameter[]),
Parameters([Parameter(a, q+1), Parameter(b, q+1)], [((a)*1//((q+1)*1//3)), ((b)*1//((q+1)*1//3)), ((a-b)*1//(q+1))]),
Parameters([Parameter(a, q^2-1)], [((a)*1//(q-1))]),
Parameters([Parameter(a, q^2-q+1)], [((a)*1//((q^2-q+1)*1//3))])
]

charparams=[
Parameters(Parameter[]),
Parameters(Parameter[]),
Parameters(Parameter[]),
Parameters([Parameter(n, (q+1)*1//3)], [((n)*1//((q+1)*1//3))]),
Parameters([Parameter(n, (q+1)*1//3)], [((n)*1//((q+1)*1//3))]),
Parameters([Parameter(n, (q^2-1)*1//3)], [((n)*1//((q-1)*1//3))]),
Parameters(Parameter[]),
Parameters(Parameter[]),
Parameters(Parameter[]),
Parameters([Parameter(n, (q+1)*1//3), Parameter(m, q+1)], [((3*n-m)*1//(q+1)), ((m)*1//(q+1)), ((3*n-2*m)*1//(q+1))]),
Parameters([Parameter(n, (q^2-q+1)*1//3)], [((n)*1//((q^2-q+1)*1//3))])
]

classparamindex=var_index.([a,b])
charparamindex=var_index.([m,n])


information = raw"""- Information about the generic character table of $PSU_3(q)$,
  $q$ congruent to $2$ modulo $3$
  (See ``SU3.n2`` for the generic character table of $SL_3(q)$,
  $q$ not congruent to $2$ modulo $3$.)

- CHEVIE-name of the table: ``PSU3.2``

- The table was first computed in:
  {\sc W.~A.~Simpson and J.~S.~Frame}, The character tables for
  $SL(3,q)$, $SU(3,q^2)$, $PSL(3,q)$, $PSU(3,q^2)$,
  {\em Can.\ J.~Math.} {\bf 25} (1973), 486--494.

- The table was constructed by Jochen Gruber from the generic character
  table of $SL_3(q)$, $q$ congruent to $1$ modulo $3$.

- Corrections in:
  {\sc M.~Geck}, Irreducible Brauer characters of the $3$-dimensional
  special unitary groups in non-defining characteristic,
  {\em Comm.\ Algebra}, {\bf 18} (1990), 563--584.

- See also:
  {\sc V.~Ennola}, On the characters of the finite unitary groups,
  {\em Ann.\ Acad.\ Sci.\ Fenn.\ Ser.~A~I} {\bf 323} (1963), 1--36.

  {\sc V.~Ennola}, On the conjugacy classes of the finite unitary groups,
  {\em Ann.\ Acad.\ Sci.\ Fenn.\ Ser.~A~I} {\bf 313} (1962), 1--14.

  {\sc R.~Steinberg}, The representations of $GL(3,q)$, $GL(4,q)$, $PGL(3,q)$
  and $PGL(4,q)$, {\em Can.\ J.~Math.} {\bf 3} (1951), 225--235.
"""

TABLE=CharTable(order,permutedims(table),classinfo,classlength,charinfo,chardegree,
	classsums,charsums,classparamindex,charparamindex,classparams,charparams,S,information,splitext(basename(@__FILE__))[1])