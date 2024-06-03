using ..GenericCharacterTables
import ..GenericCharacterTables: SimpleCharTable
using Oscar
R, q = polynomial_ring(QQ, "q")

order = q^10*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*(q-1)^2*(q+1)^5
table = QQPolyRingElem[[
	q^10,
	R(0),
	R(0),
	R(0),
	R(0),
	R(0),
	R(0),
	q^6,
	R(0),
	R(0),
	R(0),
	R(0),
	q^4,
	R(0),
	R(0),
	R(0),
	R(0),
	R(0),
	-q^3,
	R(0),
	R(0),
	q^3,
	R(0),
	R(0),
	q^2,
	R(0),
	R(0),
	R(0),
	q^2,
	R(0),
	-q,
	R(0),
	q,
	R(0),
	-q,
	R(0),
	R(1),
	R(-1),
	R(1),
	R(1),
	R(-1),
	R(-1),
	R(1)] [
	q^6*(q-1)*(q^2+1),
	-q^6,
	R(0),
	R(0),
	R(0),
	R(0),
	R(0),
	-(q-1)*(q^2+1)*q^3,
	q^3,
	R(0),
	R(0),
	R(0),
	-q^2*(q-1)^2,
	q^3,
	-q^2,
	R(0),
	R(0),
	R(0),
	q*(2*q^2-q+1),
	q,
	R(0),
	(q-1)*q,
	-q,
	R(0),
	-2*(q-1)*q,
	q,
	q,
	R(0),
	R(0),
	R(0),
	3*q-1,
	R(-1),
	-q+1,
	R(1),
	R(-1),
	R(-1),
	R(-4),
	R(2),
	R(0),
	R(-1),
	R(-1),
	R(0),
	R(1)] [
	(q^4-q^3+q^2-q+1)*q^4,
	-q^4*(q-1),
	q^4,
	R(0),
	R(0),
	R(0),
	R(0),
	-q^2*(q^3-2*q^2+q-1),
	-q^2*(q-1),
	q^2,
	R(0),
	R(0),
	q*(q^3-q^2+2*q-1),
	(q-1)*q,
	(q-1)*q,
	-q,
	R(0),
	R(0),
	-q*(q^2-2*q+2),
	-2*q,
	R(0),
	q^3,
	R(0),
	R(0),
	2*q^2-2*q+1,
	-q+1,
	-q+1,
	R(1),
	R(1),
	R(1),
	-3*q+2,
	R(2),
	q,
	R(0),
	R(-1),
	R(-1),
	R(5),
	R(-1),
	R(1),
	R(-1),
	R(-1),
	R(1),
	R(0)] [
	q^3*(q^2-q+1)*(q^2+1),
	(q^2-q+1)*q^3,
	q^3,
	q^3,
	R(0),
	R(0),
	R(0),
	q*(q^2-q+1)*(q^2+1),
	q*(q^2-q+1),
	q,
	q,
	R(0),
	2*q*(q^2-q+1),
	q*(q^2-q+1),
	-q*(q-2),
	q,
	q,
	R(0),
	(q-1)*(q^2-q+1),
	2*q-1,
	R(-1),
	(q+1)*(q^2-q+1),
	R(1),
	R(1),
	-q^2+4*q-1,
	2*q-1,
	2*q-1,
	R(-1),
	q^2+1,
	R(1),
	3*q-3,
	R(-3),
	q+1,
	R(1),
	R(0),
	R(0),
	R(-6),
	R(0),
	R(2),
	R(0),
	R(0),
	R(0),
	R(-1)] [
	q^2*(q^4-q^3+q^2-q+1),
	(q^2-q+1)*q^2,
	-q^2*(q-1),
	q^2,
	q^2,
	R(0),
	R(0),
	q*(q^3-q^2+2*q-1),
	q*(2*q-1),
	(q-1)*q,
	-q,
	R(0),
	-q^3+2*q^2-q+1,
	q^2-q+1,
	q^2-q+1,
	-q+1,
	R(1),
	R(1),
	2*q^2-2*q+1,
	-2*q+1,
	R(1),
	R(1),
	R(1),
	R(1),
	q^2-2*q+2,
	-q+2,
	-q+2,
	R(2),
	q^2,
	R(0),
	-2*q+3,
	R(3),
	R(1),
	R(1),
	q,
	R(0),
	R(5),
	R(1),
	R(1),
	R(-1),
	R(1),
	R(-1),
	R(0)] [
	q*(q-1)*(q^2+1),
	-q*(q^2-q+1),
	(q-1)*q,
	(q-1)*q,
	-q,
	-q,
	R(0),
	-(q-1)*(q^2+1),
	q^2-q+1,
	-q+1,
	-q+1,
	R(1),
	(q-1)^2,
	q^2-q+1,
	-2*q+1,
	-q+1,
	-q+1,
	R(1),
	q^2-q+2,
	-q+2,
	R(2),
	(q-1)*q,
	-q,
	R(0),
	-2*q+2,
	-q+2,
	-q+2,
	R(2),
	R(0),
	R(0),
	-q+3,
	R(3),
	-q+1,
	R(1),
	-q,
	R(0),
	R(4),
	R(2),
	R(0),
	R(1),
	R(-1),
	R(0),
	R(-1)] [
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1),
	R(1)]]
classinfo = Vector{Any}[
	["^2A_4(q)",q+1,[[1,1,1,1,1]]],
	["^2A_4(q)",q+1,[[2,1,1,1]]],
	["^2A_4(q)",q+1,[[2,2,1]]],
	["^2A_4(q)",q+1,[[3,1,1]]],
	["^2A_4(q)",q+1,[[3,2]]],
	["^2A_4(q)",q+1,[[4,1]]],
	["^2A_4(q)",q+1,[[5]]],
	["^2A_3(q)",(q+1)^2,[[1,1,1,1]]],
	["^2A_3(q)",(q+1)^2,[[2,1,1]]],
	["^2A_3(q)",(q+1)^2,[[2,2]]],
	["^2A_3(q)",(q+1)^2,[[3,1]]],
	["^2A_3(q)",(q+1)^2,[[4]]],
	["^2A_2(q)+A_1(q)",(q+1)^2,[[1,1,1],[1,1]]],
	["^2A_2(q)+A_1(q)",(q+1)^2,[[1,1,1],[2]]],
	["^2A_2(q)+A_1(q)",(q+1)^2,[[2,1],[1,1]]],
	["^2A_2(q)+A_1(q)",(q+1)^2,[[2,1],[2]]],
	["^2A_2(q)+A_1(q)",(q+1)^2,[[3],[1,1]]],
	["^2A_2(q)+A_1(q)",(q+1)^2,[[3],[2]]],
	["^2A_2(q)",(q+1)^3,[[1,1,1]]],
	["^2A_2(q)",(q+1)^3,[[2,1]]],
	["^2A_2(q)",(q+1)^3,[[3]]],
	["^2A_2(q)",(q-1)*(q+1)^2,[[1,1,1]]],
	["^2A_2(q)",(q-1)*(q+1)^2,[[2,1]]],
	["^2A_2(q)",(q-1)*(q+1)^2,[[3]]],
	["A_1(q)+A_1(q)",(q+1)^3,[[1,1],[1,1]]],
	["A_1(q)+A_1(q)",(q+1)^3,[[1,1],[2]]],
	["A_1(q)+A_1(q)",(q+1)^3,[[2],[1,1]]],
	["A_1(q)+A_1(q)",(q+1)^3,[[2],[2]]],
	["A_1(q^2)",(q-1)*(q+1)^2,[[1,1]]],
	["A_1(q^2)",(q-1)*(q+1)^2,[[2]]],
	["A_1(q)",(q+1)^4,[[1,1]]],
	["A_1(q)",(q+1)^4,[[2]]],
	["A_1(q)",(q-1)*(q+1)^3,[[1,1]]],
	["A_1(q)",(q-1)*(q+1)^3,[[2]]],
	["A_1(q)",(q^2-q+1)*(q+1)^2,[[1,1]]],
	["A_1(q)",(q^2-q+1)*(q+1)^2,[[2]]],
	["\\emptyset",(q+1)^5,[[[1]],1]],
	["\\emptyset",(q-1)*(q+1)^4,[[[1]],1]],
	["\\emptyset",(q-1)^2*(q+1)^3,[[[1]],1]],
	["\\emptyset",(q^2-q+1)*(q+1)^3,[[[1]],1]],
	["\\emptyset",(q-1)*(q^2-q+1)*(q+1)^2,[[[1]],1]],
	["\\emptyset",(q-1)*(q^2+1)*(q+1)^2,[[[1]],1]],
	["\\emptyset",(q+1)*(q^4-q^3+q^2-q+1),[[[1]],1]]]
classlength = R.([
	1,
	q^8-q^7+q^6-q^5+q^3-q^2+q-1,
	q^12-q^11+q^10-q^8+2*q^7-2*q^6+q^5-q^3+q^2-q,
	q^14-q^13+q^12-q^10+2*q^9-2*q^8+q^7-q^5+q^4-q^3,
	q^16-q^15+q^13-2*q^12+2*q^11-q^10-q^9+2*q^8-2*q^7+q^6-q^4+q^3,
	q^18-q^17+q^15-2*q^14+2*q^13-q^12-q^11+2*q^10-2*q^9+q^8-q^6+q^5,
	q^20-q^18+q^17-q^16+q^14-2*q^13+q^12-q^10+q^9-q^8+q^6,
	(q^4-q^3+q^2-q+1)*q^4,
	(q^6-q^5+q^4-q^2+q-1)*(q^4-q^3+q^2-q+1)*q^4,
	(q^8+q^5-q^4-q)*(q^4-q^3+q^2-q+1)*q^4,
	(q^10-q^9+q^7-2*q^6+q^5-q^3+q^2)*(q^4-q^3+q^2-q+1)*q^4,
	(q^12-q^10+q^9-q^8-q^7+q^6-q^5+q^3)*(q^4-q^3+q^2-q+1)*q^4,
	(q^2+1)*(q^4-q^3+q^2-q+1)*q^6,
	(q^2-1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^6,
	(q^4-q^3+q-1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^6,
	(q^6-q^5-q^4+2*q^3-q^2-q+1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^6,
	(q^6-q^4+q^3-q)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^6,
	(q^8-2*q^6+q^5+q^4-2*q^3+q)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^6,
	(q-1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^7,
	(q^4-q^3+q-1)*(q-1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^7,
	(q^6-q^4+q^3-q)*(q-1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^7,
	(q+1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^7,
	(q^4-q^3+q-1)*(q+1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^7,
	(q^6-q^4+q^3-q)*(q+1)*(q^2+1)*(q^4-q^3+q^2-q+1)*q^7,
	(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^8,
	(q^2-1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^8,
	(q^2-1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^8,
	(q^4-2*q^2+1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^8,
	(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q+1)^2*q^8,
	(q^4-1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q+1)^2*q^8,
	(q-1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^9,
	(q^2-1)*(q-1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^9,
	(q+1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^9,
	(q^2-1)*(q+1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^9,
	(q-1)*(q^2+1)*(q^4-q^3+q^2-q+1)*(q+1)^2*q^9,
	(q^2-1)*(q-1)*(q^2+1)*(q^4-q^3+q^2-q+1)*(q+1)^2*q^9,
	(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*(q-1)^2*q^10,
	(q-1)*(q+1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*q^10,
	(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q^2+1)*(q+1)^2*q^10,
	(q^2+1)*(q^4-q^3+q^2-q+1)*(q-1)^2*(q+1)^2*q^10,
	(q-1)*(q^2+1)*(q^4-q^3+q^2-q+1)*(q+1)^3*q^10,
	(q-1)*(q^2-q+1)*(q^4-q^3+q^2-q+1)*(q+1)^3*q^10,
	(q^2-q+1)*(q^2+1)*(q-1)^2*(q+1)^4*q^10])
classtypeorder = R.([
	q+1,
	 q+1,
	 q+1,
	 q+1,
	 q+1,
	 q+1,
	 q+1,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 q^2+q,
	 1//2*q^3-1//2*q,
	 1//2*q^3-1//2*q,
	 1//2*q^3-1//2*q,
	 1//2*q^3-3//2*q-1,
	 1//2*q^3-3//2*q-1,
	 1//2*q^3-3//2*q-1,
	 1//2*q^3-1//2*q,
	 1//2*q^3-1//2*q,
	 1//2*q^3-1//2*q,
	 1//2*q^3-1//2*q,
	 1//2*q^3-3//2*q-1,
	 1//2*q^3-3//2*q-1,
	 1//6*q^4-1//3*q^3-1//6*q^2+1//3*q,
	 1//6*q^4-1//3*q^3-1//6*q^2+1//3*q,
	 1//2*q^4-3//2*q^2-q,
	 1//2*q^4-3//2*q^2-q,
	 1//3*q^4+1//3*q^3-1//3*q^2-1//3*q,
	 1//3*q^4+1//3*q^3-1//3*q^2-1//3*q,
	 1//120*q^5-1//24*q^4+1//24*q^3+1//24*q^2-1//20*q,
	 1//12*q^5-1//12*q^4-1//4*q^3+1//12*q^2+1//6*q,
	 1//8*q^5-1//8*q^4-7//8*q^3+1//8*q^2+7//4*q+1,
	 1//6*q^5+1//6*q^4-1//6*q^3-1//6*q^2,
	 1//6*q^5-1//6*q^4-1//2*q^3+1//6*q^2+1//3*q,
	 1//4*q^5+1//4*q^4-1//4*q^3-1//4*q^2,
	 1//5*q^5-1//5*q])
charinfo = Vector{Any}[[[1,1,1,1,1]], [[2,1,1,1]], [[2,2,1]], [[3,1,1]], [[3,2]], [[4,1]], [[5]]]
chardegree = R.([
	q^10,
	q^6*(q-1)*(q^2+1),
	(q^4-q^3+q^2-q+1)*q^4,
	q^3*(q^2-q+1)*(q^2+1),
	q^2*(q^4-q^3+q^2-q+1),
	q*(q-1)*(q^2+1),
	1])

information = "- Information about the tables of unipotent characters for GU_5(q).

- CHEVIE-name of the table: uniGU5

- This table is computed with general programs written by F.Luebeck.
  They compute the Deligne-Lusztig characters R_T^G(1) and find the
  unipotent characters as linear combinations of them.
"

TABLE=SimpleCharTable(order,permutedims(table),classinfo,classlength,classtypeorder,charinfo,chardegree,R,information)
