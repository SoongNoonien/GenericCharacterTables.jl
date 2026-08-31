# GenericCharacterTables
# Copyright (C) 2022 Andreas Leim <leim@rhrk.uni-kl.de>, Martin Wagner <marwagne@rhrk.uni-kl.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

module GenericCharacterTables

using Oscar
import Oscar.AbstractAlgebra.Generic

const NfPoly = Union{PolyRingElem{QQFieldElem},PolyRingElem{AbsSimpleNumFieldElem}}

include("Imports.jl")
include("Types.jl")
include("GenericCyclotomics.jl")
include("ParameterExceptions.jl")
include("GenericCyclotomicFractions.jl")
include("Parameter.jl")
include("CharTable.jl")
include("Shifts.jl")
include("SumProc.jl")
include("GenericConjugacyClasses.jl")
include("GenericCharacter.jl")
include("Iteration.jl")
include("Congruence.jl")
include("Exports.jl")
include("Deprecations.jl")
include("CyclotomicFac.jl")
include("Symmetry.jl")

end
