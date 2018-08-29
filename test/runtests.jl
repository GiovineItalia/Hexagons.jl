using Hexagons
using Test

convert(HexagonOffsetOddR, HexagonAxial(2, 4))

x, y = center(HexagonAxial(2, 3))
@test isapprox(x,7.06217782649107)
@test isapprox(y,5.5)

h = cube_round(23.5, 4.67)

collect(vertices(HexagonAxial(2, 3)))
