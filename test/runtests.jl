using Hexagons
using Test
using Random: seed!
using Statistics: mean

# Test a few identities for the hexagon containing the point x, y
function run_point_test(x, y)
    hex_cubic = cube_round(x, y)
    verts = collect(vertices(hex_cubic))
    # x,y should be in the bounding box of the hexagon vertices
    @test minimum(v[1] for v in verts) <= x <= maximum(v[1] for v in verts)
    @test minimum(v[2] for v in verts) <= y <= maximum(v[2] for v in verts)
    # the center of the hexagon should be near its vertex mean
    mean_vert_x = mean(v[1] for v in verts)
    mean_vert_y = mean(v[2] for v in verts)
    hex_center = center(hex_cubic)
    @test isapprox(hex_center[1], mean_vert_x; atol = 1e-6)
    @test isapprox(hex_center[2], mean_vert_y; atol = 1e-6)
    # a string of type conversions should recover hex_cubic
    hex_axial = convert(HexagonAxial, hex_cubic)
    hex_offset = convert(HexagonOffsetOddR, hex_axial)
    other_hex_cubic = convert(HexagonCubic, hex_offset)
    @test other_hex_cubic == hex_cubic
end
    
test_points = [
    (0, 0),
    (1, 1),
    (1, -1),
    (0, 47),
    (-4.7, 0),
    (1.234, 5.678),
    (1e6, 2e6),
]
# bunch of random test points
seed!(1234)
for _ in 1:1000
    push!(test_points, (rand() * 100, rand() * 100))
end
for point in test_points
    run_point_test(point...)
end
