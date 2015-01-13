
module Hexagons

using Compat

import Base: convert, start, next, done, length

export HexagonAxial, HexagonCubic, HexagonOffsetOddR, HexagonOffsetEvenR,
       center, points, HexPoints, hexpoints, pointhex, neighbors


# Various ways to index hexagons in a grid
# ----------------------------------------

abstract Hexagon

immutable HexagonAxial <: Hexagon
    q::Int
    r::Int
end


immutable HexagonCubic <: Hexagon
    x::Int
    y::Int
    z::Int
end


immutable HexagonOffsetOddR <: Hexagon
    q::Int
    r::Int
end


immutable HexagonOffsetEvenR <: Hexagon
    q::Int
    r::Int
end


# Neighbor hexagon iterator
# -------------------------

immutable HexagonCubicNeighborIterator
    hex::HexagonCubic
end

const cubic_hex_neighbor_offsets =
    [ 1 -1  0
      1  0 -1
      0  1 -1
     -1  1  0
     -1  0  1
      0 -1  1 ]


function neigbors(hex::HexagonCubic)
    HexagonCubicNeighborIterator(hex)
end


function start(it::HexagonCubicNeighborIterator)
    return 1
end


function next(it::HexagonCubicNeighborIterator, state)
    neighbor = HexagonCubic(it.hex.x + cubic_hex_neighbor_offsets[state, 1],
                            it.hex.y + cubic_hex_neighbor_offsets[state, 2],
                            it.hex.z + cubic_hex_neighbor_offsets[state, 3])
    return (neighbor, state+1)
end


function done(it::HexagonCubicNeighborIterator, state)
    return state > size(cubic_hex_neighbor_offsets, 1)
end


# Convert between hexagon indexing
# --------------------------------


function convert(::Type{HexagonAxial}, hex::HexagonCubic)
    return HexagonAxial(hex.x, hex.z)
end


function convert(::Type{HexagonAxial}, hex::HexagonOffsetOddR)
    return convert(HexagonAxial, convert(HexagonCubic, hex))
end


function convert(::Type{HexagonAxial}, hex::HexagonOffsetEvenR)
    return convert(HexagonAxial, convert(HexagonCubic, hex))
end


function convert(::Type{HexagonCubic}, hex::HexagonAxial)
    return HexagonCubic(hex.q, hex.r, -hex.q - hex.r)
end


function convert(::Type{HexagonCubic}, hex::HexagonOffsetOddR)
    x = hex.q - div(hex.r - (hex.r & 1), 2)
    z = hex.r
    y = -x - z
    return HexagonCubic(x, y, z)
end


function convert(::Type{HexagonCubic}, hex::HexagonOffsetEvenR)
    x = hex.q - div(hex.r + (hex.r & 1), 2)
    z = hex.r
    y = -x - z
    return HexagonCubic(x, y, z)
end


function convert(::Type{HexagonOffsetOddR}, hex::HexagonCubic)
    q = hex.x + div(hex.z - (hex.z & 1), 2)
    r = hex.z
    return HexagonOffsetOddR(q, r)
end


function convert(::Type{HexagonOffsetOddR}, hex::HexagonAxial)
    return convert(HexagonOffsetOddR, convert(HexagonCubic, hex))
end


# Conversion between euclidean and hexagon coordinates
# ----------------------------------------------------

function center(hex::HexagonAxial, xsize=1.0, ysize=1.0, xoff=1.0, yoff=1.0)
    (xoff + xsize * sqrt(3) * (hex.q + hex.r/2), yoff + ysize * (3/2) * hex.r)
end


function center(hex::Hexagon, xsize=1.0, ysize=1.0, xoff=1.0, yoff=1.0)
    center(convert(HexagonAxial, hex), xsize, ysize, xoff, yoff)
end


immutable HexPoints
    x_center::Float64
    y_center::Float64
    xsize::Float64
    ysize::Float64

    function HexPoints(x, y, xsize=1.0, ysize=1.0)
        return new(float64(x), float64(y), float64(xsize), float64(ysize))
    end
end


function length(it::HexPoints)
    return 6
end


function points(hex::Hexagon, xsize=1.0, ysize=1.0, xoff=0.0, yoff=0.0)
    c = center(hex, size, xoff, yoff)
    return HexPoints(c[1], c[2], size)
end


function hexpoints(x, y, xsize=1.0, ysize=1.0)
    collect((Float64, Float64),
            HexPoints(float64(x), float64(y),
                      float64(xsize), float64(ysize)))
end


function start(it::HexPoints)
    return 1
end


function next(it::HexPoints, state)
    theta = 2*pi/6 * (state-1+0.5)
    x_i = it.x_center + it.xsize * cos(theta)
    y_i = it.y_center + it.ysize * sin(theta)
    return ((x_i, y_i), state+1)
end


function done(it::HexPoints, state)
    return state > 6
end


# Find the nearest hexagon in cubic coordinates.
function nearest_cubic_hexagon(x::Real, y::Real, z::Real)
    rx, ry, rz = @compat round(Integer, x), round(Integer, y), round(Integer, z)
    x_diff, y_diff, z_diff = abs(rx - x), abs(ry - y), abs(rz - z)

    if x_diff > y_diff && x_diff > z_diff
        rx = -ry - rz
    elseif y_diff > z_diff
        ry = -rx - rz
    else
        rz = -rx - ry
    end

    return HexagonCubic(rx, ry, rz)
end


# Return the index (in axial coordinates) of the hexagon containing the
# point x, y
function pointhex(x, y, xsize=1.0, ysize=1.0)
    x /= xsize
    y /= ysize
    q = sqrt(3)/3 * x - y/3
    r = 2 * y / 3
    h = nearest_cubic_hexagon(q, -q - r, r)
    #return h

    x0, y0 = center(h)
    d0 = (x0-x)^2 + (y0-y)^2
    h_best = h
    d_best = d0
    for neighbor in neigbors(h)
        xn, yn = center(neighbor)
        dn = (xn-x)^2 + (yn-y)^2
        if dn < d_best
            d_best = dn
            h_best = neighbor
        end
    end
    return h_best
end



end # module Hexagons

