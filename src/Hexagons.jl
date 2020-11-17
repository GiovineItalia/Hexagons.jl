module Hexagons

import Base: convert, length, collect, iterate

export HexagonAxial, HexagonCubic, HexagonOffsetOddR, HexagonOffsetEvenR,
       hexagon, center, hexpoints, cube_round, cube_linedraw, neighbor,
       HexagonVertexIterator, vertices,
       HexagonNeighborIterator, neighbors,
       HexagonDiagonalIterator, diagonals,
       HexagonDistanceIterator, hexagons_within,
       HexagonRingIterator, ring,
       HexagonSpiralIterator, spiral


# Various ways to index hexagons in a grid
# ----------------------------------------

abstract type Hexagon end

struct HexagonAxial <: Hexagon
    q::Int
    r::Int
end

struct HexagonCubic <: Hexagon
    x::Int
    y::Int
    z::Int
end

struct HexagonOffsetOddR <: Hexagon
    q::Int
    r::Int
end

struct HexagonOffsetEvenR <: Hexagon
    q::Int
    r::Int
end


# Basic constructors
# ------------------

hexagon(x::Int, y::Int, z::Int) = HexagonCubic(x, y, z)
hexagon(q::Int, r::Int) = HexagonAxial(q, r)


# Convert between hexagon indexing
# --------------------------------

function convert(::Type{HexagonAxial}, hex::HexagonCubic)
    HexagonAxial(hex.x, hex.z)
end

function convert(::Type{HexagonAxial}, hex::HexagonOffsetOddR)
    convert(HexagonAxial, convert(HexagonCubic, hex))
end

function convert(::Type{HexagonAxial}, hex::HexagonOffsetEvenR)
    convert(HexagonAxial, convert(HexagonCubic, hex))
end

function convert(::Type{HexagonCubic}, hex::HexagonAxial)
    HexagonCubic(hex.q, -hex.q - hex.r, hex.r)
end

function convert(::Type{HexagonCubic}, hex::HexagonOffsetOddR)
    x = hex.q - (hex.r >> 1)
    z = hex.r
    y = -x - z
    HexagonCubic(x, y, z)
end

function convert(::Type{HexagonCubic}, hex::HexagonOffsetEvenR)
    x = hex.q - (hex.r >> 1) + int(isodd(hex.r))
    z = hex.r
    y = -x - z
    HexagonCubic(x, y, z)
end

function convert(::Type{HexagonOffsetOddR}, hex::HexagonCubic)
    q = hex.x + (hex.z >> 1)
    r = hex.z
    HexagonOffsetOddR(q, r)
end

function convert(::Type{HexagonOffsetOddR}, hex::HexagonAxial)
    convert(HexagonOffsetOddR, convert(HexagonCubic, hex))
end

function convert(::Type{HexagonOffsetEvenR}, hex::HexagonCubic)
    q = hex.x + (hex.z >> 1) + int(isodd(hex.z))
    r = hex.z
    HexagonOffsetEvenR(q, r)
end

function convert(::Type{HexagonOffsetEvenR}, hex::HexagonAxial)
    convert(HexagonOffsetEvenR, convert(HexagonCubic, hex))
end


# Neighbor hexagon iterator
# -------------------------

struct HexagonNeighborIterator
    hex::HexagonCubic
end

const CUBIC_HEX_NEIGHBOR_OFFSETS = [
     1 -1  0;
     1  0 -1;
     0  1 -1;
    -1  1  0;
    -1  0  1;
     0 -1  1;
]

neighbors(hex::Hexagon) = HexagonNeighborIterator(convert(HexagonCubic, hex))

length(::HexagonNeighborIterator) = 6

function iterate(it::HexagonNeighborIterator, state=1)
    state > 6 && return nothing
    dx = CUBIC_HEX_NEIGHBOR_OFFSETS[state, 1]
    dy = CUBIC_HEX_NEIGHBOR_OFFSETS[state, 2]
    dz = CUBIC_HEX_NEIGHBOR_OFFSETS[state, 3]
    neighbor = HexagonCubic(it.hex.x + dx, it.hex.y + dy, it.hex.z + dz)
    return (neighbor, state + 1)
end


# Diagonal hexagon iterator
# -------------------------

struct HexagonDiagonalIterator
    hex::HexagonCubic
end

const CUBIC_HEX_DIAGONAL_OFFSETS = [
    +2 -1 -1;
    +1 +1 -2;
    -1 +2 -1;
    -2 +1 +1;
    -1 -1 +2;
    +1 -2 +1;
]

diagonals(hex::Hexagon) = HexagonDiagonalIterator(convert(HexagonCubic, hex))

length(::HexagonDiagonalIterator) = 6

function iterate(it::HexagonDiagonalIterator, state=1)
    state > 6 && return nothing
    dx = CUBIC_HEX_DIAGONAL_OFFSETS[state, 1]
    dy = CUBIC_HEX_DIAGONAL_OFFSETS[state, 2]
    dz = CUBIC_HEX_DIAGONAL_OFFSETS[state, 3]
    diagonal = HexagonCubic(it.hex.x + dx, it.hex.y + dy, it.hex.z + dz)
    return (diagonal, state + 1)
end


# Iterator over the vertices of a hexagon
# ---------------------------------------

struct HexagonVertexIterator
    x_center::Float64
    y_center::Float64
    xsize::Float64
    ysize::Float64

    function HexagonVertexIterator(x, y, xsize=1.0, ysize=1.0)
        new((Float64(x)), (Float64(y)),
            (Float64(xsize)), (Float64(ysize)))
    end

    function HexagonVertexIterator(hex::Hexagon,
                                   xsize=1.0, ysize=1.0, xoff=0.0, yoff=0.0)
        c = center(hex, xsize, ysize, xoff, yoff)
        new((Float64(c[1])), (Float64(c[2])),
            (Float64(xsize)), (Float64(ysize)))
    end
end

function vertices(hex::Hexagon, xsize=1.0, ysize=1.0, xoff=0.0, yoff=0.0)
    c = center(hex, xsize, ysize, xoff, yoff)
    HexagonVertexIterator(c[1], c[2], xsize, ysize)
end

# TODO: remove this function?
function hexpoints(x, y, xsize=1.0, ysize=1.0)
    collect(Tuple{Float64, Float64},
            HexagonVertexIterator(Float64(x), Float64(y),
                                  Float64(xsize), Float64(ysize)))
end

length(::HexagonVertexIterator) = 6

function iterate(it::HexagonVertexIterator, state=1)
    state > 6 && return nothing
    theta = 2*pi/6 * (state-1+0.5)
    x_i = it.x_center + it.xsize * cos(theta)
    y_i = it.y_center + it.ysize * sin(theta)
    return ((x_i, y_i), state + 1)
end

struct HexagonDistanceIterator
    hex::HexagonCubic
    n::Int
end

function hexagons_within(n::Int, hex::Hexagon = hexagon(0, 0))
    cubic_hex = convert(HexagonCubic, hex)
    HexagonDistanceIterator(hex, n)
end
hexagons_within(hex::Hexagon, n::Int) = hexagons_within(n, hex)

length(it::HexagonDistanceIterator) = it.n * (it.n + 1) * 3 + 1

function iterate(it::HexagonDistanceIterator, state=(-it.n, 0))
    x, y = state
    x > it.n && return nothing
    z = -x-y
    hex = HexagonCubic(x, y, z)
    y += 1
    if y > min(it.n, it.n-x)
        x += 1
        y = max(-it.n, -it.n - x)
    end
    hex, (x, y)
end


collect(it::HexagonDistanceIterator) = collect(HexagonCubic, it)

# Iterator over a ring of hexagons
# ---------------------------------------

struct HexagonRingIterator
    hex::HexagonCubic
    n::Int
end

function ring(n::Int, hex::Hexagon = hexagon(0, 0))
    # println("New hexring with center $hex and n $n")
    cubic_hex = convert(HexagonCubic, hex)
    HexagonRingIterator(cubic_hex, n)
end
ring(hex::Hexagon, n::Int) = ring(n, hex)

length(it::HexagonRingIterator) = it.n * 6

function iterate(it::HexagonRingIterator,
                 state::(Tuple{Int, HexagonCubic})=(1, neighbor(it.hex, 5, it.n)))
    hex_i, cur_hex = state
    hex_i > length(it) && return nothing
    # println("HexagonRingIterator: at position $hex_i ($cur_hex)")
    ring_part = div(hex_i - 1, it.n) + 1
    next_hex = neighbor(cur_hex, ring_part)
    cur_hex, (hex_i + 1, next_hex)
end

collect(it::HexagonRingIterator) = collect(HexagonCubic, it)

# Iterator over all hexes within a certain distance
# -------------------------------------------------

struct HexagonSpiralIterator
    hex::HexagonCubic
    n::Int
end

struct HexagonSpiralIteratorState
    hexring_i::Int
    hexring_it::HexagonRingIterator
    hexring_it_i::Int
    hexring_it_hex::HexagonCubic
end

function spiral(n::Int, hex::Hexagon = hexagon(0, 0))
    cubic_hex = convert(HexagonCubic, hex)
    HexagonSpiralIterator(cubic_hex, n)
end
spiral(hex::Hexagon, n::Int) = spiral(n, hex)

length(it::HexagonSpiralIterator) = it.n * (it.n + 1) * 3

# The state of a HexagonSpiralIterator consists of
# 1. an Int, the index of the current ring
# 2. a HexagonRingIterator and its state to keep track of the current position
#    in the ring.

function iterate(it::HexagonSpiralIterator)
    first_ring = ring(it.hex, 1)
    iterate(it, HexagonSpiralIteratorState(1, first_ring, start(first_ring)...))
end

function iterate(it::HexagonSpiralIterator, state::HexagonSpiralIteratorState)
    state.hexring_i > it.n && return nothing
    # Get current state
    hexring_i, hexring_it, hexring_it_i, hexring_it_hex =
        state.hexring_i, state.hexring_it, state.hexring_it_i, state.hexring_it_hex
    # Update state of inner iterator
    hexring_it_hex, (hexring_it_i, hexring_it_hex_next) =
                next(hexring_it, (hexring_it_i, hexring_it_hex))
    # Check if inner iterator is done, and update if necessary
    if done(hexring_it, (hexring_it_i, hexring_it_hex_next))
        hexring_i += 1
        hexring_it = ring(it.hex, hexring_i)
        hexring_it_i, hexring_it_hex_next = start(hexring_it)
        # println("In new ring $hexring_it")
    end

    # println("Currently at $hexring_it_hex, hexring is $hexring_it, state is $((hexring_i, (hexring_it_i, hexring_it_hex)))")
    hexring_it_hex, HexagonSpiralIteratorState(hexring_i, hexring_it,
                                               hexring_it_i, hexring_it_hex_next)
end

collect(it::HexagonSpiralIterator) = collect(HexagonCubic, it)

# Utilities
# ---------

function distance(a::Hexagon, b::Hexagon)
    hexa = convert(HexagonCubic, a)
    hexb = convert(HexagonCubic, b)
    max(abs(hexa.x - hexb.x),
        abs(hexa.y - hexb.y),
        abs(hexa.z - hexb.z))
end

function center(hex::Hexagon, xsize=1.0, ysize=1.0, xoff=0.0, yoff=0.0)
    axh = convert(HexagonAxial, hex)
    (xoff + xsize * sqrt(3) * (axh.q + axh.r/2), yoff + ysize * (3/2) * axh.r)
end

# TODO: Split up in two functions for performance (distance)?
function neighbor(hex::HexagonCubic, direction::Int, distance::Int = 1)
    dx = CUBIC_HEX_NEIGHBOR_OFFSETS[direction, 1] * distance
    dy = CUBIC_HEX_NEIGHBOR_OFFSETS[direction, 2] * distance
    dz = CUBIC_HEX_NEIGHBOR_OFFSETS[direction, 3] * distance
    HexagonCubic(hex.x + dx, hex.y + dy, hex.z + dz)
end

function cube_linedraw(a::Hexagon, b::Hexagon)
    hexa = convert(HexagonCubic, a)
    hexb = convert(HexagonCubic, b)
    N = distance(hexa, hexb)
    dx, dy, dz = hexb.x - hexa.x, hexb.y - hexa.y, hexb.z - hexa.z
    ax, ay, az = hexa.x + 1e-6, hexa.y + 1e-6, hexa.z - 2e-6
    map(i -> nearest_cubic_hexagon(ax + i*dx, ay + i*dy, az + i*dz), 0:(1/N):1)
end

# Find the nearest hexagon in cubic coordinates.
function nearest_cubic_hexagon(x::Real, y::Real, z::Real)
    rx, ry, rz = round(Integer, x), round(Integer, y), round(Integer, z)
    x_diff, y_diff, z_diff = abs(rx - x), abs(ry - y), abs(rz - z)

    if x_diff > y_diff && x_diff > z_diff
        rx = -ry - rz
    elseif y_diff > z_diff
        ry = -rx - rz
    else
        rz = -rx - ry
    end

    HexagonCubic(rx, ry, rz)
end

# Return the index (in cubic coordinates) of the hexagon containing the
# point x, y
function cube_round(x, y, xsize=1.0, ysize=1.0)
    x /= xsize
    y /= ysize
    q = sqrt(3)/3 * x - y/3
    r = 2 * y / 3
    h = nearest_cubic_hexagon(q, -q - r, r)
    return h
end

end # module Hexagons
