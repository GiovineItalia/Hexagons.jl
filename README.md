
# Hexagons

[![Build
Status](https://travis-ci.org/GiovineItalia/Hexagons.jl.svg?branch=master)](https://travis-ci.org/GiovineItalia/Hexagons.jl)

This package provides some basic utilities for working with hexagonal grids. It
is largely works from Amit Patel's [terrific
refererence](http://www.redblobgames.com/grids/hexagons/).

## Synopsis

Hexagonal grids can be indexed a number of different ways. Indexes are
represented with one of the Hexagon types. The following are currently provided:

```julia
HexagonAxial(q, r)
HexagonCubic(x, y, z)
HexagonOffsetOddR(q, r)
HexagonOffsetEvenR(q, r)
```

One indexing system can be converted to another with `convert`.

```julia
convert(HexagonOffsetOddR, HexagonAxial(2, 4))
```

The six points (in cartesian space) of a hexagon can be iterated through with
`points`.

```julia
for (x, y) in vertices(HexagonAxial(2, 3))
    # ...
end
```

The center point can be obtained with `center`

```julia
x, y = center(HexagonAxial(2, 3))
```

A point in cartesian space can be mapped to the index of the hexagon that
contains it with the `cube_round` function.

```julia
h = cube_round(23.5, 4.67)
```

## Status

This library is not mature or complete, but provides just enough to implement
hexagonal bin visualizations. If your hexagon project requires something
that's not provided, file bug or pull request.


