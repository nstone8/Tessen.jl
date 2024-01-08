module Tessen

"""
```julia
HatchLine(p,v)
```
A `HatchLine` is defined by a point `p` and a direction vector `v`.
Both are `Vector`s with length 2 (slices are always in x-y plane)
"""
struct HatchLine
    p :: Vector{<:Number}
    v :: Vector{<:Number}
    function HatchLine(p,v)
        @assert length(p) == 2 "point used to define a HatchLine must have 2 coordinates"
        @assert length(v) == 2 "vector used to define a HatchLine must have 2 coordinates"
        new(p,v)
    end
end

"""
```julia
pointalong(hl,t)
```
Get the point along a `HatchLine` corresponding to parametric coordinate `t`
"""
pointalong(hl::HatchLine,t::Number) = hl.p + (t*v)

#supertype for Edges
abstract type Edge end

"""
```julia
intersections(edge,hatchline)
```
Get all values of `t` such that `pointalong(hatchline,t)` lies on `edge`.
Returns `nothing` if there is no intersection.
"""
function intersections end

"""
```julia
LineEdge(p1,p2)
```
A linear edge defined by two endpoints `p1` and `p2`, which should both be
supplied as `Vector`s of length 2.
"""
struct LineEdge <: Edge
    p1 :: Vector{<:Number}
    p2 :: Vector{<:Number}
    function LineEdge(p1,p2)
        @assert all([p1,p2]) do p
            length(p) == 2
        end
        new(p1,p2)
    end
end

function intersections(le::LineEdge,hl::HatchLine)
    #I worked this out on paper
    p1 = le.p1 #point on first line
    v1 = le.p2 - le.p1 #direction of first line
    p2 = hl.p #point on second line
    v2 = hl.v #direction of second line
    
    A = hcat(v1-v2, -v2)
    if iszero(det(A))
        return nothing
    end
    b = p2 - p1
    (t,lambda) = A\b
    #we now have t and lambda such that (p1 + v1*t) == pointalong(hl,t+lambda)
    #this intersection point lies between our endpoints if 0 <= t <= 1
    return (0 <= t <= 1) ? [t+lambda] : nothing
end

#slices consist of a list of edges
struct Slice
    edges :: Vector{<:Edge}
end

function hatch(s::Slice,hatchdistance::Number,hatchdirection::Vector{<:Number})
    #fill me out
end

end # module Tessen
