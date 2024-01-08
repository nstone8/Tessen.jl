module Tessen

using LinearAlgebra, Unitful, RecipesBase

"""
```julia
HatchLine(p,v)
```
A `HatchLine` is defined by a point `p` and a direction vector `v`.
Both are `Vector`s with length 2 (slices are always in x-y plane)
with units of `Unitful.Length`
"""
struct HatchLine
    #internally all lengths will be stored in microns
    p :: Vector{<:Number}
    v :: Vector{<:Number}
    function HatchLine(p::Vector{<:Unitful.Length},v::Vector{<:Unitful.Length})
        @assert length(p) == 2 "point used to define a HatchLine must have 2 coordinates"
        @assert length(v) == 2 "vector used to define a HatchLine must have 2 coordinates"
        #convert to microns and strip the units
        new(ustrip.(u"μm",p),ustrip.(u"μm",v))
    end
end

"""
```julia
pointalong(hl,t)
```
Get the point along a `HatchLine` corresponding to parametric coordinate `t`
"""
pointalong(hl::HatchLine,t::Number) = hl.p + (t*hl.v)

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
boundpoints(edge,[numpoints=100])
```
Get points along `edge` for plotting. If `edge` isn't straight, get `numpoints` points.
"""
function boundpoints end

"""
```julia
LineEdge(p1,p2)
```
A linear edge defined by two endpoints `p1` and `p2`, which should both be
supplied as `Vector`s of length 2 with units of `Unitful.Length`
"""
struct LineEdge <: Edge
    #all coordinates in microns internally
    p1 :: Vector{<:Number}
    p2 :: Vector{<:Number}
    function LineEdge(p1::Vector{<:Unitful.Length},p2::Vector{<:Unitful.Length})
        @assert all([p1,p2]) do p
            length(p) == 2
        end
        #convert to microns, strip units
        new(ustrip.(u"μm",p1),ustrip.(u"μm",p2))
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

boundpoints(le::LineEdge) = [le.p1,le.p2]

"""
```julia
ArcEdge(c,r,startangle,stopangle)
```
A constant radius arc edge defined by `c`, a vector with length 2 giving the center
of curvature, `r`, the radius of curvature as well as the angles at which the arc
starts and stops. These angles should be supplied in radians and be in the interval
``[0,2pi]``, a value of zero corresponds to a point with coordinates `c + [r,0]` and
a value of `pi/2` corresponds to a point with coordinates `c + [0,r]`.
"""
struct ArcEdge <: Edge
    #all length dimensions in microns internally
    c :: Vector{<:Number}
    r :: Number
    startangle :: Number
    stopangle :: Number
    function ArcEdge(c::Vector{<:Unitful.Length},r::Unitful.Length,startangle::Number,stopangle::Number)
        @assert length(c) == 2 "all coordinates must have 2 entries"
        @assert startangle != stopangle
        #convert to microns and strip units
        new(ustrip.(u"μm",c),ustrip(u"μm",r),startangle,stopangle)
    end
end

"""
```julia
zrotate(v,rotation)
```
rotate a `Vector` with length 2 (defined on xy plane) around the z axis
"""
function zrotate(v,rotation)
    @assert length(v) == 2
    rotmatrix=[cos(rotation) -sin(rotation)
               sin(rotation) cos(rotation)]
    return rotmatrix*v
end

function intersections(ae::ArcEdge,hl::HatchLine)
    #first rotate c around hl.p such that hl.v would be parallel to [1,0]
    rotamount = -atan(reverse(hl.v)...)
    #we're subtracting so we can treat hl.p as the origin
    crot = zrotate(ae.c-hl.p,rotamount)

    #abs(crot[2]) is now the closest distance between c and the hatch line. If this
    #value is greater than ae.r, there is no intersection with a full circle centered
    #on c. If this value is equal to r, there is one intersection, if it is less than
    #r there are two. We will calculate these intersections and then filter for points
    #on the arc

    #relrotcircintersections give the coordinates of potential intersections in a frame
    #with ae.c at its origin that has its x axis along hl.v
    relrotcircintersections = if abs(crot[2]) > ae.r
        #no intersections with full circle
        nothing
    elseif abs(crot[2]) == ae.r
        #one intersection at [crot[1],0] in the frame with hl.p at the origin, subtract
        #off crot to get what we want
        [[crot[1],0] - crot] #we want a vector of coordinates
    else
        #In our rotated frame, the x distance between ae.c and the intersections are
        xoffset = sqrt((ae.r^2) + (crot[2]^2))
        map([1,-1]) do s
            #the y coordinate of both intersections in the rotated frame with ae.c at
            #its center is -crot[2]
            [s*xoffset,-crot[2]]
        end
    end

    #if our hatchline doesn't intersect the circle we're done
    if isnothing(relrotcircintersections)
        return nothing
    end
    #now we need to test if intersections with the circle lie on our arc

    #find the acceptable interval in this rotated frame, take the mod so that all
    #angles are between 0 and 2pi
    rotangleinterval=([ae.startangle, ae.stopangle] + rotamount) .% 2pi
    filter!(relrotcircintersections) do rrci
        #angle of the intersection
        iangle = atan(reverse(rrci)...)
        #atan gives answers in the range [-pi,pi], we want a range of [0,2pi]
        iangle = (iangle < 2pi) ? iangle + 2pi : iangle
        #if rotangleinterval[2] > rotangleinterval[1] then iangle must lie between them to be on the arc
        if rotangleinterval[2] > rotangleinterval[1]
            return rotangleinterval[1] <= iangle <= rotangleinterval[2]
        else
            #iangle must be greater than rotangleinterval[1] or less than rotangleinterval[2]
            return !(rotangleinterval[2] <= iangle <= rotangleinterval[1])
        end
    end
    #if relrotcircintersections is empty we're done
    if isempty(relrotcircintersections)
        return nothing
    end

    #now we need to get the parametric coordinate t such that pointalong(hl,t) gives us each intersection
    return map(relrotcircintersections) do rrci
        #the parametric coordinate is the distance along hl.v (crot[1] + rrci[1]) divided by the magnitude of hl.v
        (crot[1] + rrci[1]) ./ sqrt(sum(hl.v .^ 2))
    end
end

function boundpoints(ae::ArcEdge,numpoints=100)
    angleinterval = [ae.startangle, ae.stopangle]
    if angleinterval[2] < angleinterval[1]
        angleinterval[2] += 2pi
    end
    map(range(angleinterval...,numpoints)) do a
        ae.c + ae.r*[cos(a),sin(a)]
    end
end

@recipe function plotedge(e::Edge)
    points=boundpoints(e)
    pointsmat = vcat(permutedims.(points)...)
    :legend --> false
    :aspect_ratio --> 1
    (pointsmat[:,1],pointsmat[:,2])
end

#Contours are a list of edges
struct Contour
    edges :: Vector{<:Edge}
end

#slices consist of a list of contours
struct Slice
    contours :: Vector{Contour}
end

#blocks are a list of slices with corresponding z coordinates

#jobs are a list of blocks with corresponding [x,y,z] offsets

#=====Hatching function should probably operate on Blocks======================
function hatch(s::Slice,hatchdistance::Number,hatchdirection::Vector{<:Number})
    #fill me out
end
==============================================================================#

end # module Tessen
