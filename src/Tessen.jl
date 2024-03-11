module Tessen

using LinearAlgebra, Unitful, RecipesBase, Statistics

export LineEdge, ArcEdge, Contour, Slice, translate, rotate
export Block, SuperBlock, blocks, slices, hatch, HatchedSlice
export origin, rotation

#absolute tolerance for isapprox
atol = 1E-12

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
intersections(obj,hatchline)
```
Get all values of `t` such that `pointalong(hatchline,t)` lies on `obj`.
Returns `nothing` if there is no intersection. This function is defined
for `Edge`, `Contour` and `Slice`.
"""
function intersections end

"""
```julia
boundbox(obj)
```
Get a box (as a length 2 vector of corner coordinates) guarenteed to enclose a
`Tessen` object. The provided box may be slightly larger than the `element` due
to my laziness. Always returns a numeric type, which should be interpreted as having
units of microns.
"""
function boundbox end

"""
```julia
boundpoints(edge,[numpoints=100])
```
Get points along `edge` for plotting. If `edge` isn't straight, get `numpoints` points.
"""
function boundpoints end

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

"""
```julia
rotate(obj,amount,[point])
```
Rotate a `Tessen` object by `amount` radians around `point` in the xy plane.
`point` should be a `Vector` of `Unitful.Length`. If `point` is not provided,
rotate around the origin.
"""
function rotate end

#define the two argument version here
rotate(obj,amount) = rotate(obj,amount,[0,0] * u"μm")

"""
```julia
translate(obj,displacement)
```
Translate a `Tessen` object in the xy plane. `displacement` should be a `Vector`
of `Unitful.Length`.
"""
function translate end

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

function boundbox(le::LineEdge)
    xcoords = [le.p1[1], le.p2[1]]
    ycoords = [le.p1[2], le.p2[2]]
    [[minimum(xcoords),minimum(ycoords)], [maximum(xcoords),maximum(ycoords)]]
end

function rotate(le::LineEdge,amount,pointunits::Vector{<:Unitful.Length})
    #get the point coordinates in um
    point = ustrip.(u"μm",pointunits)
    originalendpoints = [le.p1,le.p2]
    #first translate everything so point lies at the origin
    transpoints = [ep - point for ep in originalendpoints]
    #do the rotation
    rotpoints = zrotate.(transpoints,amount)
    #translate back, add back on units
    finalpoints = [rp + point for rp in rotpoints] .* u"μm"
    LineEdge(finalpoints...)
end

function translate(le::LineEdge,displacement::Vector{<:Unitful.Length})
    #need to add units to the current points (all stored in microns)
    oldpoints = [le.p1, le.p2] .* u"μm"
    LineEdge((op + displacement for op in oldpoints)...)
end

"""
```julia
approxeqle(x,y)
```
Test if x is less than or approximately equal to y
"""
function approxeqle(x,y)
    #need to supply atol because we want to compare values near zero
    #I believe this value is insanely conservative
    (x<y) || isapprox(x,y;atol)
end

function intersections(le::LineEdge,hl::HatchLine)
    #I worked this out on paper
    p1 = le.p1 #point on first line
    v1 = le.p2 - le.p1 #direction of first line
    p2 = hl.p #point on second line
    v2 = hl.v #direction of second line

    #first check to see if le is colinear to hl, if so we should just return our endpoints
    #checking the slope by division doesn't work well if the line is vertical due to numerical
    #precision issues, we will call things parallel if either the slope or inverse slope is the
    #same
    if (isapprox(/(v1...),/(v2...);atol) || isapprox(/(reverse(v1)...),/(reverse(v2)...);atol))
        #the lines are parallel, these are colinear if either p1==p2
        samepoints = all(isapprox.(p1,p2;atol))
        #or if (p1-p2) is parallel to v2
        #same numerical precision thing here
        pvec=p1-p2
        parallel = (isapprox(/(pvec...),/(v2...);atol) ||
            isapprox(/(reverse(pvec)...),/(reverse(v2)...);atol))
        if samepoints || parallel
            #need to convert these endpoints into a parametric coordinate
            return map([le.p1,le.p2]) do ep
                #get a vector going from hl.p to ep
                d = ep - hl.p
                #we know everything is colinear already, so we only need to look at one coordinate
                #to get our parametric distance. We'll use the one that is largest to minimize
                #numerical issues
                (m,im) = findmax(d)
                m/hl.v[im]
            end
        end
    end
    
    A = hcat(v1-v2, -v2)
    if iszero(det(A))
        return nothing
    end
    b = p2 - p1
    (t,lambda) = A\b
    #we now have t and lambda such that (p1 + v1*t) == pointalong(hl,t+lambda)
    #this intersection point lies between our endpoints if 0 <= t <= 1
    #we will use approxeqle (defined above) to handle numerical precision issues
    return (approxeqle(0,t) && approxeqle(t,1)) ? [t+lambda] : nothing
end

boundpoints(le::LineEdge) = [le.p1,le.p2]

"""
```julia
ArcEdge(c,r,startangle,stopangle)
```
A constant radius arc edge defined by `c`, a vector with length 2 giving the center
of curvature, `r`, the radius of curvature as well as the angles at which the arc
starts and stops (the arc will be traced between these point in a counterclockwise
direction). These angles should be supplied in radians.
A value of zero corresponds to a point with coordinates `c + [r,0]` and
a value of `pi/2` corresponds to a point with coordinates `c + [0,r]`.
"""
struct ArcEdge <: Edge
    #all length dimensions in microns internally
    c :: Vector{<:Number}
    r :: Number
    #these will be stored internally in the interval [0,2pi]
    startangle :: Number
    stopangle :: Number
    function ArcEdge(c::Vector{<:Unitful.Length},r::Unitful.Length,startangle::Number,stopangle::Number)
        @assert length(c) == 2 "all coordinates must have 2 entries"
        @assert startangle != stopangle
        #make sure all angles are between 0 and 2pi
        angles = map([startangle,stopangle]) do a
            if abs(a) > 2pi
                a %= 2pi #make sure abs(a) is less than 2pi
            end
            #make sure a is positive
            (a >= 0) ? a : a + 2pi
        end
        #convert to microns and strip units
        new(ustrip.(u"μm",c),ustrip(u"μm",r),angles...)
    end
end

#I'm going to be lazy and just give the bounding box of the whole circle
boundbox(ae::ArcEdge) = map([-1,1]) do corner
    ae.c + ae.r*repeat([corner],2)
end

function rotate(ae::ArcEdge,amount,pointunits::Vector{<:Unitful.Length})
    point = ustrip.(u"μm",pointunits) #convert to microns and strip units
    #translate c such that `point` lies at the origin
    transc = ae.c - point
    #rotate around point
    rotc = zrotate(transc,amount)
    #translate back, add units back on
    newc = (rotc + point) * u"μm"
    #add rotation to startangle and stopangle
    newangles = [ae.startangle,ae.stopangle] .+ amount
    ArcEdge(newc,ae.r*u"μm",newangles...)
end

function translate(ae::ArcEdge,displacement::Vector{<:Unitful.Length})
    #all we have to do is translate our center point
    oldc = ae.c * u"μm" #need to add back on units
    newc = oldc + displacement
    ArcEdge(newc,ae.r*1u"µm",ae.startangle,ae.stopangle)
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
    elseif isapprox(abs(crot[2]),ae.r;atol)
        #one intersection at [crot[1],0] in the frame with hl.p at the origin, subtract
        #off crot to get in untransformed coordinates.
        #[[crot[1],0] - crot] #we want a vector of coordinates

        #We should return nothing. A tangent point doesn't affect hatching
        nothing
    else
        #In our rotated frame, the x distance between ae.c and the intersections are
        xoffset = sqrt((ae.r^2) - (crot[2]^2))
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
    #make sure we're in the range [-2pi, 2pi]
    rotangleinterval=([ae.startangle, ae.stopangle] .+ rotamount) .% 2pi
    #make sure we're positive
    rotangleinterval = map(rotangleinterval) do ri
        ri >= 0 ? ri : ri+2pi
    end
    filter!(relrotcircintersections) do rrci
        #angle of the intersection
        iangle = atan(reverse(rrci)...)
        #atan gives answers in the range [-pi,pi], we want a range of [0,2pi]
        iangle = (iangle >= 0) ? iangle : iangle + 2pi
        #if rotangleinterval[2] > rotangleinterval[1] then iangle must lie between them to be on the arc
        if rotangleinterval[2] > rotangleinterval[1]
            return rotangleinterval[1] <= iangle <= rotangleinterval[2]
        else
            #iangle must be greater than rotangleinterval[1] or less than rotangleinterval[2]
            return !(rotangleinterval[2] < iangle < rotangleinterval[1])
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
    legend --> false
    aspect_ratio --> 1
    (pointsmat[:,1],pointsmat[:,2])
end

#Contours are a list of edges
struct Contour
    edges :: Vector{<:Edge}
end

#abstract supertype for Slices
abstract type AbstractSlice end

#slices consist of a list of contours
"""
```julia
Slice(contours)
```
A `Slice` represents a set of coplanar contours
"""
struct Slice <: AbstractSlice
    contours :: Vector{Contour}
end


"""
```julia
children(obj)
```
Get the children of a `Contour` or `Slice`
"""
function children end

children(c::Contour) = c.edges

#rotating and translating compound objects is just rotating/translating all the children
function rotate(c::CO,args...) where CO <: Union{Contour,Slice}
    CO([rotate(ch,args...) for ch in children(c)])
end

function translate(c::CO,args...) where CO <: Union{Contour,Slice}
    CO([translate(ch,args...) for ch in children(c)])
end


children(s::Slice) = s.contours
function intersections(x::Union{Contour,Slice},hl::HatchLine)
    childintersections = [intersections(cx,hl) for cx in children(x)]
    #get rid of nothing entries and vcat into a flat Vector
    inters = vcat(filter(childintersections) do ci
                      !isnothing(ci)
                  end...)
    if isempty(inters)
        return nothing
    end
    #try to get only unique values (we can intersect multiple entities at points
    #like vertices
    uniqueinters = Vector{eltype(inters)}()
    for thisinter in inters
        if isempty(uniqueinters)
            push!(uniqueinters,thisinter)
        elseif !any(isapprox.(thisinter,uniqueinters;atol))
            #we don't yet have something like thisinter
            push!(uniqueinters,thisinter)
        end
    end
    return uniqueinters
end

function boundbox(co::Union{Contour,Slice})
    childboxes = boundbox.(children(co))
    map(zip([minimum,maximum],[1,2])) do (m,corner)
        map([1,2]) do dim
            m(b[corner][dim] for b in childboxes)
        end
    end
end

@recipe function plotco(co::Union{Contour,Slice})
    legend --> false
    aspect_ratio --> 1
    for c in children(co)
        @series begin
            c
        end
    end
    return nothing
end

#HatchedSlice structs will have a series of line segments, guarenteed to be in a rational zigzag
#order
"""
```julia
HatchedSlice(hatchlines)
```
`struct` representing a hatched `Slice`. The `hatchlines` field of this struct
is a vector of point vectors describing line segments presented in order to
neatly hatch the `Slice` in a zigzag pattern. In between top level `Vectors` (i.e. between
`hatchlines[1]` and `hatchlines[2]` the tool path must be broken to ensure that the shape is
represented faithfully. This constructor is not intended to be used directly. Use the `hatch`
function instead.
"""
struct HatchedSlice <: AbstractSlice
    hatchlines :: Vector{Vector{Vector{<:Number}}}
end

function rotate(hs::HatchedSlice,amount::Number,pointunits::Vector{<:Unitful.Length})
    point = ustrip.(u"µm",pointunits)
    map(hs.hatchlines) do line
        map(line) do p
            transp = p - point
            zrotate(transp,amount) + point            
        end
    end |> HatchedSlice
end

function translate(hs::HatchedSlice,displacement::Vector{<:Unitful.Length})
    disp = ustrip.(u"µm",displacement)
    map(hs.hatchlines) do line
        map(line) do p
            p + disp
        end
    end |> HatchedSlice
end

function hatch(s::Slice,dhatchunits::Unitful.Length,hatchdir::Number)::HatchedSlice
    #=============================Hatching strategy=========================
    1) get the bounding box of our slice
    2) trace this box with a contour
    3) place the first hatchline on the center of our bounding box
    4) place hatchlines in each direction until we place a hatchline that doesn't intersect our box
    5) done
    ======================================================================#
    bboxcorners = boundbox(s) .* 1u"µm" #need units for the LineEdge constructor
    bboxedges = [LineEdge(bboxcorners[1],
                          [bboxcorners[2][1],bboxcorners[1][2]]),
                 LineEdge([bboxcorners[2][1],bboxcorners[1][2]],
                          bboxcorners[2]),
                 LineEdge(bboxcorners[2],
                          [bboxcorners[1][1],bboxcorners[2][2]]),
                 LineEdge([bboxcorners[1][1],bboxcorners[2][2]],
                          bboxcorners[1])]
    #hatchlines will be written such that they are perpendicular to hatchdir
    vhatch = dhatchunits*[cos(hatchdir+pi/2),sin(hatchdir+pi/2)]

    startpoint = mean(bboxcorners) #center of bounding box
    #get the distance between startpoints for neighboring hatchlines. these points will be
    #dhatch apart along hatchdir
    poffset = dhatchunits*[cos(hatchdir),sin(hatchdir)]

    #lay down our hatchlines. We will collect each direction from the center in different
    #vectors so we can rearrange into a continuous order
    hlinesvec = map([false,true]) do rev
        curpoint = startpoint
        #if we've already gone 'forward' reverse direction and don't redo start point
        if rev
            poffset *= -1
            curpoint += poffset
        end
        #lay down hatchlines in this direction
        thesehlines = HatchLine[]
        while true
            hl = HatchLine(curpoint,vhatch)
            #if this doesn't intersect our bounding box, we're done
            inters = [intersections(be,hl) for be in bboxedges]
            if all(isnothing,inters)
                break
            end
            #otherwise pop it into thesehlines and move curpoint
            push!(thesehlines,hl)
            curpoint += poffset
        end
        #need to explicitly return thesehlines so they end up in hlinesvec
        return thesehlines
    end
    #make one list of hatch lines in the correct order
    hlines = vcat(reverse(hlinesvec[2]),hlinesvec[1])
    #now build an array of point pairs alternating directions
    #build a generator that gives alternating true false to do the zigzag
    alternator = (isodd(i) for i in 1:length(hlines))
    #this will be a vector of vectors ordered such that vcatting and reshaping
    #can be used to get the point pairs we want
    intervec = map(zip(hlines,alternator)) do (hl,rev)
        #get the even number of points where `hl` intersects `s`
        inters = intersections(s,hl)
        #if there are no intersections, just add nothing
        if isnothing(inters)
            return nothing
        end
        #if inters has only one entry we are perfectly clipping the corner of a polygon
        #these cases (like tangent points on arcs) don't matter for hatching
        @assert length(inters) != 1 "guess i was wrong"
        if length(inters) == 1
            return nothing
        end
        @assert iseven(length(inters))
        #we will sort these into an order based on `rev`
        sort!(inters;rev)
        #now need to turn this into points rather than parametric coords
        [pointalong(hl,i) for i in inters]
    end
    #remove nothing entries (using the non-mutating filter so the type can change)
    filter!(intervec) do iv
        !isnothing(iv)
    end
    #====================Old version======================================
    #split the intersections into pairs via reshape        
    intermat = reshape(vcat(intervec...),2,:)
    #change the matrix into a vector of vectors of vectors
    hatchlines = map(1:size(intermat)[2]) do i
        intermat[:,i]
    end
    #now filter out any zero-length hatchlines
    filter!(hatchlines) do (p1,p2)
        !all(isapprox.(p1,p2;atol))
    end
    ====================================================================#
    #intervec now contains a vector of vectors of vectors. The innermost vectors are coordinates
    #in 2D, grouped into vectors of points which occur on the same hatch line. We only need to
    #pick up our pen (so to speak) to avoid filling in a desired hole if a second-level vector
    #has length > 2 (as this would imply intersection with a hatch line more than twice. Our goal
    #here will be to make a new vector of vectors of vectors where the second level vector is a
    #list of points that can be written continuously without lifting our pen
    numtype = eltype(intervec[1][1][1])
    hatchlines = Vector{Vector{Vector{numtype}}}()
    #the current path we are building
    curpath = Vector{Vector{numtype}}()
    for iv in intervec
        if length(iv) == 2
            #throw them on
            push!(curpath,iv...)
        else
            #put the first two points on
            push!(curpath,iv[1:2]...)
            for j in 2:Int(length(iv)/2)
                #each extra pair of points should be made into its own path, except for the last
                #two, which can be part of a new one
                push!(hatchlines,curpath)
                curpath = collect(iv[(2j-1):2j])
            end
        end
    end
    push!(hatchlines,curpath)
    HatchedSlice(hatchlines)
end

#ignore slices which are already hatched
hatch(hs::HatchedSlice,args...) = hs

"""
```julia
hatch(slice; dhatch [,hatchdir])
```
Hatch an `AbstractSlice` object with the provided hatch distance (as `Unitful.Length`)
and optional hatch direction (default is `0` radians). The hatch lines themselves will
be perpendicular to `hatchdir`.
"""
hatch(s::AbstractSlice;dhatch,hatchdir=0) = hatch(s,dhatch,hatchdir)

@recipe function ploths(hs::HatchedSlice)
    legend --> false
    aspect_ratio --> 1
    #intersperse hs.hatchlines with nothing pairs to break up the lines
    broken = [vcat(l,[[missing,missing]]) for l in hs.hatchlines]
    mat = hcat(vcat(broken...)...)
    (mat[1,:],mat[2,:])
end

"""
abstract type for Block and SuperBlock
we will assume there are methods for origin and rotation
as well as an inner constructor T(t::T,translation,rotation)
"""
abstract type LocalFrame{T} end

"""
```julia
slicetype(lf)
```
Get the slice type of a `LocalFrame`. (get `T` for `lf::LocalFrame{T}`)
"""
function slicetype(::LocalFrame{T}) where {T}
    T
end


"""
```julia
changeslicetype(T,lf)
```
Change the slice type of a `LocalFrame` to `T`
"""
function changeslicetype end

"""
```julia
translate(block, displacement[; preserveframe=false])
```
Translate a `Block` or `SuperBlock`. If preserveframe=true is passed, the local coordinate system
of `block` is not modified (this is accomplished by recursively translating all contained `Slice`
objects)
"""
function translate(lf::LocalFrame,displacement;preserveframe=false)
    if !preserveframe
        npftranslate(lf,displacement)
    else
        pftranslate(lf,displacement)
    end
end

"""
```julia
rotate(block, amount, [point; preserveframe=false])
```
Rotate a `Block` or `SuperBlock`. If preserveframe=true is passed, the local coordinate system
of `block` is not modified (this is accomplished by recursively moving all contained `Slice`
objects)
"""
function rotate(lf::LocalFrame,amount::Number,point::Vector{<:Unitful.Length};preserveframe=false)
    if !preserveframe
        npfrotate(lf,amount,point)
    else
        pfrotate(lf,amount,point)
    end
end

#make point optional
rotate(lf::LocalFrame,amount;preserveframe=false) = rotate(lf,amount,[0u"µm",0u"µm"];preserveframe)

"""
```julia
pftranslate(lf,displacement)
```
Translate a `LocalFrame` without modifying any local coordinate systems. This function
is called when `translate` is used on a `LocalFrame` with `preserveframe=true`
"""
function pftranslate end

"""
```julia
pfrotate(lf,amount,point)
```
Rotate a `LocalFrame` without modifying any local coordinate systems. This function
is called when `rotate` is used on a `LocalFrame` with `preserveframe=true`
"""
function pfrotate end

"""
```julia
npftranslate(lf,displacement)
```
Translate a `LocalFrame` by modifying the local coordinate system. This function
is called when `translate` is used on a `LocalFrame` with `preserveframe=false`
"""
function npftranslate(lf::LF,displacement::Vector{<:Unitful.Length}) where {LF <: LocalFrame}
    LF(lf,displacement,0)
end

"""
```julia
npfrotate(lf,amount,point)
```
Rotate a `LocalFrame` by modifying the local coordinate system. This function
is called when `rotate` is used on a `LocalFrame` with `preserveframe=false`.
"""
function npfrotate(lf::LF,amount,point::Vector{<:Unitful.Length}) where {LF <: LocalFrame}
    #point = ustrip.(u"µm",pointunits)
    #get the translation associated with this rotation
    o = origin(lf)[1:2] #xy coordinates
    #translate o such that point lies at the origin
    transo = o - point
    #rotate about o and translate back
    newo = zrotate(transo,amount) + point
    translation = newo - o
    #do the translation and rotation about the local origin
    LF(lf,translation,amount)
end

"""
```julia
Block(z1 => slice1, z2 => slice2...; origin, rotation)
```
Assemble a series of `Slice`s into a `Block`. `Slice` objects should be
provided as a series of `z => slice` `Pair`s where `z` is the elevation of
the slice with units of `Unitful.Length`. The optional `origin` and `rotation`
keywords define a local reference frame in which all the component `Slice`
geometry is interpreted.
"""
struct Block{T} <: LocalFrame{T}
    #origin position (in 3d) in microns
    origin::Vector{<:Number}
    #local frame rotation in radians
    rotation::Number
    #a dict where the keys are slice elevation and the values are a
    #vector of all the slices to be written at that elevation
    slices::Dict{<:Number,Vector{T}}
    function Block(origin::Vector{<:Unitful.Length},rotation::Number,
                   slicepairs::Pair...)
        #origin needs to be 3D
        @assert length(origin) == 3 "LocalFrame coordinate origins must have 3 coordinates"
        #strip units, we will represent everything internally in microns
        rawslice = [ustrip(u"µm",sp.first) => sp.second for sp in slicepairs]
        #change our slice pairs into a dict
        allz = Set([rs.first for rs in rawslice])
        slicetype = eltype([rs.second for rs in rawslice])
        @assert slicetype <: AbstractSlice "Blocks are built from slices"
        #initialize all required slice vectors in the dict
        slices = Dict(z => Vector{slicetype}() for z in allz)
        #now push all the slices onto the correct vector
        for (z,s) in rawslice
            push!(slices[z],s)
        end
        new{slicetype}(ustrip.(u"µm",origin),rotation,slices)
    end
    #LocalFrame innner constructor to make translation/rotation easier
    function Block{T}(b::Block{T},translation::Vector{<:Unitful.Length},rotation) where {T}
        #if translation isn't in 3D, go ahead and assume we mean movement in xy
        if length(translation) == 2
            #create a new vector so we don't modify `translation` in the calling scope
            translation=vcat(translation,0u"µm")
        end
        new{T}(ustrip.(u"µm",(b.origin * 1u"µm") + translation),
               b.rotation + rotation, b.slices)
    end
    #inner constructor to allow for conversions of Blocks{T} to Block{AbstractSlice}
    function Block{AbstractSlice}(b::Block)
        numtype=promote_type(typeof.(keys(b.slices))...)
        new{AbstractSlice}(b.origin,b.rotation,
                           convert(Dict{numtype,Vector{AbstractSlice}},b.slices))
    end
end

function Base.convert(::Type{Block{AbstractSlice}},b::Block)
    Block{AbstractSlice}(b)
end

function changeslicetype(T::Type{<:AbstractSlice},b::Block)
    convert(Block{T},b)
end

#make kwargs optional
Block(slices...;origin=[0u"µm",0u"µm",0u"µm"],rotation=0) = Block(origin,rotation,slices...)
#required localframe methods
origin(b::Block) = b.origin * 1u"µm"
rotation(b::Block) = b.rotation

"""
```julia
slices(block)
```
Get all the `Slice`s and `HatchedSlice`s which make up a `Block`. Returned
as a `Vector` of `z => slice` pairs (i.e. the format expected by the `Block`
constructor).
"""
slices(b::Block) = [(z*u"µm") => slice for (z,slicevec) in b.slices for slice in slicevec]
function pftranslate(b::Block,displacement::Vector{<:Unitful.Length})
    #if displacement has length 2, assume translation in xy
    if length(displacement) == 2
        displacement = vcat(displacement,0u"µm")
    end
    @assert length(displacement) == 3
    #translate every slice individually in xy
    #need to convert xy displacement into the local coordinate system
    localdisp = zrotate(displacement[1:2],-rotation(b))
    flatslices = [(z + displacement[3]) => translate(ts,localdisp) for (z,ts) in slices(b)]
    Block(flatslices...,origin=origin(b),rotation=rotation(b))
end

"""
```julia
merge(blocks...)
```
Create a `Block` which contains all of the slices present in `blocks`. The source `Blocks` must
all share a common local coordinate system (i.e. `origin(b)` and `rotation(b)` must be the same
for all arguments
"""
function Base.merge(b1::Block,blks::Block...)
    #if we're only given one block just bounce it back
    if length(blks) == 0
        return b1
    end
    @assert all(origin(bi) == origin(b1) for bi in blks) &&
        all(rotation(bi) == rotation(b1) for bi in blks) "all arguments must share a coordinate system"
    #build a vector of vectors of all the slices
    allslices = [slices(b) for b in vcat(b1,blks...)]
    #build a Block containing all of them
    Block(vcat(allslices...)...,origin=origin(b1),rotation=rotation(b1))        
end

function pfrotate(b::Block,amount::Number,point::Vector{<:Unitful.Length})
    #figure out the coordinates of `point` in our local coordinate system
    pvec = point - origin(b)[1:2] #vector pointing from the origin of b to point in our global frames xy plane
    #minus sign on rotation because rotation(b) is the rotation required to convert from the
    #local frame to the global frame, we're going in the opposite direction
    plocal = zrotate(pvec,-rotation(b))
    #now we can just do the rotations about plocal in the local frame
    rotslices = [z => rotate(s,amount,plocal) for (z,s) in slices(b)]
    Block(rotslices...,origin=origin(b),rotation=rotation(b))
end

function hatch(b::Block, dhatch::Unitful.Length, bottomdir::Number,diroffset::Number)::Block{HatchedSlice}
    #get the elevation of every slice in order
    oldslices = slices(b)
    sortedz = Set(z for (z,_) in oldslices) |> collect |> sort
    #get a vector of the same length giving the hatch direction of each slice
    dirvec = range(start=bottomdir, step=diroffset, length=length(sortedz))
    #make a dict so we can grab direction by z
    dirdict = Dict(z => d for (z,d) in zip(sortedz,dirvec))
    #build a vector of z => hatchedslice pairs
    pairvec = map(oldslices) do (z,s)
        z => hatch(s;dhatch,hatchdir=dirdict[z])
    end
    Block(pairvec..., origin=origin(b), rotation=rotation(b))
end

"""
```julia
hatch(block; dhatch [,bottomdir, diroffset])
```
Hatch a `Block` with uniform hatching distance `dhatch`. Any `Slice`s in the block which are
already hatched will not be modified. `dhatch` is the uniform slicing distance, the optional
keyword arguments `bottomdir` and `diroffset` control the hatching direction of the bottommost
slice and the direction offset between slices.
"""
function hatch(b::Block; dhatch, bottomdir=0, diroffset=pi/2)
    hatch(b,dhatch,bottomdir,diroffset)
end

"""
```julia
SuperBlock(blocks...;origin,rotation)
```
Combine `Block`s and `SuperBlock`s into a single entity. The optional `origin` and
`rotation` keyword arguments define a local coordinate system in which the enclosed
geometry is interpreted. The `blocks` argument is treated such that order matters,
earlier arguments will be printed first.
"""
struct SuperBlock{T} <: LocalFrame{T}
    #origin position (in 3d) in microns
    origin::Vector{<:Number}
    #local frame rotation in radians
    rotation::Number
    blocks::Vector{LocalFrame{<:T}}

    function SuperBlock{T}(origin::Vector{<:Unitful.Length},rotation::Number,
                        blocks::Vector{<:LocalFrame{<:T}}) where {T}
        @assert length(origin) == 3 "LocalFrame coordinate origins must have 3 coordinates"
        new{T}(ustrip.(u"µm",origin),rotation,collect(blocks))
    end

    #localframe constructor
    function SuperBlock{T}(sb::SuperBlock{T},translation::Vector{<:Unitful.Length},rotation) where {T}
        
        #if translation isn't in 3D, go ahead and assume we mean movement in xy
        if length(translation) == 2
            translation = vcat(translation,0u"µm")
        end
        
        new{T}(ustrip.(u"µm",(sb.origin * 1u"µm") + translation),
               sb.rotation + rotation, sb.blocks)
    end
end

function changeslicetype(T::Type{<:AbstractSlice},sb::SuperBlock)
    b = sb.blocks
    convb = map(b) do bi
        changeslicetype(T,bi)
    end
    #need to add the units back on
    SuperBlock{T}(sb.origin*u"µm",sb.rotation,convb)
end

function SuperBlock(blocks::LocalFrame...;origin=[0u"µm",0u"µm",0u"µm"],rotation=0)
    blockvec=collect(blocks)
    T=promote_type(slicetype.(blockvec)...)
    convvec = map(blockvec) do bv
        changeslicetype(T,bv)
    end
    SuperBlock{T}(origin,rotation,convvec)
end
origin(sb::SuperBlock) = sb.origin * 1u"µm"
rotation(sb::SuperBlock) = sb.rotation

"""
```julia
blocks(sb)
```
Get all of the blocks which make up a `SuperBlock`
"""
blocks(sb::SuperBlock) = sb.blocks

function pftranslate(sb::SuperBlock,displacement::Vector{<:Unitful.Length})
    #if displacement has length 2, assume translation in xy
    if length(displacement) == 2
        displacement = vcat(displacement,0u"µm")
    end
    #need to convert xy displacement into local coordinate system
    xylocaldisp = zrotate(displacement[1:2],-rotation(sb))
    localdisp = vcat(xylocaldisp,displacement[3])
    translatedblocks = map(blocks(sb)) do b
        translate(b,localdisp,preserveframe=true)
    end
    SuperBlock(translatedblocks...,origin=origin(sb),rotation=rotation(sb))
end

function pfrotate(sb::SuperBlock,amount::Number,point::Vector{<:Unitful.Length})
    #just have to apply the rotation to every block
    newblocks = map(blocks(sb)) do b
        rotate(b,amount,point,preserveframe=true)
    end
    SuperBlock(newblocks...,origin=origin(sb),rotation=rotation(sb))
end

function hatch(sb::SuperBlock, dhatch::Unitful.Length, bottomdir::Number, diroffset::Number)
    SuperBlock(hatch.(blocks(sb);dhatch,bottomdir,diroffset)...,
               origin=origin(sb),
               rotation=rotation(sb))
end

"""
```julia
hatch(superblock; dhatch [,bottomdir, diroffset])
```
Hatch all `Block`s contained in a `SuperBlock` with uniform hatching distance `dhatch`.
Any `Slice`s in the block which are already hatched will not be modified. `dhatch` is
the uniform slicing distance, the optional keyword arguments `bottomdir` and `diroffset`
control the hatching direction of the bottommost slice and the direction offset between slices.
"""
function hatch(sb::SuperBlock; dhatch, bottomdir=0, diroffset=pi/2)
    hatch(sb,dhatch,bottomdir,diroffset)
end
#plot recipes for blocks and superblocks
@recipe function plotblock(b::Block)
    legend --> false
    aspect_ratio --> 1
    for z in keys(b.slices)
        for rawslice in b.slices[z]
            #do coordinate transformation in xy here
            slice = translate(rotate(rawslice,rotation(b)),origin(b)[1:2])
            #add on our z coordinate, add coordinate transformation
            transz = z + ustrip(u"µm",origin(b)[3])
            if slice isa Slice
                for contour in children(slice)
                    for edge in children(contour)
                        points = boundpoints(edge)
                        pointsmat = vcat(permutedims.(points)...)
                        @series begin
                            seriestype --> :path3d
                            (pointsmat[:,1],pointsmat[:,2],repeat([transz],size(pointsmat)[1]))
                        end
                    end
                end
            else
                @assert slice isa HatchedSlice
                pointsmat = vcat(map(slice.hatchlines) do hl
                                     vcat(permutedims.(hl)...,[missing missing])
                                 end...)
                @series begin
                    seriestype --> :path3d
                    (pointsmat[:,1],pointsmat[:,2],repeat([transz],size(pointsmat)[1]))
                end
            end
        end
    end
    return nothing
end

@recipe function plotsuperblock(sb::SuperBlock)
    legend --> false
    aspect_ratio --> 1
    for b in sb.blocks
        #do coordinate tranformation
        transblock = translate(rotate(b,rotation(sb)),origin(sb))
        @series begin
            transblock
        end
    end
    return nothing
end
                         
end # module Tessen
