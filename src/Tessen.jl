module Tessen

using LinearAlgebra, Unitful, RecipesBase

export LineEdge, ArcEdge, Contour, Slice, translate, rotate
export Block, SuperBlock, blocks, hatch

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
    if isapprox((v1./v2)...;atol)
        #the lines are parallel, these are colinear if either p1==p2
        samepoints = all(isapprox.(p1,p2;atol))
        #or if (p1-p2) is parallel to v2
        parallel = isapprox(((p1-p2)./v2)...;atol)
        if samepoints || parallel
            return [le.p1,le.p2]
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
    ArcEdge(newc,ae.r,ae.startangle,ae.stopangle)
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

#make an abstract supertype for objects which are just lists of child objects
abstract type CompoundObject end

"""
```julia
children(obj)
```
Get the children of a `CompoundObject`
"""
function children end

#rotating and translating compound objects is just rotating/translating all the children
function rotate(c::CO,args...) where {CO <: CompoundObject}
    CO([rotate(ch,args...) for ch in children(c)])
end

function translate(c::CO,args...) where {CO <: CompoundObject}
    CO([translate(ch,args...) for ch in children(c)])
end

#Contours are a list of edges
struct Contour <: CompoundObject
    edges :: Vector{<:Edge}
end

children(c::Contour) = c.edges

#slices consist of a list of contours
"""
```julia
Slice(contours)
```
A `Slice` represents a set of coplanar contours
"""
struct Slice <: CompoundObject
    contours :: Vector{Contour}
end

children(s::Slice) = s.contours
function intersections(x::Union{Contour,Slice},hl::HatchLine)
    childintersections = [intersections(cx,hl) for cx in children(x)]
    #get rid of nothing entries and vcat into a flat Vector
    inters = vcat(filter(childintersections) do ci
                      !isnothing(ci)
                  end...)
    isempty(inters) ? nothing : inters
end

function boundbox(co::CompoundObject)
    childboxes = boundbox.(children(co))
    map(zip([minimum,maximum],[1,2])) do (m,corner)
        map([1,2]) do dim
            m(b[corner][dim] for b in childboxes)
        end
    end
end

@recipe function plotco(co::CompoundObject)
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
HatchedSlice(slice,dhatch,hatchdir)
```
Create a `HatchedSlice` object representing `slice` hatched with a hatching distance `dhatch`.
`hatchdir` is the direction in which hatch lines will be laid down, the lines themselves will
be perpendicular to this direction (which should be provided in radians). The `hatchlines`
field of this struct is a vector of point pairs describing line segments presented in order to
neatly hatch the provided `slice` in a zigzag pattern.
"""
struct HatchedSlice
    hatchlines :: Vector{Vector{Vector{<:Number}}}
    function HatchedSlice(s::Slice,dhatchunits::Unitful.Length,hatchdir::Number)
        #=============================Hatching strategy=========================
        1) get the bounding box of our slice
        2) trace this box with a contour
        3) place the first hatchline on the corner of the bounding box where intersections between
           this line and the bounding box are closest together (if these intersections are far
           apart we're missing part of the shape
        4) continue placing hatchlines until we place a hatchline that doesn't intersect our box
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
        #try placing initial hatch lines at every corner and pick the corner where the 2
        #intersections with our bounding box are closest together
        (_,startpointindex) = findmin(bboxcorners) do bc
            hl = HatchLine(bc,vhatch)
            inters = [intersections(be,hl) for be in bboxedges]
            #remove nothing entries corresponding to edges with no intersection
            filter!(inters) do i
                !isnothing(i)
            end
            inters=vcat(inters...)
            #we can end up with more than two intersections since we are directly on a corner
            #but this will never be a good starting point
            if length(inters) > 2
                return Inf
            end
            #could also freak out if we're hatching parallel to an edge, but I think we can
            #just commit to not doing that.
            #intersections now containts two parametric coordinates for the two intersections
            #along hl, we want to compare the distances
            @assert length(inters) == 2
            abs(-(inters...))
        end
        startpoint = bboxcorners[startpointindex]
        #get the distance between startpoints for neighboring hatchlines. these points will be
        #dhatch apart along hatchdir
        poffset = dhatchunits*[cos(hatchdir),sin(hatchdir)]
        #it is possible this poffset is in the wrong direction (going away from the bounding box)
        #check that there are intersections for the second hatchline, if there are none, reverse
        #direction
        secondpoint = startpoint + poffset
        secondpointints = [intersections(be,HatchLine(secondpoint,vhatch)) for be in bboxedges]
        correctdir = any(secondpointints) do spi
            !isnothing(spi)
        end
        if !correctdir
            poffset *= -1
        end
        #lay down our hatchlines
        curpoint = startpoint
        hlines = HatchLine[]
        while true
            hl = HatchLine(curpoint,vhatch)
            #if this doesn't intersect our bounding box, we're done
            inters = [intersections(be,hl) for be in bboxedges]
            if all(isnothing,inters)
                break
            end
            #otherwise pop it into hlines and move curpoint
            push!(hlines,hl)
            curpoint += poffset
        end
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
            #we're having an issue when we're tangent to arcs, I think we can just remove
            #tangent points from intersection(arcedge,hl)
            @assert iseven(length(inters))
            #we will sort these into an order based on `rev`
            sort!(inters;rev)
            #now need to turn this into points rather than parametric coords
            [pointalong(hl,i) for i in inters]
        end
        #remove nothing entries
        filter!(intervec) do iv
            !isnothing(iv)
        end
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
        new(hatchlines)
    end
end

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

function translate(lf::LF,displacement::Vector{<:Unitful.Length}) where {LF <: LocalFrame}
    LF(lf,displacement,0)
end

function rotate(lf::LF,amount,point::Vector{<:Unitful.Length}) where {LF <: LocalFrame}
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
Assemble a series of `Slices` into a `Block`. `Slice` objects should be
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
        #strip units, we will represent everything internally in microns
        rawslice = [ustrip(u"µm",sp.first) => sp.second for sp in slicepairs]
        #change our slice pairs into a dict
        allz = [rs.first for rs in rawslice]
        slicetype = eltype([rs.second for rs in rawslice])
        @assert slicetype <: Union{Slice,HatchedSlice} "Blocks are built from slices"
        #initialize all required slice vectors in the dict
        slices = Dict(z => Vector{slicetype}() for z in allz)
        #now push all the slices onto the correct vector
        for (z,s) in rawslice
            push!(slices[z],s)
        end
        new{slicetype}(ustrip.(u"µm",origin),rotation,slices)
    end
    #LocalFrame innner constructor to make translation/rotation easier
    function Block{T}(b::Block{T},translation,rotation) where {T}
        #if translation isn't in 3D, go ahead and assume we mean movement in xy
        if length(translation) == 2
            push!(translation,0u"µm")
        end
        new{T}(ustrip.(u"µm",(b.origin * 1u"µm") + translation),
               b.rotation + rotation, b.slices)
    end
end

#make kwargs optional
Block(slices...;origin=[0u"µm",0u"µm",0u"µm"],rotation=0) = Block(origin,rotation,slices...)
#required localframe methods
origin(b::Block) = b.origin * 1u"µm"
rotation(b::Block) = b.rotation

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
    blocks::Vector{LocalFrame{T}}

    function SuperBlock(origin::Vector{<:Unitful.Length},rotation::Number,
                        blocks::LocalFrame{T}...) where {T}
        new{T}(ustrip.(u"µm",origin),rotation,collect(blocks))
    end

    #localframe constructor
    function SuperBlock{T}(sb::SuperBlock{T},translation,rotation) where {T}
        #if translation isn't in 3D, go ahead and assume we mean movement in xy
        if length(translation) == 2
            push!(translation,0)
        end
        new{T}(ustrip.(u"µm",(sb.origin * 1u"µm") + translation),
               sb.rotation + rotation, sb.blocks)
    end
end
SuperBlock(blocks::LocalFrame...;origin=[0u"µm",0u"µm"],rotation=0) = SuperBlock(origin,rotation,blocks...)
origin(sb::SuperBlock) = sb.origin * 1u"µm"
rotation(sb::SuperBlock) = sb.rotation

#plot recipes for blocks and superblocks
@recipe function plotblock(b::Block)
    legend --> false
    aspect_ratio --> 1
    for z in keys(b.slices)
        for rawslice in b.slices[z]
            #do coordinate transformation in xy here
            slice = translate(rotate(rawslice,rotation(b)),origin(b)[1:2])
            for contour in children(slice)
                for edge in children(contour)
                    points = boundpoints(edge)
                    pointsmat = vcat(permutedims.(points)...)
                    @series begin
                        seriestype --> :path3d
                        #add on our z coordinate, add coordinate transformation
                        transz = z + ustrip(u"µm",origin(b)[3])
                        (pointsmat[:,1],pointsmat[:,2],repeat([transz],size(pointsmat)[1]))
                    end
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
