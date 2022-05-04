module ConvexPolygons

using Meshes
#import Meshes: Polygon, Chain, Point, Point2, Point3, normal, centroid, area
#import Meshes: coordinates, vertices, nvertices, measure
#import Meshes: isclosed

export ConvexPolygon, area, centroid, normal, coordinates, convex2poly, vertices
export pnpoly, nvertices, measure, issimple, hasholes, chains


"""
`ConvexPolygon(pts)`

`ConvexPolygon(x, y)`

`ConvexPolygon(x, y, z)`


Creates a convex polygon from the vertices. It does not check anything! 
It just stores the vertices as a contour.
"""
struct ConvexPolygon{Dim,T,C<:Chain{Dim,T}} <: Polygon{Dim,T}
    contour::C
    function ConvexPolygon(contour::Chain{Dim,T}) where {Dim,T}
        @assert isclosed(contour) "invalid outer chain"
        new{Dim,T,Chain{Dim,T}}(contour)
    end
end

    

ConvexPolygon(x::AbstractVector{T}, y::AbstractVector{T}) where {T} =
    ConvexPolygon([Point{2,T}(xx, yy) for (xx,yy) in zip(x,y)])

ConvexPolygon(x::AbstractVector{T},
              y::AbstractVector{T}, z::AbstractVector{T}) where {T} =
    ConvexPolygon([Point{3,T}(xx, yy,zz) for (xx,yy,zz) in zip(x,y,z)])

ConvexPolygon(pts::AbstractVector{P}) where {P<:Point} =
    ConvexPolygon(Chain(pts))

ConvexPolygon(pts::AbstractVector{TP}) where {TP<:Tuple} =
    ConvexPolygon(Chain(Point.(pts)))


"""
`vertices(p::ConvexPolygon)`

Returns the coordinates of the vertices of a [`ConvexPolygon`].
"""
Meshes.vertices(p::ConvexPolygon) = vertices(p.contour)
Meshes.nvertices(p::ConvexPolygon) = nvertices(p.contour)

Meshes.issimple(p::ConvexPolygon) = true
Meshes.hasholes(p::ConvexPolygon) = false
Meshes.chains(p::ConvexPolygon) = [p.contour]

"""
`normal(p::ConvexPolygon{2,T})`

Computes the surface area of a 2D [`ConvexPolygon`](@ref).

If the area is positive, the normal is pointed upwards and the vertices
are in counter clockwise direction.
"""
function Meshes.normal(p::ConvexPolygon{2,T}) where {T}
    nv = nvertices(p)
    verts = vertices(p)
    x₁,y₁ = coordinates(verts[end])
    x₂,y₂ = coordinates(verts[begin])
    
    A = x₁*y₂ - x₂*y₁
    
    for i in 2:nv
        x₁,y₁ = coordinates(verts[i-1])
        x₂,y₂ = coordinates(verts[i])
        A += x₁*y₂ - x₂*y₁
    end

    return A/2
end

import Meshes.area

"""
`normal(p::ConvexPolygon{2,T})`

Computes the surface area of a [`ConvexPolygon`](@ref).

For 2d polygons it returns a positive area  if the contour
is oriented in counter-clockwise direction and negative area if
clockwise.

"""
Meshes.area(p::ConvexPolygon{2,T}) where {T} = abs(normal(p))


crossprod(u, v) = (u[2]*v[3] - u[3]*v[2],
                   u[3]*v[1] - u[1]*v[3],
                   u[1]*v[2] - u[2]*v[1])
                                   
function Meshes.normal(p::ConvexPolygon{3,T}) where {T}
    v = vertices(p)
    nv = length(v)
    x,y,z = crossprod(coordinates(v[end]), coordinates(v[begin]))

    for i in 2:nv
        prd = crossprod(coordinates(v[i-1]), coordinates(v[i]))
        x += prd[1]
        y += prd[2]
        z += prd[3]
    end

    return Vec(x/2, y/2, z/2)
    
end

"""
`area(p::ConvexPolygon{2,T})`

Computes the surface area of a 2D [`ConvexPolygon`](@ref).
"""
Meshes.area(p::ConvexPolygon{3}) = hypot(normal(p)...)

"""
`centroid(p)`

Computes the centroid of a [`ConvexPolygon`](@ref).
"""                      
function Meshes.centroid(p::ConvexPolygon{2,T}) where {T}

    A = normal(p)

    v = vertices(p)
    nv = length(v)

    # xᵢyᵢ₊₁ - xᵢ₊₁yᵢ  First index = last index
    x₁, y₁ = coordinates(v[end])
    x₂, y₂ = coordinates(v[begin])
    tmp = x₁*y₂ - x₂*y₁
    
    Cx = (x₁ + x₂) * tmp
    Cy = (y₁ + y₂) * tmp

    for i in firstindex(v):lastindex(v)-1 
        x₁, y₁ = coordinates(v[i])
        x₂, y₂ = coordinates(v[i+1])
        tmp = x₁*y₂ - x₂*y₁
        Cx += (x₁ + x₂) * tmp
        Cy += (y₁ + y₂) * tmp
    end

    return Point{2,T}(Cx/6A, Cy/6A)
end

function Meshes.centroid(p::ConvexPolygon{3,T}) where {T}
    v = vertices(p)
    nv = length(v)
    nrm = normal(p)
    A = hypot(nrm...)
    n⃗ = nrm ./ A
    
    C = Vec{3,T}(0,0,0)
    v₀ = coordinates(v[begin])
    for i = (firstindex(v)+1):(lastindex(v)-1)
        v₁ = coordinates(v[i])
        v₂ = coordinates(v[i+1])
                          
        centr = (v₀ + v₁ + v₂) ./ 3
        Aᵢ⃗ = crossprod(v₁-v₀, v₂-v₀) ./ 2
        s = sign(sum(Aᵢ⃗ .* n⃗))
        C += s * centr .* hypot(Aᵢ⃗...)
    end
    return C
    
end

convex2poly(p::ConvexPolygon) = Polygon(coordinates(p))


"""
`pnpoly(p, poly)`

Checks whether a point is inside a polygon.
"""
function pnpoly(p::Point{2,T}, poly::ConvexPolygon{2,T}) where {T}
    verts = vertices(poly.contour)
    j = lastindex(verts)
    test = false
    x, y = coordinates(p)
    for i in eachindex(verts)
        vi = coordinates(verts[i])
        vj = coordinates(verts[j])
        if ((vi[2] > y) != (vj[2] > y)) &&
            ( x ≤ (vj[1]-vi[1]) * (y - vi[2]) / (vj[2] - vi[2]) + vi[1])
            test = !test
        end
        j = i
    end
    return test
end

Base.in(p::Point, poly::ConvexPolygon) where {T} = pnpoly(p, poly)

function project_point(p::Point{3}, p₀::Point{3}, P::Plane)
    w = p-p₀
    return Point(w⋅P.u, w⋅P.v)
end

import CircularArrays: CircularVector
function pnpoly(p::Point{3,T}, poly::ConvexPolygon{3,T}) where {T}
    # 3d case
    # First we needto determine if the point lies in the same plane
    # as the polygon

    n = normal(poly)  # Normal To the polygon
    A = hypot(n...)  # Area
    L = sqrt(A)
    verts = vertices(poly)
    v₀ = verts[1]
    δ = abs(n⋅(p-v₀)) / A
    e = sqrt(eps(L))
    if δ > e
        return false
    end

    # The point is in the same plane as the polygon
    # Project the polygon on the plane
    # Since the polygon is complex, we can check each triangle
    nv = nvertices(poly)
    P = Plane(p, n)

    foundtri = false
    u₀ = Point(0.0,0.0) # reference point, first vertex in 2d
    p₂ = project_point(p, v₀, P)
    for i in 3:nv
        u₁ = project_point(verts[i-1], v₀, P)
        u₂ = project_point(verts[i],   v₀, P)
        tri = Triangle(u₀, u₁, u₂)

        if p₂ ∈ tri
            return true
        end
    end
    
    return false
    
end


"""
`pnpoly(p, tri)`

Checks whether a point is inside a Triangle.
"""
pnpoly(p::Point{2,T}, tri::Triangle{2,T}) where {T} = p ∈ tri

function pnpoly(p::Point{3,T}, tri::Triangle{3,T}) where {T}



    n = normal(tri)  # Normal To the triangle
    A = area(tri)  # Area
    L = sqrt(A)
    v₀ = tri.vertices[1]
    δ = abs(n⋅(p-v₀)) / A
    e = sqrt(eps(L))
    if δ > e
        return false
    end

    # The point is in the same plane as the polygon
    # Project the polygon on the plane
    # Since the polygon is complex, we can check each triangle
    P = Plane(p, n)
    u₀ = Point{2,T}(0,0)
    u₁ = project_point(tri.vertices[2], v₀, P)
    u₂ = project_point(tri.vertices[3], v₀, P)

    p₂ = project_point(p, v₀, P)

    tri₂ = Triangle(u₀, u₁, u₂)

    return p₂ ∈ tri₂
    

end




end
