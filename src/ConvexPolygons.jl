module ConvexPolygons

import GeometryBasics: AbstractPoint, Point, Point2, Point3, AbstractPolygon, area



export ConvexPolygon, area, centroid, normal, coordinates, convex2poly

"""
`ConvexPolygon(pts)`

`ConvexPolygon(x, y)`

`ConvexPolygon(x, y, z)`


Creates a convex polygon from the vertices. It does not check anything! 
It just stores the vertices as a contour.
"""
struct ConvexPolygon{Dim,T,P<:AbstractPoint{Dim,T},L<:AbstractVector{P}} <: AbstractPolygon{Dim,T}
    contour::L
end

ConvexPolygon(x::AbstractVector, y::AbstractVector) =
    ConvexPolygon([Point2(xx, yy) for (xx,yy) in zip(x,y)])

ConvexPolygon(x::AbstractVector, y::AbstractVector, z::AbstractVector) =
    ConvexPolygon([Point3(xx, yy,zz) for (xx,yy,zz) in zip(x,y,z)])

#ConvexPolygon(pts::E) where {E<:AbstractVector{P}} where {P<:AbstractPoint{Dim,T}} where {Dim,T} = ConvexPolygon{Dim,T,P,E}(pts)

import GeometryBasics.coordinates
"""
`coordinates(p::ConvexPolygon)`

Returns the coordinates of the vertices of a [`ConvexPolygon`].
"""
coordinates(p::ConvexPolygon) = p.contour

"""
`area(p::ConvexPolygon{2,T})`

Computes the surface area of a 2D [`ConvexPolygon`](@ref).
"""
area(p::ConvexPolygon{2,T}) where {T} = area(coordinates(p))

"""
`normal(p::ConvexPolygon{2,T})`

Computes the surface area of a [`ConvexPolygon`](@ref).
"""
normal(p::ConvexPolygon{2,T}) where {T} = area(coordinages(p))

crossprod(u::Point3, v::Point3) = (u[2]*v[3] - u[3]*v[2],
                                   u[3]*v[1] - u[1]*v[3],
                                   u[1]*v[2] - u[2]*v[1])
                                   
function normal(p::ConvexPolygon{3,T}) where {T}
    v = coordinates(p)
    nv = length(v)
    x,y,z = crossprod(v[end], v[begin])

    for i in 2:nv
        prd = crossprod(v[i-1], v[i])
        x += prd[1]
        y += prd[2]
        z += prd[3]
    end

    return Point3(x/2, y/2, z/2)
    
end

"""
`area(p::ConvexPolygon{2,T})`

Computes the surface area of a 2D [`ConvexPolygon`](@ref).
"""
area(p::ConvexPolygon{3}) = hypot(normal(p)...)

"""
`centroid(p)`

Computes the centroid of a [`ConvexPolygon`](@ref).
"""                      
function centroid(p::ConvexPolygon{2,T}) where {T}

    A = area(p)

    v = coordinates(p)
    nv = length(v)

    # xᵢyᵢ₊₁ - xᵢ₊₁yᵢ  First index = last index
    tmp = v[end][1]*v[begin][2] - v[begin][1]*v[end][2]
    
    Cx = (v[end][1] + v[begin][1]) * tmp
    Cy = (v[end][2] + v[begin][2]) * tmp

    for i in firstindex(v):lastindex(v)-1
        tmp = v[i][1]*v[i+1][2] - v[i+1][1]*v[i][2]
        Cx += (v[i][1] + v[i+1][1]) * tmp
        Cy += (v[i][2] + v[i+1][2]) * tmp
    end

    return Point2{T}(Cx/6A, Cy/6A)
end

function centroid(p::ConvexPolygon{3,T}) where {T}
    v = coordinates(p)
    nv = length(v)
    nrm = normal(p)
    A = hypot(nrm...)
    n⃗ = nrm ./ A
    
    C = Point3{T}(0,0,0)
    v₀ = v[begin]
    for i = (firstindex(v)+1):(lastindex(v)-1)
        centr = (v₀ + v[i] + v[i+1]) ./ 3
        Aᵢ⃗ = crossprod(v[i]-v₀, v[i+1]-v₀) ./ 2
        s = sign(sum(Aᵢ⃗ .* n⃗))
        C += s * centr .* hypot(Aᵢ⃗...)
    end
    return C
    
end

convex2poly(p::ConvexPolygon) = Polygon(coordinates(p))

end
