using ConvexPolygons
using Meshes
using Test

@testset "ConvexPolygons.jl" begin

    poly1 = ConvexPolygon([0.0, 1, 1, 0, 0], [0.0, 0, 1, 1, 0])
    poly2 = ConvexPolygon(reverse([0.0, 1, 1, 0, 0]), reverse([0.0, 0, 1, 1, 0]))


    @test nvertices(poly1) == 4
    @test nvertices(poly2) == 4

    @test area(poly1) ≈ 1.0
    @test area(poly2) ≈ 1.0

    @test normal(poly1) ≈ 1.0
    @test normal(poly2) ≈ -1.0

    @test centroid(poly1) ≈ Point(0.5, 0.5)
    @test centroid(poly2) ≈ Point(0.5, 0.5)

    @test Point(0.5, 0.5) ∈ poly1
    @test Point(0.5, 0.5) ∈ poly2

    @test Point(1.5, 1.5) ∉ poly1
    @test Point(1.5, 1.5) ∉ poly2

    @test Point(eps(), eps()) ∈ poly1
    @test Point(eps(), eps()) ∈ poly2

    @test Point(-eps(), eps()) ∉ poly1
    @test Point(-eps(), eps()) ∉ poly2

    @test Point(eps(), -eps()) ∉ poly1
    @test Point(eps(), -eps()) ∉ poly2

    @test Point(eps(), 0.5) ∈ poly1
    @test Point(eps(), 0.5) ∈ poly2

    
    @test vertices(poly1)==Point[(0.0,0.0),(1.0,0.0), (1.0,1.0), (0.0,1.0)]

    poly3 = ConvexPolygon([0.0, 1, 0, 0], [0.0, 0, 1, 0])
    
    @test nvertices(poly3) == 3

    @test area(poly3) ≈ 0.5
    @test normal(poly3) ≈ 0.5

    @test centroid(poly3) ≈ Point2(1/3, 1/3)
    @test centroid(poly3) ∈ poly3

    @test Point2(0.5+eps(), 0.5+eps()) ∉ poly3

    poly4 = ConvexPolygon([0.0, 1, 1, 0, 0], [0.0, 0, 1, 1, 0], [0.0, 0, 0, 0, 0])

    @test area(poly4) ≈ 1.0
    @test normal(poly4) ≈ Vec(0.0, 0.0, 1.0)
    
    poly5 = ConvexPolygon(reverse([0.0, 1, 1, 0,0]), reverse([0.0, 0, 1, 1,0]),
                          reverse([0.0, 0, 0, 0,0]))
    
    @test area(poly5) ≈ 1.0
    @test normal(poly5) ≈ Vec(0.0, 0.0, -1.0)
    
    
    @test pnpoly(Point(0.5, 0.5, 0.0), poly4)
    @test pnpoly(Point(0.5, 0.5, 0.0), poly5)


    @test Point(0.5, 0.5, 0.0) ∈ poly4
    @test Point(0.5, -0.5, 0.0) ∉ poly4
    @test Point(0.5, -eps(), 0.0) ∉ poly4
    @test Point(eps(), eps(), 0.0) ∈ poly4
    
    @test Point(0.5, 0.5, 1e-5) ∉ poly4
    @test Point(0.5, 0.5, sqrt(eps())*5) ∉ poly4
    @test Point(0.5, 0.5, sqrt(eps())/5) ∈ poly4

    poly6 = ConvexPolygon([0.0, 1, 1, 0,0], [0.0, 0, 1, 1,0],
                          [0.0, 0, 1, 1,0])
    A = area(poly6)
    @test A ≈ sqrt(2)
    n = normal(poly6)
    @test n ./ A ≈ Vec(0.0, -sqrt(2)/2, sqrt(2)/2)
                 
    
    @test Point(0.5, 0.5, 0.5) ∈ poly6
    @test Point(-eps(), 0.5, 0.5) ∉ poly6
    @test Point(0.5, -eps(), -eps()) ∉ poly6
    @test Point(eps(), eps(), eps()) ∈ poly6
    
    @test Point(0.5, 0.5, 0.5+1e-5) ∉ poly6
    @test Point(0.5, 0.5, 0.5 + sqrt(eps())*5) ∉ poly6
    @test Point(0.5, 0.5, 0.5 + sqrt(eps())/5) ∈ poly6
    
    pt = Point{3,Float64}[(0.0,0.0, 0.0), (1.0,0.0,0.0), (0.0,1.0, 1.0)]
    tri = Triangle(pt)
    
    @test pnpoly(Point(0.2, 0.2, 0.2), tri)
    @test !pnpoly(Point(0.2, 0.2, 0.21), tri)
    @test !pnpoly(Point(-eps(), -eps(), -eps()), tri)
    @test pnpoly(Point(0.2, 0.2, 0.2+0.2*√eps()), tri)
    @test !pnpoly(Point(0.2, 0.2, 0.2+5*√eps()), tri)
    

    
end
