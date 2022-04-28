using ConvexPolygons
using GeometryBasics
using Test

@testset "ConvexPolygons.jl" begin

    poly1 = ConvexPolygon([0.0, 1, 1, 0], [0.0, 0, 1, 1])
    poly2 = ConvexPolygon(reverse([0.0, 1, 1, 0]), reverse([0.0, 0, 1, 1]))


    @test length(poly1) == 4
    @test length(poly2) == 4
    @test numverts(poly1) == 4
    @test numverts(poly2) == 4

    @test area(poly1) ≈ 1.0
    @test area(poly2) ≈ -1.0

    @test centroid(poly1) ≈ Point2(0.5, 0.5)
    @test centroid(poly2) ≈ Point2(0.5, 0.5)

    @test Point2(0.5, 0.5) ∈ poly1
    @test Point2(0.5, 0.5) ∈ poly2

    @test Point2(1.5, 1.5) ∉ poly1
    @test Point2(1.5, 1.5) ∉ poly2

    @test Point2(eps(), eps()) ∈ poly1
    @test Point2(eps(), eps()) ∈ poly2

    @test Point2(-eps(), eps()) ∉ poly1
    @test Point2(-eps(), eps()) ∉ poly2

    @test Point2(eps(), -eps()) ∉ poly1
    @test Point2(eps(), -eps()) ∉ poly2

    @test Point2(eps(), 0.5) ∈ poly1
    @test Point2(eps(), 0.5) ∈ poly2

    
    @test coordinates(poly1)==Point2{Float64}[(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)]

    poly3 = ConvexPolygon([0.0, 1, 0], [0.0, 0, 1])
    
    @test length(poly3) == 3
    @test numverts(poly3) == 3

    @test area(poly3) ≈ 0.5
    @test normal(poly3) ≈ 0.5

    @test centroid(poly3) ≈ Point2(1/3, 1/3)
    @test centroid(poly3) ∈ poly3

    @test Point2(0.5+eps(), 0.5+eps()) ∉ poly3

    poly4 = ConvexPolygon([0.0, 1, 1, 0], [0.0, 0, 1, 1], [0.0, 0, 0, 0])

    @test area(poly4) ≈ 1.0
    @test normal(poly4) ≈ Point3(0.0, 0.0, 1.0)
    
    poly5 = ConvexPolygon(reverse([0.0, 1, 1, 0]), reverse([0.0, 0, 1, 1]),
                          reverse([0.0, 0, 0, 0]))
    
    @test area(poly5) ≈ 1.0
    @test normal(poly5) ≈ Point3(0.0, 0.0, -1.0)
    
    
    
    
    
    

    
end
