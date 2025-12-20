using Test
using AstroStates
using StaticArrays

@testset "CartesianState Constructors" begin
    
    # Test data
    x, y, z = 6778.0, 100.0, 200.0
    vx, vy, vz = 0.5, 7.5, 0.1
    posvel_vec = [x, y, z, vx, vy, vz]
    pos_vec = [x, y, z]
    vel_vec = [vx, vy, vz]
    
    @testset "Construction from 6-element Vector" begin
        cart = CartesianState(posvel_vec)
        @test cart isa CartesianState{Float64}
        @test cart.position[1] == x
        @test cart.position[2] == y
        @test cart.position[3] == z
        @test cart.velocity[1] == vx
        @test cart.velocity[2] == vy
        @test cart.velocity[3] == vz
    end
    
    @testset "Construction from SVector{6}" begin
        posvel_svec = SVector{6}(x, y, z, vx, vy, vz)
        cart = CartesianState(posvel_svec)
        @test cart isa CartesianState{Float64}
        @test cart.position == SVector(x, y, z)
        @test cart.velocity == SVector(vx, vy, vz)
    end
    
    @testset "Construction from two SVector{3}" begin
        pos_svec = SVector{3}(x, y, z)
        vel_svec = SVector{3}(vx, vy, vz)
        cart = CartesianState(pos_svec, vel_svec)
        @test cart isa CartesianState{Float64}
        @test cart.position === pos_svec
        @test cart.velocity === vel_svec
    end
    
    @testset "Construction from two AbstractVector" begin
        cart = CartesianState(pos_vec, vel_vec)
        @test cart isa CartesianState{Float64}
        @test cart.position == SVector(x, y, z)
        @test cart.velocity == SVector(vx, vy, vz)
    end
    
    @testset "Field access" begin
        cart = CartesianState(posvel_vec)
        
        # Direct field access
        @test cart.position isa SVector{3,Float64}
        @test cart.velocity isa SVector{3,Float64}
        @test length(cart.position) == 3
        @test length(cart.velocity) == 3
        
        # Individual elements
        @test cart.position[1] == x
        @test cart.position[2] == y
        @test cart.position[3] == z
        @test cart.velocity[1] == vx
        @test cart.velocity[2] == vy
        @test cart.velocity[3] == vz
    end
    
    @testset "Backward compatibility - posvel property" begin
        cart = CartesianState(posvel_vec)
        
        # posvel property should return Vector{Float64}
        pv = cart.posvel
        @test pv isa Vector{Float64}
        @test length(pv) == 6
        @test pv[1] == x
        @test pv[2] == y
        @test pv[3] == z
        @test pv[4] == vx
        @test pv[5] == vy
        @test pv[6] == vz
        @test pv == posvel_vec
    end
    
    @testset "to_vector function" begin
        cart = CartesianState(posvel_vec)
        vec = to_vector(cart)
        
        @test vec isa Vector{Float64}
        @test length(vec) == 6
        @test vec == posvel_vec
    end
    
    @testset "Type preservation" begin
        # Float32
        posvel_f32 = Float32[x, y, z, vx, vy, vz]
        cart_f32 = CartesianState(posvel_f32)
        @test cart_f32 isa CartesianState{Float32}
        @test cart_f32.position isa SVector{3,Float32}
        @test cart_f32.velocity isa SVector{3,Float32}
        
        # BigFloat
        posvel_big = BigFloat[x, y, z, vx, vy, vz]
        cart_big = CartesianState(posvel_big)
        @test cart_big isa CartesianState{BigFloat}
        @test cart_big.position isa SVector{3,BigFloat}
        @test cart_big.velocity isa SVector{3,BigFloat}
    end
    
    @testset "Mutability" begin
        cart = CartesianState(posvel_vec)
        
        # Should be mutable
        new_pos = SVector{3}(7000.0, 0.0, 0.0)
        new_vel = SVector{3}(0.0, 8.0, 0.0)
        cart.position = new_pos
        cart.velocity = new_vel
        
        @test cart.position == new_pos
        @test cart.velocity == new_vel
    end
    
    @testset "Constructor assertions" begin
        # Wrong length for 6-element constructor
        @test_throws AssertionError CartesianState([1.0, 2.0, 3.0])
        @test_throws AssertionError CartesianState([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        
        # Wrong length for separate vector constructors
        @test_throws AssertionError CartesianState([1.0, 2.0], [4.0, 5.0, 6.0])
        @test_throws AssertionError CartesianState([1.0, 2.0, 3.0], [4.0, 5.0])
    end
    
    @testset "Consistency across constructors" begin
        # All constructors should produce equivalent states
        cart1 = CartesianState(posvel_vec)
        cart2 = CartesianState(SVector{6}(posvel_vec...))
        cart3 = CartesianState(SVector{3}(pos_vec...), SVector{3}(vel_vec...))
        cart4 = CartesianState(pos_vec, vel_vec)
        
        @test cart1.position == cart2.position == cart3.position == cart4.position
        @test cart1.velocity == cart2.velocity == cart3.velocity == cart4.velocity
        @test cart1.posvel == cart2.posvel == cart3.posvel == cart4.posvel
        @test to_vector(cart1) == to_vector(cart2) == to_vector(cart3) == to_vector(cart4)
    end
end

nothing
