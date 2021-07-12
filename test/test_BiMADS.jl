@testset "BiMADS" begin

    @testset "B_points" begin
        cost=[1.,2.]
        coordinates=[3.,4.]
        points=B_points(cost,coordinates)
        @test points.cost == cost
        @test points.x_now==coordinates
        @test points.weight==0
    end

    @testset "BiMADS_status" begin
        # status=BiMADS_status()
        # @test s.iteration==0
        # @test s.func_evaluation==0
        # @test s.total_time==0.0
        # @test s.pareto_coverage==0.0
        # @test s.opt_status==Unoptimized
        # @test s.opt_string=="Unoptimized"
        # @test typeof(s.start_time) isa Float64
    end



end
