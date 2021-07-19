@testset "BiMADS" begin
    function test(x)
        f1(x) = (x[1] + 1) .^ 2 +(x[2] - 1) .^2 -10
        f2(x) = (x[1] + 12) .^ 2 +(x[2] - 3) .^2 + 20.
        return f1,f2
    end



    @testset "B_points" begin
        cost=[1.,2.]
        coordinates=[3.,4.]
        points=DS.B_points(cost,coordinates)
        @test points.cost == cost
        @test points.x_now==coordinates
        @test points.weight==0
    end

    @testset "BiMADS_status" begin
        s=DS.BiMADS_status()
        @test s.iteration==0
        @test s.func_evaluation==0
        @test s.total_time==0.0
        @test s.pareto_coverage==0.0
        @test s.opt_status==DS.Unoptimized
        @test s.opt_string=="Unoptimized"
        @test typeof(s.start_time)==Float64
    end

    @testset "Pareto_Evaluation" begin

    end

    @testset "Initialization" begin
        p = DSProblem(2; objective = test);
        @test p_dim(p)==2
        # if successfully splitting BOP into two SOPs
        p1, p2=DS.p_split(p)
        x = zeros(Float64, p.N)
        p.objective(x)
        f1, f2 = p.objective(x)
        @test p1.objective==f1
        @test p2.objective==f2

        @testset "Optimization for p1 and P2" begin
            SetIterationLimit(p1,1)
            @test_throws ErrorException DS.get_Bpoints(p1,p2,1)
            SetFunctionEvaluationLimit(p2,1)
            @test_throws ErrorException DS.get_Bpoints(p2,p1,2)
            #reset p1 p2 for other tests
            p1, p2=DS.p_split(DSProblem(2; objective = test))
            @test length(DS.get_Bpoints(p1,p2,1))>0
            @test length(DS.get_Bpoints(p2,p1,2))>0
        end

        # @test_throws ErrorException DS.get_Bpoints(p1,p2,1)
    end

end
