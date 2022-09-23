using FlowCytometry
using Test
using DataFrames
using JLD

# @testset "Gating" begin

#     # @test_nowarn FlowCytometryGate(("A","B"),[(0,0),(4,0),(4,4),(2,2),(0,4)])
#     gate = FlowCytometryGate(("A","B"),[(0,0),(4,0),(4,4),(2,2),(0,4)])

#     @test FlowCytometry.isInsideGate((-1,-1),gate) == false
#     @test FlowCytometry.isInsideGate((2,-1),gate) == false
#     @test FlowCytometry.isInsideGate((5,-1),gate) == false

#     @test FlowCytometry.isInsideGate((-1,0),gate) == false
#     @test FlowCytometry.isInsideGate((0,0),gate) == true
#     @test FlowCytometry.isInsideGate((2,0),gate) == true
#     @test FlowCytometry.isInsideGate((4,0),gate) == true
#     @test FlowCytometry.isInsideGate((5,0),gate) == false

#     @test FlowCytometry.isInsideGate((-1,1),gate) == false
#     @test FlowCytometry.isInsideGate((0,1),gate) == true
#     @test FlowCytometry.isInsideGate((2,1),gate) == true
#     @test FlowCytometry.isInsideGate((4,1),gate) == true
#     @test FlowCytometry.isInsideGate((5,1),gate) == false

#     @test FlowCytometry.isInsideGate((-1,2),gate) == false
#     @test FlowCytometry.isInsideGate((0,2),gate) == true
#     @test FlowCytometry.isInsideGate((2,2),gate) == true
#     @test FlowCytometry.isInsideGate((4,2),gate) == true
#     @test FlowCytometry.isInsideGate((5,2),gate) == false

#     @test FlowCytometry.isInsideGate((-1,3),gate) == false
#     @test FlowCytometry.isInsideGate((0.1,3),gate) == true
#     @test FlowCytometry.isInsideGate((0,3),gate) == true
#     @test FlowCytometry.isInsideGate((2,3),gate) == false
#     @test FlowCytometry.isInsideGate((3.9,3),gate) == true
#     @test FlowCytometry.isInsideGate((4,3),gate) == true
#     @test FlowCytometry.isInsideGate((5,3),gate) == false

#     @test FlowCytometry.isInsideGate((-1,4),gate) == false
#     @test FlowCytometry.isInsideGate((0,4),gate) == true
#     @test FlowCytometry.isInsideGate((2,4),gate) == false
#     @test FlowCytometry.isInsideGate((4,4),gate) == true
#     @test FlowCytometry.isInsideGate((5,4),gate) == false

#     @test FlowCytometry.isInsideGate((-1,5),gate) == false
#     @test FlowCytometry.isInsideGate((0,5),gate) == false
#     @test FlowCytometry.isInsideGate((2,5),gate) == false
#     @test FlowCytometry.isInsideGate((4,5),gate) == false
#     @test FlowCytometry.isInsideGate((5,5),gate) == false

# end

# @testset "FlowCytometryExperiment" begin

#     @test_nowarn FlowCytometryExperiment(ones(10,10))
#     @test_nowarn FlowCytometryExperiment(ones(10,10),obs=DataFrame(:a=>1:10))
#     @test_nowarn FlowCytometryExperiment(ones(10,10),var=DataFrame(:a=>1:10))
#     @test_throws ErrorException FlowCytometryExperiment(ones(10,10),obs=DataFrame(:a=>1:3))
#     @test_throws ErrorException FlowCytometryExperiment(ones(10,10),var=DataFrame(:a=>1:3))

#     @test_nowarn loadFCExperiment("testdata/BD-FACS-Aria-II.fcs")

#     @test_nowarn begin m = FlowCytometryExperiment(ones(10,10)); m.obs = DataFrame(:a=>ones(10)) end
#     @test_nowarn begin m = FlowCytometryExperiment(ones(10,10)); m.var = DataFrame(:a=>ones(10)) end
#     @test_nowarn begin m = FlowCytometryExperiment(ones(10,10)); m.X = zeros(10,10) end
#     @test_throws ErrorException begin m = FlowCytometryExperiment(ones(10,10)); m.obs = DataFrame(:a=>ones(9)) end
#     @test_throws ErrorException begin m = FlowCytometryExperiment(ones(10,10)); m.var = DataFrame(:a=>ones(9)) end
#     @test_throws ErrorException begin m = FlowCytometryExperiment(ones(10,10)); m.X = zeros(9,10) end
#     @test_throws ErrorException begin m = FlowCytometryExperiment(ones(10,10)); m.X = zeros(10,9) end

#     @test_nowarn begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; removeCells(m,keep) end
#     @test_nowarn begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; removeCells!(m,keep) end
#     @test_nowarn begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; removeChannels(m,keep) end
#     @test_nowarn begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; removeChannels!(m,keep) end

#     @test begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; m2 = removeCells(m,keep); size(m2.X) == (1,10) end
#     @test begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; m2 = removeCells(m,keep); size(m2.obsm["A"]) == (1,10) end
#     @test begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; removeCells!(m,keep); size(m.X) == (1,10) end
#     @test begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; removeCells!(m,keep); size(m.obsm["A"]) == (1,10) end
#     @test begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; m2 = removeChannels(m,keep); size(m2.X) == (10,1) end
#     @test begin m = FlowCytometryExperiment(ones(10,10)); m.obsm["A"] = zeros(10,10); keep = fill(false, 10); keep[1] = true; removeChannels!(m,keep); size(m.X) == (10,1) end

#     @test begin m = FlowCytometryExperiment([1 2 3].*ones(10,3),var=DataFrame(:channel=>["A","B","C"])); m["B"] == 2 .*ones(10) end
#     @test begin m = FlowCytometryExperiment([1 2 3].*ones(10,3),var=DataFrame(:channel=>[:A,:B,:C])); m[:B] == 2 .*ones(10) end
#     @test begin m = FlowCytometryExperiment([1 2 3].*ones(10,3),var=DataFrame(:channel=>[1,2,3])); m[2] == 2 .*ones(10) end

# end

@testset "DimensionalityReduction" begin

    fcs = FlowCytometryExperiment([rand(10,15);rand(10,15).+[10 0 0 0 0 0 0 0 0 0 0 0 0 0 0]])
    fcs.var[:,"useChannels"] = [true,true,true,true,true,false,false,false,false,false,false,false,false,false,false]

    @test_nowarn DimensionalityReduction.pca!(fcs)
    @test begin DimensionalityReduction.pca!(fcs); println(size(fcs.obsm["pca"])); size(fcs.obsm["pca"])==(20,14) end
    @test_nowarn DimensionalityReduction.pca!(fcs,key_used_channels="useChannels")
    @test begin DimensionalityReduction.pca!(fcs); size(fcs.obsm["pca"])==(20,4) end

end