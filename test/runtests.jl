using FlowCytometry
using Test
using DataFrames
using JLD
using CSV

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

#     @test begin m = FlowCytometryExperiment([1 2 3].*ones(10,3),channels=["A","B","C"]); m["B"] == 2 .*ones(10) end

# end

@testset "IO" begin

    @test_nowarn loadFCExperiment("testdata/BD-FACS-Aria-II.fcs")
    @test_nowarn loadFCControls(Dict([
                                    "testdata/FlowRepository_FR-FCM-Z2SS_files/Compensation Controls_APC-R700 Stained Control_007.fcs"=>"APC-R700-A",
                                    "testdata/FlowRepository_FR-FCM-Z2SS_files/Compensation Controls_APC-H7 Stained Control_008.fcs"=>"APC-H7-A"
                                    ]
                                    )
                                )

    @test_nowarn saveH5fcs(FlowCytometryExperiment(rand(10,10)),"test")
    @test_nowarn begin 
        fcs = FlowCytometryExperiment(rand(10,10))
        fcs.obs[:,"h"] .= 1
        fcs.var[:,"h"] .= 1
        fcs.obsm["f"] = rand(3,3)
        fcs.layers["5"] = rand(2,2)
        fcs.gates["5"] = FlowCytometryGate(("1","2"),[(1,2),(2,3),(3,2)])
        fcs.uns["f"] = Dict(["d"=>[(1,2)]])
        saveH5fcs(fcs,"test") 
    end
    @test_nowarn begin 
        fcs = FlowCytometryExperiment(rand(10,10))
        fcsControl = FlowCytometryControl()
        fcsControl.channels = fcs.channels
        fcsControl.controls["a"] = fcs
        fcsControl.compensationMatrix = nothing
        fcsControl.uns["s"] = "s"
        saveH5fcs(fcsControl,"test") 
    end
    @test_nowarn begin saveH5fcs(FlowCytometryExperiment(rand(10,10)),"test"); loadH5fcs("test.h5fcs") end 
    @test_nowarn begin 
        fcs = FlowCytometryExperiment(rand(10,10))
        fcs.obs[:,"h"] .= 1
        fcs.var[:,"h"] .= 1
        fcs.obsm["f"] = rand(3,3)
        fcs.layers["5"] = rand(2,2)
        fcs.gates["5"] = FlowCytometryGate(("1","2"),[(1,2),(2,3),(3,2)])
        fcs.uns["f"] = Dict(["d"=>[(1,2)]])
        saveH5fcs(fcs,"test") 
        loadH5fcs("test.h5fcs")
    end
    @test begin fcs = FlowCytometryExperiment(rand(10,10)); 
                saveH5fcs(fcs,"test"); 
                fcsLoaded = loadH5fcs("test.h5fcs"); 
                all(fcs.channels .== fcsLoaded.channels) && 
                all(fcs.obs[:,:cell] .== fcsLoaded.obs[:,:cell]) &&
                all(fcs.var[:,:channel] .== fcsLoaded.var[:,:channel])
        end 
    @test begin 
        fcs = FlowCytometryExperiment(rand(10,10))
        fcs.obs[:,"h"] .= 1
        fcs.var[:,"h"] .= 1
        fcs.obsm["f"] = rand(3,3)
        fcs.layers["5"] = rand(2,2)
        fcs.gates["5"] = FlowCytometryGate(("1","2"),[(1,2),(2,3),(3,2)])
        fcs.uns["f"] = Dict(["d"=>[(1,2)]])
        saveH5fcs(fcs,"test") 
        fcsLoaded = loadH5fcs("test.h5fcs"); 
        all(fcs.channels .== fcsLoaded.channels) && 
        all(fcs.obs[:,:h] .== fcsLoaded.obs[:,:h]) &&
        all(fcs.var[:,:h] .== fcsLoaded.var[:,:h]) &&
        all(fcs.obsm["f"] .== fcsLoaded.obsm["f"]) &&
        all(fcs.layers["5"] .== fcsLoaded.layers["5"]) &&
        all(fcs.gates["5"].channels .== fcsLoaded.gates["5"].channels) &&
        all(fcs.gates["5"].polygon .== fcsLoaded.gates["5"].polygon) &&
        all(fcs.uns["f"]["d"] .== fcsLoaded.uns["f"]["d"])
end

end

# @testset "Compensation" begin
    
#     @test_nowarn begin fcs = loadFCControls("testdata/FlowRepository_FR-FCM-Z2SS_files")
#         channelnames = CSV.read("testdata/FlowRepository_FR-FCM-Z2SS_files/attachments/fcs_control.csv",DataFrame)
#         channelnames = Dict([i=>String(j) for (i,j) in eachrow(channelnames[:,["filename","dye"]])])
#         renameControl!(fcs,channelnames)

#         Compensation.computeCompensationMatrix!(fcs)
#     end

#     @test_nowarn begin fcs = loadFCControls("testdata/FlowRepository_FR-FCM-Z2SS_files")
#         channelnames = CSV.read("testdata/FlowRepository_FR-FCM-Z2SS_files/attachments/fcs_control.csv",DataFrame)
#         channelnames = Dict([i=>String(j) for (i,j) in eachrow(channelnames[:,["filename","dye"]])])
#         renameControl!(fcs,channelnames)

#         Compensation.computeCompensationMatrix!(fcs)

#         Compensation.compensate!(fcs)
#     end

#     @test_nowarn begin 
#         fcs = loadFCExperiment("testdata/FlowRepository_FR-FCM-Z2SS_files/Compensation Controls_APC-R700 Stained Control_007.fcs")
#         fcsControl = loadFCControls("testdata/FlowRepository_FR-FCM-Z2SS_files")
#         channelnames = CSV.read("testdata/FlowRepository_FR-FCM-Z2SS_files/attachments/fcs_control.csv",DataFrame)
#         channelnames = Dict([i=>String(j) for (i,j) in eachrow(channelnames[:,["filename","dye"]])])
#         renameControl!(fcsControl,channelnames)

#         Compensation.computeCompensationMatrix!(fcsControl)

#         Compensation.compensate!(fcs,control=fcsControl)
#     end

#     @test_nowarn begin 
#         fcs = loadFCExperiment("testdata/FlowRepository_FR-FCM-Z2SS_files/Compensation Controls_APC-R700 Stained Control_007.fcs")
#         fcsControl = loadFCControls("testdata/FlowRepository_FR-FCM-Z2SS_files")
#         channelnames = CSV.read("testdata/FlowRepository_FR-FCM-Z2SS_files/attachments/fcs_control.csv",DataFrame)
#         channelnames = Dict([i=>String(j) for (i,j) in eachrow(channelnames[:,["filename","dye"]])])
#         renameControl!(fcsControl,channelnames)

#         Compensation.computeCompensationMatrix!(fcsControl)
#         Compensation.assignCompensation!(fcs,control=fcsControl)
#         #Compensation.compensate!(fcs)
#     end

# end

# @testset "DimensionalityReduction" begin

#     fcs = FlowCytometryExperiment([rand(10,15);rand(10,15).+[10 0 0 0 0 0 0 0 0 0 0 0 0 0 0]])
#     fcs.var[:,"useChannels"] = [true,true,true,true,true,false,false,false,false,false,false,false,false,false,false]

#     @test_nowarn DimensionalityReduction.pca!(fcs)
#     @test begin DimensionalityReduction.pca!(fcs); size(fcs.obsm["pca"])==(20,15) end
#     @test_nowarn DimensionalityReduction.pca!(fcs,key_used_channels="useChannels")
#     @test begin DimensionalityReduction.pca!(fcs,key_used_channels="useChannels"); size(fcs.obsm["pca"])==(20,5) end

#     @test_nowarn DimensionalityReduction.umap!(fcs)
#     @test begin DimensionalityReduction.umap!(fcs); size(fcs.obsm["umap"])==(20,2) end
#     @test_nowarn DimensionalityReduction.umap!(fcs,key_used_channels="useChannels")
#     @test begin DimensionalityReduction.umap!(fcs,key_used_channels="useChannels"); size(fcs.obsm["umap"])==(20,2) end
#     @test begin DimensionalityReduction.pca!(fcs); DimensionalityReduction.umap!(fcs, key_obsm="pca"); size(fcs.obsm["umap"])==(20,2) end

# end

# @testset "Clustering" begin
# 
    # fcs = FlowCytometryExperiment([rand(10,15);rand(10,15).+[10 0 0 0 0 0 0 0 0 0 0 0 0 0 0]])
    # fcs.var[:,"useChannels"] = [true,true,true,true,true,false,false,false,false,false,false,false,false,false,false]

    # @test_nowarn Clustering.kmeansTuning(fcs, k=[1,2,3,4,5])
    # @test begin l = Clustering.kmeansTuning(fcs, k=[1,2,3,4,5]); sum(l.>20) == 2 end
# 
    # @test_nowarn Clustering.kmeans!(fcs, k=2)
    # @test begin Clustering.kmeans!(fcs, k=2); all((fcs.obs[:,"kmeans"] .== fcs.obs[:,"kmeans"])[1:10]) && all((fcs.obs[:,"kmeans"] .== fcs.obs[end,"kmeans"])[11:end]) end
# 
    # @test_nowarn Clustering.agglomerative!(fcs, k=2)
    # @test begin Clustering.agglomerative!(fcs, k=2); all((fcs.obs[:,"kmeans"] .== fcs.obs[:,"kmeans"])[1:10]) && all((fcs.obs[:,"kmeans"] .== fcs.obs[end,"kmeans"])[11:end]) end
# 
    # @test_nowarn Clustering.gaussianMixtureEM!(fcs, k=2)
    # @test begin Clustering.gaussianMixtureEM!(fcs, k=2); all((fcs.obs[:,"gaussianMixture"] .== fcs.obs[:,"gaussianMixture"])[1:10]) && all((fcs.obs[:,"gaussianMixture"] .== fcs.obs[end,"gaussianMixture"])[11:end]) end
# 
# end

