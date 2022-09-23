module FlowCytometry

using DataFrames, FileIO, FCSFiles, MLJ, UMAP

export FlowCytometryExperiment, FlowCytometryControl, FlowCytometryGate
export loadFCExperiment, isInsideGate
export removeCells, removeCells!, removeChannels, removeChannels!
export Clustering, DimensionalityReduction

include("./structures.jl")
include("./inPolygon.jl")
include("./IO.jl")
include("./clustering.jl")
include("./dimensionalityReduction.jl")

end # module
