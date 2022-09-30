module FlowCytometry

using DataFrames, FileIO, FCSFiles, MLJ, UMAP, HDF5

export FlowCytometryExperiment, FlowCytometryControl, FlowCytometryGate
export renameControl!, checkControlNames
export loadFCExperiment, loadFCControls, saveH5fcs, loadH5fcs
export isInsideGate
export removeCells, removeCells!, removeChannels, removeChannels!
export Gating, Compensation, DimensionalityReduction, Clustering

include("./structures.jl")
include("./inPolygon.jl")
include("./IO.jl")
include("./gates.jl")
include("./compensation.jl")
include("./dimensionalityReduction.jl")
include("./clustering.jl")

end # module
