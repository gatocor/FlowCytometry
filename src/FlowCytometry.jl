module FlowCytometry

using DataFrames, FileIO, FCSFiles

export FlowCytometryExperiment, FlowCytometryControl, FlowCytometryGate
export loadFCExperiment, isInsideGate
export removeCells, removeCells!, removeChannels, removeChannels!

include("./structures.jl")
include("./inPolygon.jl")
include("./IO.jl")
include("./clustering.jl")


end # module
