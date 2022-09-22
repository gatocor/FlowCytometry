"""
    struct FlowCytometryGate

Structure that contains information from a gating process of the data.

**Elements:**
 - **channels::Tuple{String,String}** Channel names of the gate.
 - **polygon::Vector{Tuple{Real,Real}}** Points of poligon that defines the gate.
"""
struct FlowCytometryGate
    channels::Tuple{String,String}
    polygon::Vector{Tuple{Real,Real}}
end

"""
    mutable struct FlowCytometryControl
    
Structure that stores the data coming from control stainings for correcting for spillover.

**Elements**:
 - **var::DataFrame** Dataframe including the metainformation of the channels. All the controls should have the same channels measures.
 - **obs::Dict{String,DataFrame}** Dictionary of Dataframes containing the metainformation from the cells in each channel.
 - **X::Dict{String,Matrix{AbstractFloat}}** Dictionary of matrices from the control measures.

"""
mutable struct FlowCytometryControl
    var::DataFrame
    obs::Dict{String,DataFrame}
    X::Dict{String,Matrix{AbstractFloat}}
end

"""
    mutable struct FlowCytometryExperiment
    
Structure containing a Flow Cytometry experiment and all the processes applyied to it.

**Elements:**   
 - **X::Matrix{AbstractFloat}** Matrix of Cells X Channel of the experiment.
 - **obs::DataFrame** Dataframe with all the metainformation of the cells
 - **var::DataFrame** Dataframe with all the metainformation of the channels
 - **obsm::Dict{String,Matrix{AbstractFloat}}** Dictionary containing transformed matrices of the original data.
 - **layers::Dict{String,Matrix{AbstractFloat}}** Dictionary containing Cells X Channel matrices of data that are required to control (e.g. Raw matrix).
 - **controls::Dict{FlowCytometryControl}** Control object with the control staining to compute spillover matrices.
 - **gates::Dict{String,FlowCytometryGate}** List of Gate and Gate set objects
 - **uns::Dict{String,Any}** Dictionary contining all the metainformation of cells

"""
mutable struct FlowCytometryExperiment
    X::Matrix{<:AbstractFloat}
    obs::DataFrame
    var::DataFrame
    obsm::Dict{String,Matrix{<:AbstractFloat}}
    layers::Dict{String,Matrix{<:AbstractFloat}}
    controls::Union{Nothing,FlowCytometryControl}
    gates::Dict{String,FlowCytometryGate}
    uns::Dict{String,Any}

    function FlowCytometryExperiment(X::Matrix{<:AbstractFloat};obs::DataFrame=DataFrame(),var::DataFrame=DataFrame())
        
        s = size(X)

        if isempty(obs)
            obs = DataFrame(:cell=>1:s[1])
        elseif size(obs)[1] != s[1]
            error("obs dataframe must have the same number of rows than rows in X.")
        end

        if isempty(var)
            var = DataFrame(:channel=>1:s[2])
        elseif size(var)[1] != s[2]
            error("var dataframe must have the same number of rows than columns in X.")
        end

        return new(X,obs,var,
            Dict{String,Matrix{AbstractFloat}}(),
            Dict{String,Matrix{AbstractFloat}}(),
            nothing,
            Dict{String,FlowCytometryGate}(),
            Dict{String,Any}()
            )
    end
end

function Base.setproperty!(fcs::FlowCytometryExperiment, symbol::Symbol, data::DataFrame)
    if symbol == :obs 
        s = size(fcs.obs)
        s2 = size(data)
        if s[1] != s2[1]
            error("Assigned dataframe must have the same number of cells. The dataset has ", s[1], " cells and ", s2[1], " are trying to be assigned.")
        else
            setfield!(fcs,:obs,data)
        end
    elseif symbol == :var
        s = size(fcs.var)
        s2 = size(data)
        if s[1] != s2[1]
            error("Assigned dataframe must have the same number of channels. The dataset has ", s[1], " channels and ", s2[1], " are trying to be assigned.")
        else
            setfield!(fcs,:var,data)
        end
    else
        setfield!(adata,symbol,data)
    end
end

function Base.setproperty!(fcs::FlowCytometryExperiment, symbol::Symbol, data::Matrix)
    if symbol == :X 
        s1 = size(fcs.obs)[1]
        s2 = size(fcs.var)[1]
        s = size(data)
        if (s1,s2) != s
            error("Assigned matrix must have the same size as original cells x channels. The dataset has size ", (s1,s2), " and size ", s, " is trying to be assigned.")
        else
            setfield!(fcs,:X,data)
        end
    else
        setfield!(adata,symbol,data)
    end
end

"""
    function Base.getindex(fcs::FlowCytometryExperiment, s::Union{Symbol,String,<:Real})

Access a channel directly by its .var.channel property.

e.g. 
    >>> fct.var.channels
        ["FCT-Scatter-A","FCT-Scatter-B"]
    >>> fct["FCT-Scatter-A"]
        [0.1,20.33,0.45,0.28]
"""
function Base.getindex(fcs::FlowCytometryExperiment, s::Union{Symbol,String,<:Real})
    pos = findfirst(fcs.var.channel .== s)
    if pos === nothing
        error(s, " has not been found in var.channels.")
    else
        return fcs.X[:,pos]
    end
end

"""
    function removeCells(fcs::FlowCytometryExperiment, s::Vector{Bool})

Remove cells from the dataset and return a copy.

**Arguments**:
 - **fcs::FlowCytometryExperiment** Dataset to be pruned
 - **s::Vector{Bool}** Array of cells to be removed.

**Returns**
FlowCytometryExperiment with the cells not prunned.
"""
function removeCells(fcs::FlowCytometryExperiment, s::Vector{Bool})

    fcsNew = FlowCytometryExperiment(fcs.X[s,:])
    fcsNew.obs = deepcopy(fcs.obs[s,:])
    fcsNew.var = deepcopy(fcs.var)
    obsmkeys = [i for i in keys(fcs.obsm)]
    for i in obsmkeys
        fcsNew.obsm[i] = fcs.obsm[i][s,:]  
    end
    fcsNew.layers = deepcopy(fcs.layers)
    fcsNew.controls = deepcopy(fcs.controls)
    fcsNew.gates = deepcopy(fcs.gates)
    fcsNew.uns = deepcopy(fcs.uns)

    return fcsNew
end

"""
    function removeCells!(fcs::FlowCytometryExperiment, s::Vector{Bool})

Remove cells from the dataset.

**Arguments**:
 - **fcs::FlowCytometryExperiment** Dataset to be pruned
 - **s::Vector{Bool}** Array of cells to be removed.

**Returns**
Nothing
"""
function removeCells!(fcs::FlowCytometryExperiment, s::Vector{Bool})

    if size(fcs.X)[1] == length(s)
        setfield!(fcs,:X,fcs.X[s,:])
        setfield!(fcs,:obs,fcs.obs[s,:])
        obsmkeys = [i for i in keys(fcs.obsm)]
        for i in obsmkeys
            fcs.obsm[i] = fcs.obsm[i][s,:]  
        end
    end

    return
end

"""
    function removeChannels(fcs::FlowCytometryExperiment, s::Vector{Bool})

Remove channels from the dataset and return a copy.

**Arguments**:
 - **fcs::FlowCytometryExperiment** Dataset to be pruned
 - **s::Vector{Bool}** Array of channels to be removed.

**Returns**
FlowCytometryExperiment with the channels not removed.
"""
function removeChannels(fcs::FlowCytometryExperiment, s::Vector{Bool})

    fcsNew = FlowCytometryExperiment(fcs.X[:,s])
    fcsNew.obs = deepcopy(fcs.obs)
    fcsNew.var = deepcopy(fcs.var[s,:])
    fcsNew.obsm = deepcopy(fcs.obsm)
    fcsNew.layers = deepcopy(fcs.layers)
    fcsNew.controls = deepcopy(fcs.controls)
    fcsNew.gates = deepcopy(fcs.gates)
    fcsNew.uns = deepcopy(fcs.uns)

    return fcsNew
end

"""
    function removeChannels!(fcs::FlowCytometryExperiment, s::Vector{Bool})

Remove channels from the dataset.

**Arguments**:
 - **fcs::FlowCytometryExperiment** Dataset to be pruned
 - **s::Vector{Bool}** Array of channels to be removed.

**Returns**
FlowCytometryExperiment with the channels not removed.
"""
function removeChannels!(fcs::FlowCytometryExperiment, s::Vector{Bool})

    if size(fcs.X)[2] == length(s)
        setfield!(fcs,:X,fcs.X[:,s])
        setfield!(fcs,:var,fcs.var[s,:])
    end

    return
end