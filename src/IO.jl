"""
    function loadFCExperiment(file::String)

Function that loads a ftc experiment into a FlowCytometryExperiment.

**Arguments**:
 - **file::String**: Path to .ftc file.

**Returns**:
 FlowCytometryExperiment with the data from the uploaded file.
"""
function loadFCExperiment(file::String)

    flow = FileIO.load(file)

    channels = sort([i for i in keys(flow.data)])
    var = DataFrame(:channel=>channels)
    
    X = zeros(length(flow.data[channels[1]]),length(channels))
    for (i,c) in enumerate(channels)
        X[:,i] .= flow.data[c]
    end
    ftcexp = FlowCytometryExperiment(X,var=var,channels=channels)

    ftcexp.uns["ExperimentInformation"] = flow.params

    return ftcexp
end

"""
    function loadFCControls(dic::Dict{String,String})

Function to upload all the files as given by a dictionary of control files and the corresponding channels.

**Arguments**
 - **dic::Dict{String,String}** Dictionary with `name_of_file`=>`name_of_controlled_channel`.

**Returns**
 FlowCytometryControl object with all the controls uploaded.
"""
function loadFCControls(dic::Dict{String,String})

    object = FlowCytometryControl()

    channels = []
    notCommonChannels = []
    for (file,control) in pairs(dic)
        fcs = loadFCExperiment(file)
        if isempty(channels)
            channels = copy(fcs.channels)
        else
            channels = sort(unique([channels;fcs.channels]))
            notCommonChannels = sort(unique([notCommonChannels;[i for i in channels if !(i in fcs.channels)]]))
        end
        
        if !(control in channels)
            error("Control channel name ", control ," not found in channels. Closest candidates: ", [i for i in channels if occursin(uppercase(control), uppercase(i))])
        end
        
        object.controls[control] = deepcopy(fcs)
    end
    
    object.channels = channels
    
    if length(notCommonChannels) > 0
       println("Not all control files have the same number of channels. These channels will not be compensated for controls that does not have them.\n Channels not present in all files: ",notCommonChannels) 
    end

    return object

end

function saveObject!(saveObject::Union{HDF5.File,HDF5.Group},data::DataFrame,name::String)
    create_group(saveObject,name)
    saveObject[name]["names"] = names(data)
    saveObject[name]["type"] = "DataFrame"
    create_group(saveObject[name],"data")
    for i in names(data)
        saveObject[name]["data"][i] = data[:,i]
    end

    return
end

function saveObject!(saveObject::Union{HDF5.File,HDF5.Group},data::Dict,name::String)
    create_group(saveObject,name)
    saveObject[name]["keys"] = [i[1] for i in pairs(data)]
    saveObject[name]["type"] = "Dict"
    create_group(saveObject[name],"data")
    for (key,d) in pairs(data)
        type = typeof(d)
        if type <: DataFrame || type <: Dict || type <: FlowCytometryGate || type <: Vector{<:Tuple} || type <: FlowCytometryExperiment
            saveObject!(saveObject[name]["data"],d,key)
        else
            saveObject[name]["data"][key] = d
        end
    end

    return
end

function saveObject!(saveObject::Union{HDF5.File,HDF5.Group},data::Vector{<:Tuple},name::String)

    create_group(saveObject,name)
    saveObject[name]["type"] = "Vector{Tuple}"
    saveObject[name]["x"] = [i[1] for i in data]
    saveObject[name]["y"] = [i[2] for i in data]

    return
end

function saveObject!(saveObject::Union{HDF5.File,HDF5.Group},data::FlowCytometryGate,name::String)

    create_group(saveObject,name)
    saveObject[name]["type"] = "FlowCytometryGate"
    saveObject[name]["channels"] = [i for i in data.channels]
    saveObject!(saveObject[name],data.polygon,"polygon")

    return
end

function saveObject!(saveObject::Union{HDF5.File,HDF5.Group},fcs::FlowCytometryExperiment,name::String)

    create_group(saveObject,name)

    saveObject[name]["type"] = "FlowCytometryExperiment"
    saveObject[name]["channels"] = fcs.channels
    saveObject[name]["X"] = fcs.X
    saveObject!(saveObject[name],fcs.obs,"obs")
    saveObject!(saveObject[name],fcs.var,"var")
    saveObject!(saveObject[name],fcs.obsm,"obsm")
    saveObject!(saveObject[name],fcs.layers,"layers")
    saveObject!(saveObject[name],fcs.gates,"gates")
    saveObject!(saveObject[name],fcs.uns,"uns")

    return
end

"""
    function saveH5fcs(fcs::FlowCytometryExperiment,file::String)

Function to save in .h5fcs a FlowCytometryExperiment.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment to save.
 - **file::String** File name where to store the data.

**Returns** Nothing
"""
function saveH5fcs(fcs::FlowCytometryExperiment,file::String)

    if !occursin(".h5fcs",file)
        file = string(file,".h5fcs")
    end
    saveObject = h5open(file,"w")

    saveObject["type"] = "FlowCytometryExperiment"
    saveObject["channels"] = fcs.channels
    saveObject["X"] = fcs.X
    saveObject!(saveObject,fcs.obs,"obs")
    saveObject!(saveObject,fcs.var,"var")
    saveObject!(saveObject,fcs.obsm,"obsm")
    saveObject!(saveObject,fcs.layers,"layers")
    saveObject!(saveObject,fcs.gates,"gates")
    saveObject!(saveObject,fcs.uns,"uns")

    close(saveObject)

    return

end

"""
    function saveH5fcs(fcs::FlowCytometryControl,file::String)

Function to save in .h5fcs a FlowCytometryControl.

**Arguments**
 - **fcs::FlowCytometryControl** FlowCytometryControl to save.
 - **file::String** File name where to store the data.

**Returns** Nothing
"""
function saveH5fcs(fcs::FlowCytometryControl,file::String)

    if !occursin(".h5fcs",file)
        file = string(file,".h5fcs")
    end
    saveObject = h5open(file,"w")

    saveObject["type"] = "FlowCytometryControl"
    saveObject["channels"] = fcs.channels
    saveObject!(saveObject,fcs.controls,"controls")
    if fcs.compensationMatrix === nothing
        saveObject["compensationMatrix"] = "nothing"
    else
        saveObject["compensationMatrix"] = fcs.compensationMatrix
    end
    saveObject!(saveObject,fcs.uns,"uns")

    close(saveObject)

    return

end

function loadObject(object::Union{HDF5.File,HDF5.Group})

    type = read(object["type"])
    if type == "DataFrame"

        d = DataFrame()
        for i in read(object["names"])
            d[:,i] = read(object["data"][i])
        end

        return d

    elseif type == "Dict"

        d = Dict()
        for i in read(object["keys"])
            if typeof(object["data"][i]) <: HDF5.Dataset 
                d[i] = read(object["data"][i])
            else
                d[i] = loadObject(object["data"][i])
            end
        end

        return d

    elseif type == "Vector{Tuple}"

        return [i for i in zip(read(object["x"]),read(object["y"]))]

    elseif type == "FlowCytometryGate"

        return FlowCytometryGate((read(object["channels"])[1],read(object["channels"])[2]),loadObject(object["polygon"]))

    elseif type == "FlowCytometryExperiment"

        fcs = FlowCytometryExperiment(read(object["X"]))
        fcs.channels = read(object["channels"])
        fcs.obs = loadObject(object["obs"])
        fcs.var = loadObject(object["var"])
        fcs.obsm = loadObject(object["obsm"])
        fcs.layers = loadObject(object["layers"])
        fcs.gates = loadObject(object["gates"])
        fcs.uns = loadObject(object["uns"])

        return fcs

    elseif type == "FlowCytometryControl"

        fcs = FlowCytometryControl()
        fcs.channels = read(object["channels"])
        fcs.controls = loadObject(object["controls"])
        compensation = read(object["compensationMatrix"])
        if compensation == "nothing"
            saveObject["compensationMatrix"] = nothing
        else
            saveObject["compensationMatrix"] = compensation
        end
        fcs.uns = loadObject(object["uns"])
    
        return fcs

    end
end

"""
    function loadH5fcs(file::String)

Function to load .h5fcs files into data.

**Arguments**
 - **file::String** File name of the datafile.

**Returns** FlowCytometryControl or FlowCytometryExperiment loaded from file.
"""
function loadH5fcs(file::String)

    object = h5open(file,"r")

    objectLoaded = loadObject(object)

    close(object)

    return objectLoaded

end
