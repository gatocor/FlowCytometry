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

    channels = [i for i in keys(flow.data)]
    var = DataFrame(:channel=>channels)
    
    X = zeros(length(flow.data[channels[1]]),length(channels))
    for (i,c) in enumerate(channels)
        X[:,i] .= flow.data[c]
    end
    ftcexp = FlowCytometryExperiment(X,var=var)

    ftcexp.uns["ExperimentInformation"] = flow.params

    return ftcexp
end

"""
    function loadFCControls(folder::String)

Function to upload all the files in a folder corresponding to control files.

**Arguments**
 - **folder::String** Address of folder where all the controls are.

**Returns**
 FlowCytometryControl object with all the controls uploaded.
"""
function loadFCControls(folder::String)

    controllist = readdir(folder)
    controllist = [i for i in controllist if occursin(".fcs",i)]

    object = FlowCytometryControl()

    for control in controllist
        object.controls[split(control,".fcs")[1]] = loadFCExperiment(string(folder,"/",control))
    end

    return object

end

"""
    function loadFCControls(dic::Dict{String,String})

Function to upload all the files as given by a dictionary of channels and the corresponding control files.

**Arguments**
 - **dic::Dict{String,String}** Disctionary with names of channels and file adress assotiated with them.

**Returns**
 FlowCytometryControl object with all the controls uploaded.
"""
function loadFCControls(dic::Dict{String,String})

    object = FlowCytometryControl()

    for (control,file) in pairs(dic)
        object.controls[control] = loadFCExperiment(file)
    end

    return object

end