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