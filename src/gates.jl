module Gating

    using Dash, DataFrames, JSON, PlotlyJS, FlowCytometry, Statistics

"""
    function manualGating!(fcs::FlowCytometryExperiment;debug=false)

Launch an interactive page accessible at localhost to define gates manually.

**Arguments**

 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment object where to define the gates.

**Returns**
 Nothing. Adds or removes the gates to the object in place.
"""    
    function manualGating!(fcs::FlowCytometryExperiment;debug=false)

        gates = fcs.gates

        app = dash()
    
        dropdown_options = [
            Dict("label" => i, "value" => i) for i in sort(fcs.channels)
        ]
        scale_options = [
            Dict("label" => i, "value" => i) for i in ["lin","log"]
        ]
    
        function generate_table(gates)
    
            return [Dict("Name"=>r,"ChannelX"=>gates[r].channels[1],"ChannelY"=>gates[r].channels[2]) for r = keys(gates)]
    
        end
    
        p1 = Plot(fcs[fcs.channels[1]], fcs[fcs.channels[2]], marker_size=8)
    
        app.layout = html_div() do
            html_h1("Manual Gating",
                style = Dict("color" => "#00000", "textAlign" => "center")),
            html_div("A simple to tool to define manual gates over FlowCytometryExperiment objects."),
            html_ol("1. Select channels to visualize. You can choose scales and threshold (as negative values in log scale will not show otherwise)."),
            html_ol("2. Select the area of gate using `Box select` or `Lasso Select` tools in the plot."),
            html_ol("3. Write a name and press `add gate`."),
            html_ol("4. When done to add gates, press `Close`."),
            html_div(style = Dict("margin-top" => "30px")),
            html_div(style = Dict("columnCount" => 2)) do
                dcc_graph(
                    id = "scatter",
                    figure = p1
                ),
                html_div(html_th("Visualizing gates"),style = Dict("margin-bottom" => "10px")),
                html_div(style = Dict("columnCount" => 2)) do
                    html_div("Gate X"),
                    dcc_dropdown(
                        id = "inputX",
                        options = dropdown_options,
                        value = fcs.channels[1],
                        multi = false,
                    ),
                    dcc_dropdown(
                        id = "scaleX",
                        options = scale_options,
                        value = "lin",
                        multi = false,
                    ),
                    html_label("Threshold X "),
                    dcc_input(id = "thresholdX", value = 0, type = "number"),
                    html_div("Gate Y"),
                    dcc_dropdown(
                        id = "inputY",
                        options = dropdown_options,
                        value = fcs.channels[2],
                        multi = false,
                    ),
                    dcc_dropdown(
                        id = "scaleY",
                        options = scale_options,
                        value = "lin",
                        multi = false,
                    ),
                    html_label("Threshold Y "),
                    dcc_input(id = "thresholdY", value = 0, type = "number")
                end,
                html_div(style = Dict("margin-top" => "50px")),
                html_div(style = Dict("columnCount" => 2)) do
                    html_div(html_th("Add gate")),
                    dcc_input(id = "gate-input-name", value = "Add name of gate here...", type = "text"),
                    html_button(id = "submit-gate", children = "add gate", n_clicks = 0),
                    html_div(html_th("Remove gate"),
                            style=Dict("margin-top"=>"15px")),
                    dcc_dropdown(
                        id = "gate-remove-name",
                        options = [Dict("label" => i, "value" => i) for i in keys(gates)],
                        value = "",
                        multi = false,
                    ),
                    html_button(id = "delete-gate", children = "remove gate", n_clicks = 0)
                end,
                html_div(id="relayout-data"),
                html_div(style = Dict("margin-top" => "50px")),
                html_label(html_th("Defined gates")),
                dash_datatable(
                    id = "gates",
                    data = [],
                    columns=[Dict("name" =>c, "id" => c) for c in ["Name","ChannelX","ChannelY"]]
                )
            end
        end
    
        callback!(app, Output("scatter", "figure"), 
                        Input("inputX", "value"), 
                        Input("inputY", "value"), 
                        Input("scaleX","value"), 
                        Input("scaleY","value"),
                        Input("thresholdX","value"), 
                        Input("thresholdY","value")) do x,y,xscale,yscale,xthreshold,ythreshold
            fig = plot(
                scatter(x=fcs[string(x)].+xthreshold, y=fcs[string(y)].+ythreshold, mode="markers", marker_size=8, yaxis_title="Hola"),
                Layout(xaxis_title=x,yaxis_title=y,height=800,width=800,xaxis_type=xscale,yaxis_type=yscale)
            )
            return fig
        end
    
        callback!(app, Output("gate-input-name","value"), 
                        Input("submit-gate", "n_clicks"), 
                        State("gate-input-name","value"), 
                        State("inputX", "value"), 
                        State("inputY", "value"), 
                        State("scaleX", "value"), 
                        State("scaleY", "value"), 
                        State("thresholdX", "value"), 
                        State("thresholdY", "value"), 
                        State("scatter", "selectedData")) do button,name,x,y,xscale,yscale,xthreshold,ythreshold,box
            if name in keys(gates)
                return "Gate name already exisiting in gates."
            elseif "Add name of gate here..." == name
                return "Add name of gate here..."
            elseif "Gate name already exisiting in gates." == name
                return "Gate name already exisiting in gates."
            else
                if "range" in keys(box)
                    gates[name] = FlowCytometryGate((x,y),[(box[:range][:x][1]-xthreshold,box[:range][:y][1]-ythreshold),
                                                            (box[:range][:x][1]-xthreshold,box[:range][:y][2]-ythreshold),
                                                            (box[:range][:x][2]-xthreshold,box[:range][:y][2]-ythreshold),
                                                            (box[:range][:x][2]-xthreshold,box[:range][:y][1]-ythreshold)])
                else
                    gates[name] = FlowCytometryGate((x,y),[(i,j) for (i,j) in zip(box[:lassoPoints][:x].-xthreshold,box[:lassoPoints][:y].-ythreshold)])
                end
                return "Add name of gate here..."
            end
        end
    
        callback!(app, Output("delete-gate", "n_clicks"), 
                        Input("delete-gate", "n_clicks"), 
                        State("gate-remove-name","value")) do button,name
            delete!(gates,name)
            return 0
        end
    
        callback!(app, Output("gates","data"), 
                        Output("gate-remove-name","options"), 
                        Input("gate-input-name", "value"), 
                        Input("delete-gate", "n_clicks")) do addbutton,delbutton
            return generate_table(gates), [Dict("label" => i, "value" => i) for i in keys(gates)]
        end
    
        run_server(app, "0.0.0.0", debug=debug)

        return
    end

    function convexHull(x,y)

        #Find initial points p
        p = findall(y.==findmin(y)[1])
        if length(p)>1 #remove ties
            p = findmin(x[p])[2]
        else
            p = p[1]
        end
    
        #Compute angles
        angle(p0,p1) = (p1[1]-p0[1])/sqrt((p1[1]-p0[1])^2+(p1[2]-p0[2])^2)
        angles = zeros(length(x))
        for i in 1:length(x)
            angles[i] = angle((x[p],y[p]),(x[i],y[i]))
        end
    
        #Sort by visit
        order = sortperm(-angles)
    
        #Make hull
        turnRight(p0,p1,p2) = ((p1[1]-p0[1])*(p2[2]-p1[2])-(p2[1]-p1[1])*(p1[2]-p0[2])) < 0
        p0 = (x[p],y[p])
        p1 = x[order[1]],y[order[1]]
        hull = [p0,p1]
        for i in order[2:end-1]
            p2 = (x[i],y[i])
            while turnRight(p0,p1,p2) #Turn right, remove previous point 
                pop!(hull)
                p1 = hull[end]
                p0 = hull[end-1]
            end
            push!(hull,p2)
            p0 = p1
            p1 = p2
        end
    
        push!(hull,hull[1])
    
        return hull
    end

"""
    function automaticQC!(fcs::FlowCytometryExperiment;
        channel1="FSC-A",channel2="SSC-A",
        trim::Tuple{<:Real,<:Real}=(0.01,0.99),
        maxtrim::Real=.05,
        densityBandwidth::Tuple{<:Real,<:Real}=(.5,.3),
        finalBandwidth=.1,
        heightFromMax=.3,
        subsample=2000,
        keyAdded="automaticQC")

Function that automatically detects the main peak corresponding to viable cells in a FCS/SSC scatterplot and adds a gate to the object. 
The method follows the heuristics of [Roca et al](https://www.nature.com/articles/s41467-021-23126-8) with some small adaptations.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment where to compute the gate

**Keyword arguments**
 - **channel1="FSC-A"** Channel name in channels corrresponding to the forwad scatter
 - **channel2="SSC-A"** Channel name in channels corrresponding to the forwad scatter
 - **trim::Tuple{<:Real,<:Real}=(0.01,0.99)** Percentile range in the channel ranges to use in the analysis before finding maximums.
 - **maxtrim::Real=.05** Minimum in the range channels where to accept maximum points. This avoids the peak present in some datasets at lower expressions (debris).
 - **densityBandwidth::Tuple{<:Real,<:Real}=(.5,.3)** Relative bandwidth to smooth the density plots.
 - **finalBandwidth=.1** Relative bandwidth to smooth the density plots in the final step.
 - **heightFromMax=.3** In the last step, heigh at global maximum where to draw the contour of the QC region. Points with more probability that that will be used to make the QC gate.
 - **subsample=2000** Subsample of points used for density steps of the system. Set to 'nothing' for using all the points (this can be extremely slow), with relative low numbers of points the gate output are already very good.
 - **keyAdded="automaticQC"** Key added to uns and gate name.

**Returns**
Nothing. A gate of automatic control is added to gates.

"""
    function automaticQC!(fcs::FlowCytometryExperiment;
                        channel1="FSC-A",channel2="SSC-A",
                        trim::Tuple=(0.01,0.99),
                        maxtrim::Real=.05,
                        densityBandwidth::Tuple=(.5,.3),
                        finalBandwidth=.1,
                        heightFromMax=.3,
                        subsample=2000,
                        keyAdded="automaticQC")

        #Check channels are present
        if !(channel1 in fcs.channels)
            error(channel1, " is no in channels.")
        elseif !(channel2 in fcs.channels)
            error(channel2, " is no in channels.")
        end

        #Make dictionary step results
        dic = Dict{Any,Any}()
    
        #Extract channels
        x = fcs[channel1]
        y = fcs[channel2]
    
        #Remove extrema
        xlims = sort(x)[[Int(round(trim[1]*length(x))),Int(round(trim[2]*length(x)))]]
        ylims = sort(y)[[Int(round(trim[1]*length(y))),Int(round(trim[2]*length(y)))]]
        keep = (x .>= xlims[1]) .& (x .<= xlims[2]) .& (y .>= ylims[1]) .& (y .<= ylims[2])
        xMaxTrim = (maximum(x[keep])-minimum(x[keep]))*maxtrim+minimum(x[keep])
        yMaxTrim = (maximum(y[keep])-minimum(y[keep]))*maxtrim+minimum(y[keep])
    
        #Scale bandwidths
        densityBandwidthScaled = [mean([std(x[keep]),std(y[keep])])*i for i in densityBandwidth]
        hullBandwidthScaled = mean([std(x[keep]),std(y[keep])])*finalBandwidth
        #densityBandwidthScaled = [(maximum(x[keep])-minimum(x[keep]))*i for i in densityBandwidth]
        #hullBandwidthScaled = (maximum(x[keep])-minimum(x[keep]))*finalBandwidth

        #Subset for the rest of the computations
        if subsample === nothing
            subset = 1:length(x)
        else
            if subsample < length(x)
                subset = unique(rand(1:length(x),subsample))
            else
                subset = 1:length(x)
            end
        end
        x = x[subset]
        y = y[subset]
    
        #Initialize arrays
        active = fill(false,length(x))
        keep = (x .>= xlims[1]) .& (x .<= xlims[2]) .& (y .>= ylims[1]) .& (y .<= ylims[2])
        active[keep] .= true
        rectKeep = fill(false,length(x))
        rectKeepOld = fill(true,length(x))
        keepMaxTrim = (x .>= xMaxTrim[1]) .& (y .>= yMaxTrim[1]) 
        pointsFromMaximum = nothing
    
        #Add to dictionary
        dic["densityBandwidth"] = densityBandwidth
        dic["finalBandwidth"] = finalBandwidth
        dic["densityBandwidthScaled"] = densityBandwidthScaled
        dic["hullBandwidthScaled"] = hullBandwidthScaled
        dic["iterations"] = length(densityBandwidth)
        dic["order"] = sortperm(subset)
        fcs.obs[!,string(keyAdded,"_subset")] .= false
        fcs.obs[subset,string(keyAdded,"_subset")] .= true
        fcs.obs[!,string(keyAdded,"_trim")] .= false
        #keep .= keepMaxTrim .& keep
        dic["trim"] = [[xlims[1],xlims[1],xlims[2],xlims[2],xlims[1]],[ylims[1],ylims[2],ylims[2],ylims[1],ylims[1]]]
        dic["trim_centers"] = [[xMaxTrim,xMaxTrim,xlims[2],xlims[2],xMaxTrim],[yMaxTrim,ylims[2],ylims[2],yMaxTrim,yMaxTrim]]
        
        dic["excluded_cells"] = keepMaxTrim
        density = zeros(length(x))
        count = 1
        for (k,band) in enumerate(densityBandwidthScaled)

            density .= 0
        
            #Weight all points
            @inbounds for i in 1:length(x)
                if active[i]
                    for j in 1:length(x)
                        if active[j]
                            density[j] += exp(-((x[i]-x[j])^2+(y[i]-y[j])^2)/(2*band^2))
                        end
                    end
                end
            end
        
            #Find local maximum
            d = 0
            maxMaximum = 0
            @inbounds for i in 1:length(x)
                if density[i] > d && keepMaxTrim[i] && k == 1
                    maxMaximum = i
                    d = density[i]
                elseif  density[i] > d && k != 1
                    maxMaximum = i
                    d = density[i]                
                end
                @inbounds for j in i+1:1:length(x)
                    if ((x[i]-x[j])^2+(y[i]-y[j])^2) < (band)^2
                        if density[i] < density[j]
                            active[i] = false
                        else
                            active[j] = false
                       end
                    end
                end
            end
        
            #Find global maximum excluding low maximums
            maximums = findall(active)
            othermaximums = [i for i in maximums if i != maxMaximum]
        
            #Tesselation
            pointsFromMaximum = fill(true,size(x)[1])
            #pointsFromMaximum[.!(keep)] .= false
            @inbounds for i in 1:length(x)
                    dMax = (x[i]-x[maxMaximum])^2+(y[i]-y[maxMaximum])^2
                    @inbounds for j in othermaximums
                    if ((x[i]-x[j])^2+(y[i]-y[j])^2) < dMax
                        pointsFromMaximum[i] = false
                    end
                end
            end
        
            #Find rectangle
            rectXmin = x[maxMaximum]-3*std(x[pointsFromMaximum])
            rectXmax = x[maxMaximum]+3*std(x[pointsFromMaximum])
            rectYmin = y[maxMaximum]-3*std(y[pointsFromMaximum])
            rectYmax = y[maxMaximum]+3*std(y[pointsFromMaximum])
            rectKeep = (x .>= rectXmin) .& (x .<= rectXmax) .& (y .>= rectYmin) .& (y .<= rectYmax)
            if k != length(densityBandwidth)
                active = rectKeep#pointsFromMaximum .& rectKeep .& keepMaxTrim .& keep .& rectKeepOld
            else
                active = pointsFromMaximum .& rectKeepOld .& keep #In the last approximation just keep the tesselation and box
            end
            rectKeepOld .= rectKeep

            dic[count] = Dict(["hull"=>convexHull(x[active],y[active]),
                                "density"=>copy(density),
                                "maximums"=>(x[othermaximums],y[othermaximums]),
                                "maximumGlobal"=>(x[maxMaximum],y[maxMaximum]),
                                "box"=>[[rectXmin,rectXmin,rectXmax,rectXmax,rectXmin],[rectYmin,rectYmax,rectYmax,rectYmin,rectYmin]]])
            count += 1
        end
    
        density2 = zeros(length(x))
        #Weight all points
        @inbounds for i in 1:length(x)
            if active[i]
                @inbounds for j in 1:length(x)
                    if active[j]
                        density2[j] += exp(-((x[i]-x[j])^2+(y[i]-y[j])^2)/(2*hullBandwidthScaled^2))
                    end
                end
            end
        end

        #Find points in the limit from maximum
        pointsInHull = (maximum(density2)*heightFromMax .<= density2) .& active
    
        hull = convexHull(x[pointsInHull],y[pointsInHull])
        fcs.gates[keyAdded] = FlowCytometryGate((channel1,channel2),hull)
    
        dic["density"] = density2
    
        fcs.uns[keyAdded] = dic
        
    end

"""
    function filterByGate!(fcs::FlowCytometryExperiment,gate::FlowCytometryGate,addKey::String)

Function that checks which cells are inside a gate.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment where to compute the gate
 - **gate::FlowCytometryGate** Gate to filter the cells.
 - **addKey::String** Kay to be added in .obs

**Returns**
Nothing. A columns with `addKey` is added to .obs.
"""
    function filterByGate!(fcs::FlowCytometryExperiment,gate::FlowCytometryGate,addKey)
       
        X = [fcs[gate.channels[1]] fcs[gate.channels[2]]]
        inside = isInsideGate(X,gate)
        
        fcs.obs[!,addKey] = inside
        
        return
    end

"""
    function filterByGate!(fcs::FlowCytometryExperiment,gates::Vector{FlowCytometryGate},addKey::String)

Function that checks which cells are inside a set of gates.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment where to compute the gate
 - **gates::Vector{FlowCytometryGate}** Set of gates to filter the cells.
 - **addKey::String** Kay to be added in .obs

**Returns**
Nothing. A columns with `addKey` is added to .obs.
"""
    function filterByGate!(fcs::FlowCytometryExperiment,gates::Vector{FlowCytometryGate},addKey)
       
        insideTotal = fill(true,size(fcs.X)[1])
        for gate in gates
            X = [fcs[gate.channels[1]] fcs[gate.channels[2]]]
            inside = isInsideGate(X,gate)
            insideTotal .&= inside
        end

        fcs.obs[!,addKey] = inside
        
        return
    end

"""
    function filterByGate!(fcs::FlowCytometryExperiment,gateName::String,addKey=string(gateName,"_gate"))

Function that checks which cells are inside a gate.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment where to compute the gate
 - **gateName::String** Name of gate inside FlowCytometryExperiment.
 - **addKey::String=string(gateName,"_gate")** Key to be added in .obs

**Returns**
Nothing. A columns with `addKey` is added to .obs.

"""
    function filterByGate!(fcs::FlowCytometryExperiment,gateName::String,addKey=string(gateName,"_gate"))
       
        X = [fcs[fcs.gates[gateName].channels[1]] fcs[fcs.gates[gateName].channels[2]]]
        inside = isInsideGate(X,fcs.gates[gateName])
        
        fcs.obs[!,addKey] = inside
        
        return
    end

"""
    function filterByGate!(fcs::FlowCytometryExperiment,gateNames::Vector{String},addKey)

Function that checks which cells are inside a set of gates.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment where to compute the gate
 - **gateNames::Vector{String}** Name of gates inside FlowCytometryExperiment.
 - **addKey::String=string(gateName,"_gate")** Key to be added in .obs

**Returns**
Nothing. A columns with `addKey` is added to .obs.

"""
    function filterByGate!(fcs::FlowCytometryExperiment,gateNames::Vector{String},addKey)
       
        insideTotal = fill(true,size(fcs.X)[1])
        for gateName in gateNames
            X = [fcs[fcs.gates[gateName].channels[1]] fcs[fcs.gates[gateName].channels[2]]]
            inside = isInsideGate(X,fcs.gates[gateName])
            insideTotal .&= inside
        end
        
        fcs.obs[!,addKey] = insideTotal
        
        return
    end

end