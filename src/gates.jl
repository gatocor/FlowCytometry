module Gating

    using Dash, DataFrames, JSON, PlotlyJS, FlowCytometry

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

    function automaticQC(fcs::FlowCytometryExperiment;channel1="FSC-A",channel2="SSC-A",trim::Tuple{<:Real,<:Real}=(0.01,0.99),maxtrim::Real=.05,bandwidth::Tuple{<:Real,<:Real}=(.05,.05),hull=.05,pBorder=.33)

        #Check channels are present
        if !(channel1 in fcs.channels)
            error(channel1, " is no in channels.")
        elseif !(channel2 in fcs.channels)
            error(channel2, " is no in channels.")
        end

        x = fcs[channel1]
        y = fcs[channel2]
    
        #Remove extrema
        xlims = sort(x)[[Int(round(trim[1]*length(x))),Int(round(trim[2]*length(x)))]]
        ylims = sort(y)[[Int(round(trim[1]*length(y))),Int(round(trim[2]*length(y)))]]
        xMaxTrim = sort(x)[Int(round(trim[1]*length(x)))]
        yMaxTrim = sort(y)[Int(round(trim[1]*length(x)))]
        keepMaxTrim = (x .>= xMaxTrim[1]) .& (y .>= yMaxTrim[1]) 
        keep = (x .>= xlims[1]) .& (x .<= xlims[2]) .& (y .>= ylims[1]) .& (y .<= ylims[2])
    
        bandwidth = [sum([maximum(x[keep])-minimum(x[keep]),maximum(y[keep])-minimum(y[keep])])/2*i for i in bandwidth]
        hull = sum([maximum(x[keep])-minimum(x[keep]),maximum(y[keep])-minimum(y[keep])])/2*hull
        active = fill(false,length(x))
        active[keep] .= true
        rectKeep = fill(false,length(x))
    
        for (k,band) in enumerate(bandwidth)
            #Find center
            density = zeros(length(x))
            #Weight all points
            @inbounds for i in 1:length(x)
                if active[i]
                    @inbounds for j in 1:length(x)
                        if active[j]
                            density[i] += exp(-((x[i]-x[j])^2+(y[i]-y[j])^2)/(2*band))
                        end
                    end
                end
            end
            #Find local maximum
            @inbounds for i in 1:(length(x)-1)
                @inbounds for j in (i+1):length(x)
                    if ((x[i]-x[j])^2+(y[i]-y[j])^2) < (3*band)^2
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
            trimmedMaximums = [i for i in maximums if keepMaxTrim[i]]
            maxMaximum = maximums[findmax(density[trimmedMaximums])[2]]
            othermaximums = [i for i in maximums if i != maxMaximum]
            pointsFromMaximum = fill(true,size(fcs.obs)[1])
            pointsFromMaximum[.!(keep)] .= false
            @inbounds for i in 1:length(x)
                    dMax = ((x[i]-x[maxMaximum])^2+(y[i]-y[maxMaximum])^2)
                    @inbounds for j in othermaximums
                    if ((x[i]-x[j])^2+(y[i]-y[j])^2) < dMax
                        pointsFromMaximum[i] = false
                    end
                end
            end
        
            if k != length(bandwidth)
                #Find rectangle
                rectXmin = x[maxMaximum]-3*std(x[pointsFromMaximum])
                rectXmax = x[maxMaximum]+3*std(x[pointsFromMaximum])
                rectYmin = y[maxMaximum]-3*std(y[pointsFromMaximum])
                rectYmax = y[maxMaximum]+3*std(y[pointsFromMaximum])
                rectKeep = (x .>= rectXmin) .& (x .<= rectXmax) .& (y .>= rectYmin) .& (y .<= rectYmax)
                active = rectKeep .& keep
            else
                active = pointsFromMaximum .& rectKeep #In the last approximation just keep the tesselation and box
            end
        
#            return maximums, density, pointsFromMaximum, rectKeep
        end
    
        density = zeros(length(x))
        #Weight all points
        @inbounds for i in 1:length(x)
            if active[i]
                @inbounds for j in 1:length(x)
                    if active[j]
                        density[i] += exp(-((x[i]-x[j])^2+(y[i]-y[j])^2)/(2*hull))
                    end
                end
            end
        end
        
        pointsInHull = (maximum(density)*pBorder .<= density) #.& active
    
        p = findall(x[pointsInHull].==findmin(y)[1])
        if length(p)>1 #remove ties
            p = findmin(x[p])[2]
        end
    
        return pointsInHull
    
    end

end