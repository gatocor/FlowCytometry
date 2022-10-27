module FCSPloting
    using Plots, FlowCytometry, Statistics

"""
    function plotQCSteps(fcs::FlowCytometryExperiment,key::String="automaticQC")

Function that returns the ploting steps from automaticQC.

**Arguments**

 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment object where to define the gates.
 - **key::String="automaticQC"** Key used to save the results from the automatic control object.

**Returns**
 Figure with the steps ploted.
"""
    function plotQCSteps(fcs::FlowCytometryExperiment,key::String="automaticQC")
        
        channel1,channel2 = fcs.gates[key].channels
        subset = fcs.obs[:,string(key,"_subset")]
        its = fcs.uns[key]["iterations"]
        order = fcs.uns[key]["order"]
        figs = []
        xMax = maximum(fcs.uns[key]["trim"][1])
        yMax = maximum(fcs.uns[key]["trim"][2])
        
        for it in 1:its
            f = scatter(fcs[channel1],fcs[channel2],xlims=(0,xMax),ylims=(0,yMax),markersize=3,markerstrokewidth=0,label="all data",colorbar=false)
            scatter!(f,fcs[channel1][subset],fcs[channel2][subset],zcolor=fcs.uns[key][it]["density"][order],xlims=(0,250E3),ylims=(0,250E3),markersize=3,markerstrokewidth=0,label="subsample with density")
            plot!(f,fcs.uns[key]["trim"][1],fcs.uns[key]["trim"][2],linewidth=2,color=:purple,style=:dash,label="trim")
            if it == 1 #Plot trimming of centers from first step
                plot!(f,fcs.uns[key]["trim_centers"][1],fcs.uns[key]["trim_centers"][2],linewidth=2,color=:pink,style=:dash,label="trim_centers")
            end
            if it != its #Plot limiting box for next step
                plot!(f,fcs.uns[key][it]["box"][1],fcs.uns[key][it]["box"][2],linewidth=2,color=:green,label="selection_box")
            end
            if it > 1 #Plot limiting box for next step
                plot!(f,fcs.uns[key][it-1]["box"][1],fcs.uns[key][it-1]["box"][2],linewidth=2,color=:green,label="cells used")
            end
            plot!(f,[i[1] for i in fcs.uns[key][it]["hull"]],[i[2] for i in fcs.uns[key][it]["hull"]],linewidth=2,color=:red,label="selected region")
            scatter!(f,fcs.uns[key][it]["maximums"][1],fcs.uns[key][it]["maximums"][2],color="red",xlims=(0,xMax),ylims=(0,yMax),markersize=5,markerstrokewidth=0,label="maximums")
            scatter!(f,[fcs.uns[key][it]["maximumGlobal"][1]],[fcs.uns[key][it]["maximumGlobal"][2]],color="red",xlims=(0,xMax),ylims=(0,yMax),markersize=10,marker=:square,markerstrokewidth=0,label="selected maximum")
            plot!(f,legend_position=:outerbottomright)
            
            push!(figs,f)
        end

        f = scatter(fcs[channel1],fcs[channel2],xlims=(0,xMax),ylims=(0,yMax),markersize=3,markerstrokewidth=0,colorbar=false,label="all data")
        scatter!(f,fcs[channel1][subset],fcs[channel2][subset],zcolor=fcs.uns[key]["density"][order],xlims=(0,xMax),ylims=(0,yMax),markersize=3,markerstrokewidth=0,label="subsample with density")
        plot!(f,[i[1] for i in fcs.gates[key].polygon],[i[2] for i in fcs.gates[key].polygon],linewidth=2,color=:red,label="final QC region")
        plot!(f,legend_position=:outerbottomright)
        
        push!(figs,f)
        
        return plot(figs...,layout=(1,its+1))
    end


"""
    function plotControlsCorrelations(fcsControl::FlowCytometryControl)

Function that displays a heatmap with the correlation between controls and channels. It can work as a summary of spillover between data.
   
**Arguments:**
 - **fcsControl::FlowCytometryControl** FlowCytometryControl object to compute the correlations between channels.

**Return**
 Heatmap figure with the control channels in rows and channels in columns.
"""
    function plotControlsCorrelations(fcsControl::FlowCytometryControl)
        
        X = zeros(length(fcsControl.controls),length(fcsControl.channels))
        for (i,m) in enumerate(pairs(fcsControl.controls))
            (channel2,fcs) = m
            for (j,channel1) in enumerate(fcsControl.channels)
               X[i,j] = cor([fcs[channel1] fcs[channel2]])[1,2]
            end
        end
        
        fig = heatmap(X,c=:redsblues)
        #,
        #    xticks=(1:length(fcsControl.channels),fcsControl.channels),
        #    yticks=(1:length([i for i in keys(fcsControl.controls)]),[i for i in keys(fcsControl.controls)]),tickfontrotation=.5)
        
        fig = plot(fig)
        
        return fig
    end


    """
    function plotControls(fcsControl::FlowCytometryControl,channels::Vector{Tuple{String,String}})

Function that makes an scatterplot of the tuples of channels among controls and channels to see the spillover. If computeCompensationMatrix! has been computed, it also plots the correction after the control.

**Arguments:**
 - **fcsControl::FlowCytometryControl** FlowCytometryControl object where to check the spillover.
 - **channels::Vector{Tuple{String,String}}** Vetor of Tuples of control sample and spilled channel that want to be seen.

**Returns**
List of figures with all the scatterplots between tuples. It can be ploted using plot(returned...,layout=(...))
"""
    function plotControls(fcsControl::FlowCytometryControl,channels::Vector{Tuple{String,String,String}})
    
        for tup in channels
            (cont,channel1,channel2) = tup
            if !(cont in keys(fcsControl.controls))
                error("Control channel name ", channel1 ," not found in controls. Closest candidates: ", [keys(fcsControl.controls)])
            elseif !(channel1 in fcsControl.controls[cont].channels)
                error("Control channel name ", channel1 ," not found in controls. Closest candidates: ", [i for i in fcsControl.controls[cont].channels if occursin(uppercase(channel1), uppercase(i))])
            elseif !(channel2 in fcsControl.controls[cont].channels)
                error("Control channel name ", channel2 ," not found in controls. Closest candidates: ", [i for i in fcsControl.controls[cont].channels if occursin(uppercase(channel2), uppercase(i))])
            end
        end
        
        figs = []
        for (cont,channel2,channel1) in channels
            fcs = fcsControl.controls[cont]
            subset = [i in fcs.channels for i in fcsControl.channels]
            xlims = (sort(fcs[channel1])[round(Int,length(fcs[channel1])*.01)],sort(fcs[channel1])[round(Int,length(fcs[channel1])*.99)])
            ylims = (sort(fcs[channel2])[round(Int,length(fcs[channel2])*.01)],sort(fcs[channel2])[round(Int,length(fcs[channel2])*.99)])
                
            if "compensation" in keys(fcsControl.uns) 
                if !fcsControl.uns["compensation"]["compensated"]
                    fig = scatter(fcs[channel1],fcs[channel2],markersize=1,xlabel=channel1,ylabel=channel2,markerstrokewidth=0,label="data") 
                    P = fcsControl.uns["compensation"]["spilloverMatrix"]
                    p1 = findfirst(fcs.channels .== channel1)
                    p2 = findfirst(fcs.channels .== channel2)
                    l = sort(copy(fcs.X[:,p2]))
                    l2 = sort(l).*P[p2,p1]
                    plot!(fig,l2,l,label="compensation trend")
                    
                    p1 = findfirst(fcs.channels.==channel1)
                    p2 = findfirst(fcs.channels.==channel2)
                    X = fcs.X*fcsControl.compensationMatrix.matrix[:,subset][subset,:]
                    scatter!(fig,X[:,p1],X[:,p2],markersize=1,markerstrokewidth=0,label="compensated")
    
                    xlims = (sort(X[:,p1])[round(Int,length(fcs[channel1])*.01)],sort(fcs[channel1])[round(Int,length(fcs[channel1])*.99)])
                    ylims = (sort(X[:,p2])[round(Int,length(fcs[channel2])*.01)],sort(fcs[channel2])[round(Int,length(fcs[channel2])*.99)])
                    plot!(xlims=xlims,ylims=ylims)
                else
                    fig = scatter(fcs[channel1],fcs[channel2],markersize=1,xlabel=channel1,ylabel=channel2,markerstrokewidth=0,label="compensated") 
                    P = fcsControl.uns["compensation"]["spilloverMatrix"]
                    p1 = findfirst(fcs.channels .== channel1)
                    p2 = findfirst(fcs.channels .== channel2)
                    
                    p1 = findfirst(fcs.channels.==channel1)
                    p2 = findfirst(fcs.channels.==channel2)
                    X = fcs.X*inv(fcsControl.compensationMatrix.matrix)[:,subset][subset,:]
                    scatter!(fig,X[:,p1],X[:,p2],markersize=1,markerstrokewidth=0,label="data")                
    
                    l = sort(copy(X[:,p2]))
                    l2 = sort(l).*P[p2,p1]
                    plot!(fig,l2,l,label="compensation trend")
    
                    xlims = (sort(fcs[channel1])[round(Int,length(fcs[channel1])*.01)],sort(X[:,p1])[round(Int,length(fcs[channel1])*.99)])
                    ylims = (sort(fcs[channel2])[round(Int,length(fcs[channel2])*.01)],sort(X[:,p2])[round(Int,length(fcs[channel2])*.99)])
                    plot!(xlims=xlims,ylims=ylims)
                end
            else
                fig = scatter(fcs[channel1],fcs[channel2],markersize=1,xlabel=channel1,ylabel=channel2,markerstrokewidth=0,label="data") 
            end
                
            
            push!(figs,fig)
        end
            
        return figs
    end
end