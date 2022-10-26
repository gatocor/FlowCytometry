module FCSPloting
    using Plots, FlowCytometry

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
end