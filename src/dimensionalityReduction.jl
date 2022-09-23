module DimensionalityReduction

    using FlowCytometry, MLJ, UMAP

    pca_ = MLJ.@load PCA pkg=MultivariateStats

"""
    function pca!(fct::FlowCytometryExperiment;
        maxoutdim::Int = 0,
        method::Symbol = :auto,
        variance_ratio::Float64 = 0.99,
        mean = nothing,
        key_added::String = "pca",
        key_used_channels::Union{Nothing,String} = nothing
        )

Perform pca analysis over the channels.

**Arguments**:
 - **fct::FlowCytometryExperiment** FlowCytometryExperiment to perform the pca.

**Keyword Arguments**
 - **maxoutdim = 0** Number of PCS to compute if 0, it will compute NChannels-1
 - **method = :auto** Method to compute pca. Choose between :svd, :cov and :auto.
 - **variance_ratio = 0.99** Variance ratio up to which compute pca.
 - **mean = nothing** If to consider the normalization as performed. nothing will consider the mean not adjusted.
 - **key_added::String = "pca"** Key to be added to object.
 - **key_used_channels::Union{Nothing,String} = nothing** Bool column from .var specifying which channels to use for clustering.

 **Returns**
 Nothing. To the FlowCytometryExperiment object provided, PCs are added to the .obsm[key_added] and information of the algorithm if included in .uns[key_added].
"""
    function pca!(fct::FlowCytometryExperiment;
                    maxoutdim::Int = 0,
                    method::Symbol = :auto,
                    variance_ratio::Float64 = 0.9999,
                    mean = nothing,
                    key_added::String = "pca",
                    key_used_channels::Union{Nothing,String} = nothing
                    )

        if key_used_channels !== nothing
            channels = fct.var[:,key_used_channels]
            if typeof(channels) != Vector{Bool}
                error("key_used_channels should be a column in var of Bool entries specifying which channels to use for clustering.")
            end
            X = fct.X[:,channels]
        else
            X = fct.X
        end
        
        model = pca_(
                    maxoutdim = maxoutdim,
                    method = method,
                    variance_ratio = variance_ratio,
                    mean = mean
                    )
        mach = machine(model,X,scitype_check_level=0)
        fit!(mach)
        fct.obsm[key_added] = mach.fitresult.proj

        fct.uns[key_added] = Dict([
                                "params"=>Dict([i=>j for (i,j) in zip(keys(params(mach)),params(mach))]),
                                "loadings"=>mach.fitresult.prinvars,
                                "variance_explained"=>mach.fitresult.prinvars./(sum(std(X,dims=1).^2)),
                                "variance_explained_norm"=>mach.fitresult.prinvars./sum(mach.fitresult.prinvars),
                                ])
        
        return
    end

"""
    function umap!(fct::FlowCytometryExperiment;
        n_components::Int = 2,
        kwargs...
        )

Function that computes the UMAP plot of the neighbors.

**Arguments**
 - **fct::FlowCytometryExperiment** FlowCytometryExperiment to perform the UMAP.

**Keyword Arguments**
 - **n_components::Int = 2** Number of dimensions for the projection.
 - **kwargs...** Keyword argumetns for the UMAP_ function. Do ?FlowCytometry.UMAP_ for more details.
"""
    function umap!(fct::FlowCytometryExperiment;
                    n_components::Int = 2,
                    key_added::String = "pca",
                    kwargs...
                    )

        if key_obsm !== nothing && key_used_channels !== nothing
            error("key_obsm and key_used_channels cannot be specified at the same time.")
        elseif key_osbm !== nothing
            if n_components !== nothing
                X = fct.obsm[key_obsm][:,1:n_components]
            else
                X = fct.obsm[key_obsm]
            end
        else
            if key_used_channels !== nothing
                channels = fct.var[:,key_used_channels]
                if typeof(channels) != Vector{Bool}
                    error("key_used_channels should be a column in var of Bool entries specifying which channels to use for clustering.")
                end

                X = fct.X[:,channels]
            else
                X = fct.X
            end
        end
        
        proj = UMAP_(X,n_components,kwargs...)

        fct.obsm[key_added] = proj.embedding

        fct.uns[key_added] = proj
        
        return
    end
end