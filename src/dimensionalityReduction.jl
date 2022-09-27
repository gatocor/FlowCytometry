module DimensionalityReduction

    using FlowCytometry, MLJ, UMAP, DataFrames

    pca_ = MLJ.@load PCA pkg=MultivariateStats verbosity=0

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
                    variance_ratio::Float64 = 1.,
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
        MLJ.fit!(mach,verbosity=0)
        fct.obsm[key_added] = MLJ.matrix(MLJ.transform(mach))

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
        n_dims::Int = 2,
        key_added::String = "umap",
        key_obsm::Union{Nothing,String} = nothing,
        n_components::Union{Nothing,Int} = nothing,
        key_used_channels::Union{Nothing,String} = nothing,
        kwargs...
        )

Function that computes the UMAP plot of the neighbors.

**Arguments**
 - **fct::FlowCytometryExperiment** FlowCytometryExperiment to perform the UMAP.

**Keyword Arguments**
 - **n_dims::Int = 2** Number of dimensions for the projection.
 - **key_added::String = "KMeans"** Key to be added to object.
 - **key_obsm::Union{Nothing,String} = nothing** Matrix in .obsm to use for the clustering. Only key_obsm or key_used_channels can be specified.
 - **n_components::Union{Nothing,Int} = nothing** Number of components to use from the .obsm.
 - **key_used_channels::Union{Nothing,String} = nothing** Bool column from .var specifying which channels to use for clustering.
 - **kwargs...** Keyword argumetns for the UMAP_ function. Do ?FlowCytometry.UMAP_ for more details on the keyword arguments.
"""
    function umap!(fct::FlowCytometryExperiment;
                    n_dims::Int = 2,
                    key_added::String = "umap",
                    key_obsm::Union{Nothing,String} = nothing,
                    n_components::Union{Nothing,Int} = nothing,
                    key_used_channels::Union{Nothing,String} = nothing,
                    kwargs...
                    )

        if key_obsm !== nothing && key_used_channels !== nothing
            error("key_obsm and key_used_channels cannot be specified at the same time.")
        elseif key_obsm !== nothing
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

        X = permutedims(X)
        proj = UMAP_(X,n_dims,kwargs...)
        X = permutedims(X)

        fct.obsm[key_added] = permutedims(proj.embedding)

        fct.uns[key_added] = proj
        
        return
    end
end