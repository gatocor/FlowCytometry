module Clustering

    using FlowCytometry, MLJ

    kmeans_ = MLJ.@load KMeans pkg=ScikitLearn
    agglomerative_ = MLJ.@load AgglomerativeClustering pkg=ScikitLearn

"""
    function kmeansTuning(fct::FlowCytometryExperiment;
                            n_clusters::Vector{Int};,
                            n_init::Int = 10, 
                            max_iter::Int = 300, 
                            tol::Real = 0.0001, 
                            verbose::Int = 0, 
                            random_state::Union{Nothing,Int} = nothing, 
                            copy_x::Bool = true, 
                            algorithm::String = "auto", 
                            init::String = "k-means++",
                            key_obsm::Union{Nothing,String} = nothing,
                            n_components::Union{Nothing,Int} = nothing,
                            key_used_channels::Union{Nothing,String} = nothing)

Function that evaluated the KMeans algorithm for different number of clusters and returns the inertia value for each of them. 
A plot of clusters vs inertia can help to choose the optimal number of clusters using the elbow heuristic.

**Arguments**
 - **fct::FlowCytometryExperiment** FlowCytometryExperiment to cluster the data.

**Keyword Arguments**
 - **n_clusters::Vector{Int}** Vector with the n of cluster partitions to be tested.
 - **n_init::Int = 10** Number of initialisations of the algorithm.
 - **max_iter::Int = 300** Maximum number of steps.
 - **tol::Real = 0.0001** Tolerance step before stopping the algorithm.
 - **verbose::Int = 0** Verbosity of the algorithm. 
 - **random_state::Union{Nothing,Int} = nothing,** Random seed to start the algorithm.
 - **copy_x::Bool = true** If to preserve a copy of X. 
 - **algorithm::String = "auto"** Algorithm for updating the KMeans. 
 - **init::String = "k-means++"** Initialisation algorithm of the KMeans.
 - **key_obsm::Union{Nothing,String} = nothing** Matrix in .obsm to use for the clustering. Only key_obsm or key_used_channels can be specified.
 - **n_components::Union{Nothing,Int} = nothing** Number of components to use from the .obsm.
 - **key_used _channels::Union{Nothing,String} = nothing** Bool column from .var specifying which channels to use for clustering.
 
**Returns**
 - Vector with the inertia coeficients for the the different number of clusters checked.
"""
    function kmeansTuning(fct::FlowCytometryExperiment;
                            n_clusters::Vector{Int},
                            n_init::Int = 10, 
                            max_iter::Int = 300, 
                            tol::Real = 0.0001, 
                            verbose::Int = 0, 
                            random_state::Union{Nothing,Int} = nothing, 
                            copy_x::Bool = true, 
                            algorithm::String = "auto", 
                            init::String = "k-means++",
                            key_obsm::Union{Nothing,String} = nothing,
                            n_components::Union{Nothing,Int} = nothing,
                            key_used_channels::Union{Nothing,String} = nothing)

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

        inertia = zeros(length(n_clusters))
        for (i,clusters) in enumerate(n_clusters)
            model = kmeans_(
                        n_clusters=clusters, 
                        n_init=n_init, 
                        max_iter=max_iter, 
                        tol=tol, 
                        verbose=verbose, 
                        random_state=random_state, 
                        copy_x=copy_x, 
                        algorithm=algorithm, 
                        init=init
                        )
            mach = machine(model,X,scitype_check_level=0)
            fit!(mach)
            inertia[i] = fitted_params(mach)[3]
        end
        
        return inertia
    end

"""
    function kmeans!(fct::FlowCytometryExperiment;
                      n_clusters::Int = 2, 
                      n_init::Int = 10, 
                      max_iter::Int = 300, 
                      tol::Real = 0.0001, 
                      verbose::Int = 0, 
                      random_state::Union{Nothing,Int} = nothing, 
                      copy_x::Bool = true, 
                      algorithm::String = "auto", 
                      init::String = "k-means++",
                      key_added::String = "KMeans",
                      key_obsm::Union{Nothing,String} = nothing,
                      n_components::Union{Nothing,Int} = nothing,
                      key_used_channels::Union{Nothing,String} = nothing)


Function that clusters the data using the KMeans algorithm.
To help asses the number of clusters you can use **KMeansTuning** function.

**Arguments**
 - **fct::FlowCytometryExperiment** FlowCytometryExperiment to cluster the data.

**Keyword Arguments**
 - **n_clusters::Int** Number of clusters to partition the data.
 - **n_init::Int = 10** Number of initialisations of the algorithm.
 - **max_iter::Int = 300** Maximum number of steps.
 - **tol::Real = 0.0001** Tolerance step before stopping the algorithm.
 - **verbose::Int = 0** Verbosity of the algorithm. 
 - **random_state::Union{Nothing,Int} = nothing,** Random seed to start the algorithm.
 - **copy_x::Bool = true** If to preserve a copy of X. 
 - **algorithm::String = "auto"** Algorithm for updating the KMeans. 
 - **init::String = "k-means++"** Initialisation algorithm of the KMeans.
 - **key_added::String = "KMeans"** Key to be added to object.
 - **key_obsm::Union{Nothing,String} = nothing** Matrix in .obsm to use for the clustering. Only key_obsm or key_used_channels can be specified.
 - **n_components::Union{Nothing,Int} = nothing** Number of components to use from the .obsm.
 - **key_used _channels::Union{Nothing,String} = nothing** Bool column from .var specifying which channels to use for clustering.
 
**Returns**
 Nothing. To the FlowCytometryExperiment object provided, clusters identities are added to the .obs[key_added] and information of the algorithm if included in .uns[key_added].
"""
    function kmeans!(fct::FlowCytometryExperiment;
                        n_clusters::Int = 2, 
                        n_init::Int = 10, 
                        max_iter::Int = 300, 
                        tol::Real = 0.0001, 
                        verbose::Int = 0, 
                        random_state::Union{Nothing,Int} = nothing, 
                        copy_x::Bool = true, 
                        algorithm::String = "auto", 
                        init::String = "k-means++",
                        key_added::String = "kmeans",
                        key_obsm::Union{Nothing,String} = nothing,
                        n_components::Union{Nothing,Int} = nothing,
                        key_used_channels::Union{Nothing,String} = nothing)

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
        
        model = kmeans_(
                    n_clusters=n_clusters, 
                    n_init=n_init, 
                    max_iter=max_iter, 
                    tol=tol, 
                    verbose=verbose, 
                    random_state=random_state, 
                    copy_x=copy_x, 
                    algorithm=algorithm, 
                    init=init
                    )
        mach = machine(model,X,scitype_check_level=0)
        fit!(mach)
        fct.obs[:,key_added] = predict(mach)

        fct.uns[key_added] = Dict([
                                "params"=>Dict([i=>j for (i,j) in zip(keys(params(mach)),params(mach))]),
                                "cluster_centers"=>fitted_params(mach)[1],
                                "inertia"=>fitted_params(mach)[3],
                                ])
        
        return
    end

"""
    function agglomerative!(fct::FlowCytometryExperiment;
        n_clusters::Int = 2, 
        affinity::String = "euclidean", 
        compute_full_tree::Union{String,Bool} = "auto", 
        linkage::String = "ward", 
        distance_threshold::Union{Nothing,Real} = nothing,
        key_added::String = "KMeans",
        key_obsm::Union{Nothing,String} = nothing,
        n_components::Union{Nothing,Int} = nothing,
        key_used_channels::Union{Nothing,String} = nothing)


Function that clusters the data using the AgglomerativeClustering algorithm.

**Arguments**
 - **fct::FlowCytometryExperiment** FlowCytometryExperiment to cluster the data.

**Keyword Arguments**
 - **n_clusters::Int = 2** Number of clusters to parition the data.
 - **affinity::String = "euclidean"** Metric to use for partitioning the data.
 - **compute_full_tree::Union{String,Bool} = "auto"** Stop early the construction of the tree at n_clusters. This is useful to decrease computation time if the number of clusters is not small compared to the number of samples.
 - **linkage::String = "ward"** Linkage criterion.
 - **distance_threshold::Union{Nothing,Real} = nothing** The linkage distance threshold above which, clusters will not be merged. If not None, n_clusters must be None and compute_full_tree must be True.
 - **key_added::String = "KMeans"** Key to be added to object.
 - **key_obsm::Union{Nothing,String} = nothing** Matrix in .obsm to use for the clustering. Only key_obsm or key_used_channels can be specified.
 - **n_components::Union{Nothing,Int} = nothing** Number of components to use from the .obsm.
 - **key_used _channels::Union{Nothing,String} = nothing** Bool column from .var specifying which channels to use for clustering.
 
**Returns**
 Nothing. To the FlowCytometryExperiment object provided, clusters identities are added to the .obs[key_added] and information of the algorithm if included in .uns[key_added].
"""
    function agglomerative!(fct::FlowCytometryExperiment;
                            n_clusters = 2, 
                            affinity = "euclidean", 
                            memory = nothing, 
                            connectivity = nothing, 
                            compute_full_tree = "auto", 
                            linkage = "ward", 
                            distance_threshold = nothing,
                            key_added::String = "KMeans",
                            key_obsm::Union{Nothing,String} = nothing,
                            n_components::Union{Nothing,Int} = nothing,
                            key_used_channels::Union{Nothing,String} = nothing)

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
        
        model = agglomerative_(
                    n_clusters = n_clusters, 
                    affinity = affinity, 
                    memory = memory, 
                    connectivity = connectivity, 
                    compute_full_tree = compute_full_tree, 
                    linkage = linkage, 
                    distance_threshold = distance_threshold
                    )
        mach = machine(model,X,scitype_check_level=0)
        fit!(mach)
        fct.obs[:,key_added] = fitted_params(mach)[2]

        fct.uns[key_added] = Dict([
                                "params"=>Dict([i=>j for (i,j) in zip(keys(params(mach)),params(mach))]),
                                "n_clusters"=>fitted_params(mach)[1],
                                "n_leaves"=>fitted_params(mach)[3],
                                "n_connected_components"=>fitted_params(mach)[4],
                                "children"=>fitted_params(mach)[5],
                                ])
        
        return
    end

end