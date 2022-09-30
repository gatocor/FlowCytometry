module Clustering

    using FlowCytometry, MLJ, Statistics, LinearAlgebra, Distances

    kmeans_ = MLJ.@load KMeans pkg=Clustering verbosity=0
    # agglomerative_ = MLJ.@load AgglomerativeClustering pkg=ScikitLearn verbosity=0

"""
    function kmeans!(fct::FlowCytometryExperiment;
        k::Int = 2, 
        metric = Distances.SqEuclidean(0.0),
        init::Symbol = :kmpp,
        key_added::String = "kmeans",
        key_obsm::Union{Nothing,String} = nothing,
        n_components::Union{Nothing,Int} = nothing,
        key_used_channels::Union{Nothing,String} = nothing)


Function that clusters the data using the KMeans algorithm.
To help asses the number of clusters you can use **KMeansTuning** function.

**Arguments**
 - **fct::FlowCytometryExperiment** FlowCytometryExperiment to cluster the data.

**Keyword Arguments**
 - **k::Int** Number of clusters to partition the data.
 - **metric = Distances.SqEuclidean(0.0)** Metric to compute distances. 
 - **init::Symbol = :kmpp** Initialisation algorithm of the KMeans. Choose between :kmpp, :rand or :kmenc. ?Clustering.kmeans_ for more information.
 - **key_added::String = "KMeans"** Key to be added to object.
 - **key_obsm::Union{Nothing,String} = nothing** Matrix in .obsm to use for the clustering. Only key_obsm or key_used_channels can be specified.
 - **n_components::Union{Nothing,Int} = nothing** Number of components to use from the .obsm.
 - **key_used _channels::Union{Nothing,String} = nothing** Bool column from .var specifying which channels to use for clustering.
 
**Returns**
 Nothing. To the FlowCytometryExperiment object provided, clusters identities are added to the .obs[key_added] and information of the algorithm if included in .uns[key_added].
"""
    function kmeans!(fct::FlowCytometryExperiment;
                        k::Int = 2, 
                        metric = Distances.SqEuclidean(0.0),
                        init::Symbol = :kmpp,
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
                    k=k, 
                    metric=metric, 
                    init=init
                    )
        mach = machine(model,X,scitype_check_level=0)
        fit!(mach,verbosity=0)
        fct.obs[:,key_added] = mach.report[1]

        fct.uns[key_added] = Dict([
                                "params"=>Dict([i=>j for (i,j) in zip(keys(params(mach)),params(mach))]),
                                "cluster_centers"=>permutedims(fitted_params(mach)[1])
                                ])
        
        return
    end

"""
    function agglomerative!(fct::FlowCytometryExperiment;
        k::Int = 2, 
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
 - **k::Int = 2** Number of clusters to parition the data.
 - **affinity::String = "euclidean"** Metric to use for partitioning the data.
 - **compute_full_tree::Union{String,Bool} = "auto"** Stop early the construction of the tree at k. This is useful to decrease computation time if the number of clusters is not small compared to the number of samples.
 - **linkage::String = "ward"** Linkage criterion.
 - **distance_threshold::Union{Nothing,Real} = nothing** The linkage distance threshold above which, clusters will not be merged. If not None, k must be None and compute_full_tree must be True.
 - **key_added::String = "KMeans"** Key to be added to object.
 - **key_obsm::Union{Nothing,String} = nothing** Matrix in .obsm to use for the clustering. Only key_obsm or key_used_channels can be specified.
 - **n_components::Union{Nothing,Int} = nothing** Number of components to use from the .obsm.
 - **key_used _channels::Union{Nothing,String} = nothing** Bool column from .var specifying which channels to use for clustering.
 
**Returns**
 Nothing. To the FlowCytometryExperiment object provided, clusters identities are added to the .obs[key_added] and information of the algorithm if included in .uns[key_added].
"""
    function agglomerative!(fct::FlowCytometryExperiment;
                            k = 2, 
                            affinity = "euclidean", 
                            memory = nothing, 
                            connectivity = nothing, 
                            compute_full_tree = "auto", 
                            linkage = "ward", 
                            distance_threshold = nothing,
                            key_added::String = "agglomerative",
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
                    k = k, 
                    affinity = affinity, 
                    memory = memory, 
                    connectivity = connectivity, 
                    compute_full_tree = compute_full_tree, 
                    linkage = linkage, 
                    distance_threshold = distance_threshold
                    )
        mach = machine(model,X,scitype_check_level=0)
        fit!(mach,verbosity=0)
        fct.obs[:,key_added] = fitted_params(mach)[2]

        fct.uns[key_added] = Dict([
                                "params"=>Dict([i=>j for (i,j) in zip(keys(params(mach)),params(mach))]),
                                "k"=>fitted_params(mach)[1],
                                "n_leaves"=>fitted_params(mach)[3],
                                "n_connected_components"=>fitted_params(mach)[4],
                                "children"=>fitted_params(mach)[5],
                                ])
        
        return
    end

    function closerPoints(X::Matrix,centers::Matrix)
        if size(X)[2] != size(centers)[2]
            error("Second dimension")
        end
        d = zeros(size(X)[1],size(centers)[1])

        for i in 1:size(X)[1]
            x = @views X[i,:]
            for j in 1:size(centers)[1]
                c = @views centers[j,:]
                d[i,j] = sum((x.-c).^2)
            end
        end

        return Matrix([i[2] for i in argmin(d,dims=2)])[:,1]
    end

    function gmloglikelihood!(p::Matrix,X::Matrix,centers::Matrix,covariances::Vector,weights::Vector)

        for j in 1:size(centers)[1]
            c = @views centers[j,:]
            sigma = covariances[j]
            sigmaInv = inv(sigma)
            determinant = abs(det(sigma))
            w = weights[j]
            for i in 1:size(X)[1]
                x = @views X[i,:]
                m = (x-c).^2
                p[i,j] = -(transpose(m)*sigmaInv*m/2)[1]-log(determinant)/2+log(w)
            end
        end

        return
    end

"""
    function gaussianMixture!(fct::FlowCytometryExperiment;
        k::Int,
        initialization::Union{String,Matrix} = "kmeans",
        maximumSteps::Int = 10000,
        key_added::String = "gaussianMixture",
        key_obsm::Union{Nothing,String} = nothing,
        n_components::Union{Nothing,Int} = nothing,
        key_used_channels::Union{Nothing,String} = nothing                                
        )

Clustering with multivariate gaussian distributions.

**Arguments**
 - **fct::FlowCytometryExperiment** FlowCytometryExperiment to cluster the data.

**Keyword Arguments**
 - **k::Int** Number of clusters of the model.
 - **initialization::Union{String,Matrix} = "kmeans"** Initialization. Select between "kmeans", "random" or give a matrix of (k,dimensions)
 - **maximumSteps::Int = 10000** Number of maximum of steps before stopping the algorithm.
 - **key_added::String = "KMeans"** Key to be added to object.
 - **key_obsm::Union{Nothing,String} = nothing** Matrix in .obsm to use for the clustering. Only key_obsm or key_used_channels can be specified.
 - **n_components::Union{Nothing,Int} = nothing** Number of components to use from the .obsm.
 - **key_used _channels::Union{Nothing,String} = nothing** Bool column from .var specifying which channels to use for clustering.

 **Returns**
 Nothing. To the FlowCytometryExperiment object provided, clusters identities are added to the .obs[key_added] and information of the algorithm if included in .uns[key_added].
"""
    function gaussianMixture!(fct::FlowCytometryExperiment;
                                k::Int,
                                initialization::Union{String,Matrix} = "kmeans",
                                maximumSteps::Int = 10000,
                                key_added::String = "gaussianMixture",
                                key_obsm::Union{Nothing,String} = nothing,
                                n_components::Union{Nothing,Int} = nothing,
                                key_used_channels::Union{Nothing,String} = nothing                                
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

        #Check shape of initial centers
        if typeof(initialization) == Matrix
            if size(initialisation) != (k,size(X)[2])
                error("If initialization is a matrix of given initial positions, the dimensions must be the same size as (k,variables).")
            end
        end

        nCells = size(X)[1]

        #Initialization centers
        centers = zeros(k,size(X)[2])
        if initialization == "kmeans"
            model = kmeans_(k=k)
            mach = machine(model,X,scitype_check_level=0)
            fit!(mach,verbosity=0)
            centers = permutedims(fitted_params(mach)[1])
        elseif initialization == "random"
            centers = rand(k,size(X)[2])
            centers .= centers .* (maximum(X,dims=1).-minimum(X,dims=1)) .-minimum(X,dims=1)
        elseif typeof(initialization) == Matrix
            centers .= initialization
        else
            error("initialization must be 'kmeans', 'random' or a matrix of specifying the centers of the gaussian.")
        end
        #Initialization of identities
        identities = closerPoints(X,centers)
        #Initialization of covariances
        covariances = []
        for id in 1:k
            push!(covariances,cov(@views X[identities.==id,:]))
        end
        #Initialization of weights
        weights = [sum(identities.==id)/nCells for id in 1:k]

        #Loop
        p = zeros(nCells,k)
        steps = 0
        identitiesNew = fill(0,nCells)
        while !all(identities.==identitiesNew) && steps < maximumSteps
            identities .= identitiesNew
            #Maximization step
            gmloglikelihood!(p,X,centers,covariances,weights)
            identitiesNew .= Matrix([i[2] for i in argmax(p,dims=2)])[:,1]
            #Expectation step
            for i in 1:k
                ids = identitiesNew.==i

                weights[i] = sum(ids)/nCells
                if weights[i] > 0            
                    centers[i,:] .= mean(X[ids,:],dims=1)[1,:]
                    covariances[i] .= cov(X[ids,:])
                end
            end
            steps += 1
        end

        fct.uns[key_added] = Dict([
            "k" => k,
            "maximumSteps" => maximumSteps,
            "stepsBeforeConvergence" => steps,
            "initialization" => initialization,
            "weights" => p
        ])

        fct.obs[!,key_added] = identitiesNew

        return

    end

end