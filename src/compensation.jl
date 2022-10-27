module Compensation

    using MLJ, FlowCytometry, LinearAlgebra, MLJLinearModels, ProgressMeter, Random
    
    linearregressor_ = @load LinearRegressor pkg=MultivariateStats verbosity=0
    huberregressor_ = @load HuberRegressor pkg=MLJLinearModels verbosity=0

"""
    function computeCompensationMatrix!(fcs::FlowCytometryControl; error::AbstractFloat = 10E-5, nMaxiterations::Int=10)

Compute the compensation matrix as proposed by [Roca et al](https://www.nature.com/articles/s41467-021-23126-8).

**Arguments**
 - **fcs::FlowCytometryControl** FlowCytometryControl where to compute the spillover matrix.

**Keyword arguments**
 - **error::AbstractFloat = 10E-5** Maximum error to stop the refinement of the spillover matrix.
 - **nMaxIterations::Int=10** Maximum number of iterations of the refinement matrix. Set to 0 to compute the classical spillover matrix.

**Results**
 Nothing. Adds the compensation matrix to the FlowCytometryControl object.
"""
    function computeCompensationMatrix!(fcs::FlowCytometryControl; error::AbstractFloat = 10E-5, nMaxIterations::Int=10, channelsCompensate = nothing, subsample = 5000, seed = 0)

#        checkControlNames(fcs)
        if channelsCompensate === nothing
            channelsCompensate == fcs.channels
        else
            l = [i for i in channelsCompensate if !(i in fcs.channels)]
            if !isempty(l)
                println(l," not found in channels of FlowCytometryControl object.")
            
                return
            end
        end

        if "compensation" in keys(fcs.uns)
            if "compensated" in keys(fcs.uns["compensation"])
                if fcs.uns["compensation"]["compensated"]
                    println("FlowCytometryControl already compensated.")        
                    
                    return
                end
            end
        end

        #Find channels
        channels = fcs.channels

        #Initialize Spillover matrix
        S = Matrix{Float64}(I,length(channels),length(channels))
    
        fcs.uns["compensation"] = Dict{String,Any}()

        #Create regressor object
        regressor = huberregressor_(fit_intercept=false)
        #Make S(0)
        @showprogress 1 "Computing S(0)..." for (channel1,control) in pairs(fcs.controls)
            Random.seed!(seed)
            sub = shuffle(1:size(control.X)[1])[1:min(size(control.X)[1],subsample)]
            i = findfirst(control.channels .== channel1)
            iS = findfirst(channels .== channel1)
            x = MLJ.table(fcs.controls[channel1].X[sub,i:i])
            for channel2 in channelsCompensate
                j = findfirst(control.channels .== channel2)
                jS = findfirst(channels .== channel2)
                if j !== nothing
                    if i != j
                        y = fcs.controls[channel1].X[sub,j]
                        mach = machine(regressor,x,y,scitype_check_level=0)
                        fit!(mach,verbosity=0)
                        S[iS,jS] = mach.fitresult[1][1] 
                    end
                end
            end
        end
        fcs.uns["compensation"]["spilloverMatrix"] = copy(S)
        #Make correction matrix
        Sinv = inv(S)

        #Make refinements S(t)
        eMax = 0
        E = zeros(length(channels),length(channels))
        count = 0
        @showprogress 1 "Refinement Iterations..." for i in 1:1:nMaxIterations
            E .= 0
            for (channel1,control) in pairs(fcs.controls)
                Random.seed!(seed)
                sub = shuffle(1:size(control.X)[1])[1:min(size(control.X)[1],subsample)]
                i = findfirst(control.channels .== channel1)
                iS = findfirst(channels .== channel1)
                channelsPresent = [i in control.channels for i in channels]
                X = control.X * Sinv[:,channelsPresent][channelsPresent,:]
                x = MLJ.table(X[sub,i:i])
                for channel2 in channelsCompensate
                    j = findfirst(control.channels .== channel2)
                    jS = findfirst(channels .== channel2)
                    if j !== nothing
                        if i != j
                            y = X[sub,j]

                            mach = machine(regressor,x,y,scitype_check_level=0)
                            fit!(mach,verbosity=0)
                            E[iS,jS] = mach.fitresult[1][1] 
                        end
                    end
                end
            end

            #Break is converged
            eMax = maximum(E)
            if eMax < error 
                break
            end
        
            #Add correction
            S .+= E*S
            #Normalize
            for j in 1:length(channels)
                for k in 1:length(channels)
                    if j != k
                        S[j,k] /= S[j,j]
                    end
                end
                S[j,j] = 1
            end
            fcs.uns["compensation"]["spilloverMatrix"] = copy(S)
            Sinv = inv(S)
            count += 1
        end

        fcs.compensationMatrix = CompensationMatrix(fcs.channels,Sinv)
        fcs.uns["compensation"]["error"] = error
        fcs.uns["compensation"]["errorMax"] = eMax
        fcs.uns["compensation"]["nMaxIterations"] = Int(nMaxIterations)
        fcs.uns["compensation"]["convergedIterations"] = Int(count)
        fcs.uns["compensation"]["compensated"] = false

        return

    end

"""
    function compensate!(fcs::FlowCytometryControl)

Compensate all the control experiments with the compensationMatrix of the object

**Arguments**
 - **fcs::FlowCytometryControl** FlowCytometryControl object where to compute the compensation.

**Returns**
Nothing
"""
    function compensate!(fcs::FlowCytometryControl)
   
        if fcs.compensationMatrix !== nothing 
            if !(fcs.uns["compensation"]["compensated"])
                for (_,control) in pairs(fcs.controls)
                    l = [i in control.channels for i in fcs.channels]
                    control.X = control.X*fcs.compensationMatrix.matrix[:,l][l,:]
                end
                
                fcs.uns["compensation"]["compensated"] = true
            else
                println("FlowCytometryControl object already compensated. Ignoring.")
            end
        else
            error("computeCompensation! should have been called before compensating.") 
        end
       
        return
        
    end
    
"""
    function compensate!(fcs::FlowCytometryExperiment,compensation::CompensationMatrix)

Compensate FlowCytometryExperiment experiment with the compensationMatrix from an FlowCytometryControl object.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment object where to compute the compensation.
 - **compensation::CompensationMatrix** CompensationMatrix object

**Returns**
Nothing, compensated the matrix X from fcs.
"""
    function compensate!(fcs::FlowCytometryExperiment,compensation::CompensationMatrix)
        
        if "compensation" in keys(fcs.uns)
            compensated = true
        else
            compensated = false
        end
        
        if control.compensationMatrix !== nothing && !compensated
            l = [i in fcs.channels for i in compensation.channels]
            fcs.X = fcs.X*compensation.matrix[:,l][l,:] 
            
            fcs.uns["compensation"] = Dict("compensated"=>true)
        elseif compensated
            println("FlowCytometryExperiment is already compensated. Ignoring.")
        else
           error("computeCompensation! should have been called on the Control object before compensating.") 
        end
        
        return
        
    end
    
end