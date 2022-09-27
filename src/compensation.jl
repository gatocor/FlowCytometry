module Compensation

    using MLJ, FlowCytometry, LinearAlgebra
    
    linearregressor_ = @load LinearRegressor pkg=MultivariateStats verbosity=0

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
    function computeCompensationMatrix!(fcs::FlowCytometryControl; error::AbstractFloat = 10E-5, nMaxIterations::Int=10)

        checkControlNames(fcs)

        if "compensation" in keys(fcs.uns)
            if fcs.uns["compensation"]["compensated"]
                error("FlowCytometryControl already compensated.")
            end
        end

        #Find channels
        channels = fcs.controls[iterate(fcs.controls)[1][1]].var.channel
        #controls
        controls = [i for (i,_) in pairs(fcs.controls)]

        #Initialize Sipillover matrix
        S = Matrix{Float64}(I,length(channels),length(channels))

        #Create regressor object
        regressor = linearregressor_(bias=false)
        #Make S(0)
        for channel1 in controls
            i = findfirst(channels .== channel1)
            x = fcs.controls[channel1].X[:,i:i]
            for channel2 in controls
                j = findfirst(channels .== channel2)
                if i != j
                    y = fcs.controls[channel1].X[:,j]

                    mach = machine(regressor,x,y,scitype_check_level=0)
                    fit!(mach,verbosity=0)
                    S[i,j] = mach.fitresult.sol_matrix[1] 
                end
            end
        end
        #Make correction matrix
        Sinv = inv(S)

        #Make refinements S(t)
        eMax = 0
        E = zeros(length(channels),length(channels))
        count = 0
        for i in 1:1:nMaxIterations
            E .= 0
            for channel1 in controls
                i = findfirst(channels .== channel1)
                X = fcs.controls[channel1].X *Sinv
                x = X[:,i:i]
                for channel2 in controls
                    j = findfirst(channels .== channel2)
                    if i != j
                        y = X[:,j]

                        mach = machine(regressor,x,y,scitype_check_level=0)
                        fit!(mach,verbosity=0)
                        E[i,j] = mach.fitresult.sol_matrix[1] 
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
            Sinv = inv(S)
            count += 1
        end

        fcs.compensationMatrix = Sinv
        fcs.uns["compensation"] = Dict{String,Any}("error"=>error,"errorMax"=>eMax,"nMaxIterations"=>Int(nMaxIterations),"convergedIterations"=>Int(count),"compensated"=>false)

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
                for (control,_) in pairs(fcs.controls)
                    fcs.controls[control].X = fcs.controls[control].X*fcs.compensationMatrix 
                end
                
                fcs.uns["compensation"]["compensated"] = true
            end
        elseif fcs.uns["compensation"]["compensated"]
            error("FlowCytometryControl object already compensated.")
        else
            error("computeCompensation! should have been called before compensating.") 
        end
       
        return
        
    end
    
"""
    function compensate!(fcs::FlowCytometryExperiment;control::FlowCytometryControl)

Compensate FlowCytometryExperiment experiment with the compensationMatrix from an FlowCytometryControl object.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment object where to compute the compensation.

**Keyword Arguments**
 - **control::FlowCytometryControl** FlowCytometryControl to use the compensationMatrix

**Returns**
Nothing
"""
    function compensate!(fcs::FlowCytometryExperiment;control::FlowCytometryControl)
        
        if "compensation" in keys(fcs.uns)
            compensated = true
        else
            compensated = false
        end
        
        if control.compensationMatrix !== nothing && !compensated
            fcs.X = fcs.X*control.compensationMatrix 
            
            fcs.uns["compensation"] = Dict("compensated"=>true)
        elseif compensated
            error("FlowCytometryExperiment is already compensated.")
        else
           error("computeCompensation! should have been called on the Control object before compensating.") 
        end
        
        return
        
    end
    
"""
    function compensate!(fcs::FlowCytometryExperiment)

Compensate FlowCytometryExperiment experiment with an assigned compensationMatrix from a FlowCytometryControl object. See assignCompensation!.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment object where to compute the compensation.

**Returns**
Nothing
"""
    function compensate!(fcs::FlowCytometryExperiment)
       
        if !("compensation" in keys(fcs.uns))
            error("No compensation matrix in the object. Assign a compensation matrix before calling the function.")
        elseif !(fcs.uns["compensation"]["compensated"])
            fcs.X = fcs.X*fcs.uns["compensation"]["compensationMatrix"] 
            
            fcs.uns["compensation"]["compensated"] = true
        elseif  fcs.uns["compensation"]["compensated"]
            error("FlowCytometryExperiment is already compensated.")
        else
           error("computeCompensation! should have been called on the Control object before compensating.") 
        end
        
        return
    
    end
    
"""
    function assignCompensation!(fcs::FlowCytometryExperiment;control::FlowCytometryControl)

Assign a compensation mattrix to a FlowCytometryExperiment experiment from the compensationMatrix of a FlowCytometryControl object.

**Arguments**
 - **fcs::FlowCytometryExperiment** FlowCytometryExperiment object where to compute the compensation.

**Keyword Arguments**
 - **control::FlowCytometryControl** FlowCytometryControl to use the compensationMatrix

**Returns**
Nothing
"""
    function assignCompensation!(fcs::FlowCytometryExperiment;control::FlowCytometryControl)
        
        if control.compensationMatrix !== nothing 
            
            if "compensation" in keys(fcs.uns)
                if fcs.uns["compensation"]["compensated"]
                    error("FlowCytometryExperiment is already compensated.")
                end
            end
            
            fcs.uns["compensation"] = deepcopy(control.uns["compensation"])
            fcs.uns["compensation"]["compensationMatrix"] = control.compensationMatrix
            fcs.uns["compensation"]["compensated"] = false
            
        elseif control.compensationMatrix === nothing

            error("FlowCytometryControl control object does not have any compesation matrix computed. Compute it first before assignement.")

        else

            error("FlowCytometryExperiment already compensated, cannot assign a new compensation matrix.")
        
        end
        
        return
        
    end

end