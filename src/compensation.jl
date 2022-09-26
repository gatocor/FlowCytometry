module Compensation

    using MLJ, FlowCytometry, LinearAlgebra
    
    linearregressor_ = @load LinearRegressor pkg=MultivariateStats

"""
    function computeSpilloverMatrix!(fcs::FlowCytometryControl; error::AbstractFloat = 10E-5, nMaxiterations::Int=10)

Compute the spillover matrix as proposed by [Roca et al](https://www.nature.com/articles/s41467-021-23126-8).

**Arguments**
 - **fcs::FlowCytometryControl** FlowCytometryControl where to compute the spillover matrix.

**Keyword arguments**
 - **error::AbstractFloat = 10E-5** Maximum error to stop the refinement of the spillover matrix.
 - **nMaxIterations::Int=10** Maximum number of iterations of the refinement matrix. Set to 0 to compute the classical spillover matrix.

**Results**
 Nothing. Adds the spillover matrix to the FlowCytometryControl object.
"""
    function computeSpilloverMatrix!(fcs::FlowCytometryControl; error::AbstractFloat = 10E-5, nMaxIterations::Int=10)

        checkControlNames(fcs)

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
            X = fcs.controls[channel1].X[:,i:i]
            for channel2 in controls
                j = findfirst(channels .== channel2)
                if i != j
                    Y = fcs.controls[channel1].X[:,j]

                    mach = machine(regressor,X,Y,scitype_check_level=0)
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
        for i in 1:nMaxIterations
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
            #Add correction
            S .+= E
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

            #Break is converged
            eMax = maximum(E)
            if eMax < error 
                break
            end
        end

        fcs.spillover = Sinv
        fcs.uns["spillover"] = Dict(["error"=>error,"errorMax"=>eMax,"nMaxIterations"=>nMaxIterations])

        return

    end

end