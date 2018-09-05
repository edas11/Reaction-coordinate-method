classdef RCMEdensitySystem < RCMEdensity
    %RCMEDENSITYOPERATOR Density operator results (with time) for system
    %without RCs
    
    properties
        dimSystem
        densitySystem
    end
    
    methods
        
        function this = RCMEdensitySystem(RCMEbundle, t, density)
            this.dimSystem = RCMEbundle.getDim();
            this.densitySystem = zeros(0,this.dimSystem,this.dimSystem);
            
            if exist('t') && exist('density')
                this.addDensity(t, density);
            end
        end
        
        function densitySystem = getDensitySystem(this)
            densitySystem = this.densitySystem;
        end
        
        function addDensityMatrix(this, densityToAdd)
            if this.densityIsPossibleProductOfSystem(size(densityToAdd,2))
                densitySystemToAdd = this.traceLargerToSystem(densityToAdd);
            else
                densitySystemToAdd = densityToAdd;
            end
            this.densitySystem = cat(1, this.densitySystem, densitySystemToAdd);
        end
        
        function densitySystem = traceLargerToSystem(this, densityLarger)
            dimLarger = size(densityLarger, 2);
            densitySystem = RCMEmath.traceOutEnd(this.dimSystem, dimLarger/this.dimSystem, densityLarger);
        end
        
        function isBadDims = badMatrixDims(this, possibleDim1, possibleDim2)
            if possibleDim1~=possibleDim2
                isBadDims = true;
                return
            end
            if possibleDim1==this.dimSystem
                isBadDims = false;
                return
            end
            if this.densityIsPossibleProductOfSystem(possibleDim1)
                isBadDims = false;
                return
            end
            isBadDims = true;
        end
        
        function isLargerDensity = densityIsPossibleProductOfSystem(this, densityDim)
            isLargerDensity = round(densityDim/this.dimSystem)==(densityDim/this.dimSystem);
        end
        
    end
    
end

