classdef RCMEdensityFull < RCMEdensity
    %RCMEDENSITYOPERATOR Density operator results (with time) for system + all RCs.
    
    properties
        dimSystem
        dimRC
        dimRCs
        dimFull
        n
        densityFull
    end
    
    methods
        
        function this = RCMEdensityFull(RCMEbundle, t, density)
            this.dimSystem = RCMEbundle.getDim();
            this.dimRC = RCMEbundle.getDimRC();
            this.dimRCs = RCMEbundle.getDimRCs();
            this.dimFull = RCMEbundle.getDimAll();
            this.n = RCMEbundle.getN();
            this.densityFull = zeros(0,this.dimFull,this.dimFull);
            
            if exist('t') && exist('density')
                this.addDensity(t, density);
            end
        end
        
        function densityFull = getDensityFull(this)
            densityFull = this.densityFull;
        end
        
        function densitySystem = getDensitySystem(this)
            densitySystem = RCMEmath.traceOutEnd(this.dimSystem, this.dimRCs, this.densityFull );
        end
        
        function densityRC = getDensityRC(this, id)
            if id>this.n || id<1
                errStr = strcat({'Maziausias RC id gali buti 1, o didžiausias '}, {num2str(this.n)})
                error('densError:wrongRCid',  char(errStr));
            end
            dimToTraceBegining = this.dimSystem*this.dimRC^(id-1);
            dimReduced = this.dimFull/dimToTraceBegining;
            densityReduced = RCMEmath.traceOutBegining(dimReduced, dimToTraceBegining, this.densityFull);
            dimToTraceEnd = this.dimRC^(this.n - id);
            densityRC = RCMEmath.traceOutEnd(this.dimRC, dimToTraceEnd, densityReduced);
        end
        
        function addDensityMatrix(this, densityToAdd)
            this.densityFull = cat(1, this.densityFull, densityToAdd);
        end
        
        function isBadDims = badMatrixDims(this, possibleDim1, possibleDim2)
            isBadDims = possibleDim1~=this.dimFull || possibleDim2~=this.dimFull;
        end
        
    end
    
end

