classdef RCABmain < handle
    %RCFLMAIN Main absorption class
    
    properties (Access = private)
        params
        fi
        dataLoader
    end
    
    methods
        
        function this = RCABmain(params, fi)
            this.params = params;
            this.fi = fi;
            
            this.initDataLoader();
        end
        
        % Returns HEOM and RC absorption
        function [ freqHEOM, absorptionHEOM, freqRC, absorptionRC ] =...
                getBothAbsorption(this, paramsRC, RCfolder, HEOMfolder)
            
            [freqHEOM, absorptionHEOM] = this.getHEOMabsorption( HEOMfolder );
            [freqRC, absorptionRC ] = this.getRCabsorption( paramsRC, RCfolder );
        end
        
        function [freqHEOM, absorptionHEOM] = getHEOMabsorption( this, HEOMfolder )
            [freqHEOM, absorptionHEOM] = this.dataLoader.loadHEOMdata(HEOMfolder);
        end
        
        function [freqRC, absorptionRC ] = getRCabsorption( this, paramsRC, RCfolder )
            RCABinstance = RCABabsorption(this.dataLoader, paramsRC...
                , RCfolder, this.fi);
            [freqRC, absorptionRC ] = RCABinstance.calcAbsorptionSpectrum();
        end
        
        function initDataLoader(this)
            this.dataLoader = RCABdataLoader(this.params, this.fi);
        end
    end
    
end

