classdef RCFLmain < handle
    %RCFLMAIN Main fluorescence class
    
    properties (Access = private)
        params
        fi
        dataLoader
    end
    
    methods
        
        function this = RCFLmain(params, fi)
            this.params = params;
            this.fi = fi;
            
            this.initDataLoader();
        end
        
        % Gives HEOM and RC fluorescence
        function [ freqHEOM, fluorescenceHEOM, freqRC, fluorescenceRC ] =...
                getBothFluorescence(this, paramsRC, RCfolder, HEOMfolder)
            
            [freqHEOM, fluorescenceHEOM] = this.getHEOMfluorescence( HEOMfolder );
            [freqRC, fluorescenceRC ] = this.getRCfluorescence( paramsRC, RCfolder );
        end
        
        function [freqHEOM, fluorescenceHEOM] = getHEOMfluorescence( this, HEOMfolder )
            [freqHEOM, fluorescenceHEOM] = this.dataLoader.loadHEOMdata(HEOMfolder);
        end
        
        function [freqRC, fluorescenceRC ] = getRCfluorescence( this, paramsRC, RCfolder )
            RCFLinstance = RCFLfluorescence(this.dataLoader, paramsRC...
                , RCfolder, this.fi);
            [freqRC, fluorescenceRC ] = RCFLinstance.calcFluorescenceSpectrum();
        end
        
        function initDataLoader(this)
            this.dataLoader = RCFLdataLoader(this.params, this.fi);
        end
    end
    
end

