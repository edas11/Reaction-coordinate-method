classdef RCABdataLoader < handle
    %RCFLDATALOADER Load data
    
    properties (Access=private)
        params
        fi
    end
    
    properties (Constant, Access=private)
        RCnames = {'(T)' '(g-1)' '(d)' '(l)' '(J)' '(RCalfa)' '(RCdim)'};
        HEOMnames = {'(T)' '(g-1)' '(d)' '(l)' '(J)' '(fi)'};
    end
    
    methods
        
        function this = RCABdataLoader(params, fi)
            this.params = params;
            this.fi = fi;
        end
        
        function [ density, t ] = loadRCdata( this, RCfolder, alfa, dimRC, IC )
            [ RCfileName ] = createRCfileName( this, alfa, dimRC, IC);
            fullPath = strcat(RCfolder, RCfileName);
            RCfileContent = open(fullPath);
            density = RCfileContent.densityReturned;
            t = RCfileContent.t;
        end
        
        function [ RCfileName ] = createRCfileName( this, alfa, dimRC, IC)
            RCfileName = strcat('RC-', this.createRCparamString(alfa, dimRC), ...
                '(IC)', num2str(IC), '.mat');
        end
        
        function [ RCparamString ] = createRCparamString( this, alfa, dimRC )
            parameters = [this.params.T 1/this.params.gamma this.params.delta_e...
                this.params.lambda this.params.J alfa dimRC];
            RCparamString = this.createParamString(this.RCnames, parameters);
        end

        function [ freqHEOM, fluorescenceHEOM ] = loadHEOMdata( this, HEOMfolder )
            fullPath = strcat(HEOMfolder, this.createHEOMfileName());
            HEOMdata = importdata( fullPath,'\t' );
            freqHEOM = HEOMdata(:, 1) ;
            fluorescenceHEOM = HEOMdata(:, 5) + 1i*HEOMdata(:, 6) ;
        end
        
        function [ HEOMfileName ] = createHEOMfileName( this )
            HEOMparamString = this.createHEOMparamString();
            HEOMfileName = strcat('HEOM-', HEOMparamString, '.txt');
        end
        
        function [ HEOMparamString ] = createHEOMparamString( this )
            parameters = [this.params.T 1/this.params.gamma this.params.delta_e...
                this.params.lambda this.params.J this.fi];
            HEOMparamString = this.createParamString(this.HEOMnames, parameters);
        end

    end
    
    methods (Static)
        
        function [ paramString ] = createParamString( names, parameters )
            paramString = '';
            for i=1:length(names)
                paramString = strcat(paramString, names{i}, num2str(parameters(i)));
            end
        end
        
    end
    
end

