classdef RCFLdataLoader < handle
    %RCFLDATALOADER Load data from files
    
    properties (Access=private)
        params
        fi
    end
    
    properties (Constant, Access=private)
        RCnames = {'(T)' '(g-1)' '(d)' '(l)' '(J)' '(RCalfa)' '(RCdim)'};
        HEOMnames = {'(T)' '(g-1)' '(d)' '(l)' '(J)' '(fi)'};
    end
    
    methods
        
        function this = RCFLdataLoader(params, fi)
            this.params = params;
            this.fi = fi;
        end
        
        function [ density, t, equil ] = loadRCdata( this, RCfolder, alfa, dimRC, sample )
            [ RCfileName ] = createRCfileName( this, alfa, dimRC, sample);
            fullPath = strcat(RCfolder, RCfileName);
            RCfileContent = open(fullPath);
            density = RCfileContent.densityReturned;
            t = RCfileContent.t;
            equil = RCfileContent.equil;
        end
        
        function [ RCfileName ] = createRCfileName( this, alfa, dimRC, sample)
            RCfileName = strcat('RC-', this.createRCparamString(alfa, dimRC), ...
                '(S)', num2str(sample), '.mat');
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
            fluorescenceHEOM = HEOMdata(:, 3) + 1i*HEOMdata(:, 4) ;
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

