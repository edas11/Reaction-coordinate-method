classdef RCMEcoreBundle < handle
    %RCMEINPUT RCME parameters and basic operators
    
    properties
        DALIKLIS = 5308.837458877; % Conversion between cm^-1 and fs^-1 (cm^-1/fs^-1)
        kB = 0.6946; % Boltzmann constant cm^-1/K
        J % Resonance coupling (fs^-1)
        delta_e % Energy difference between sites (fs^-1)
        n % Number of terms interracting with RCs (dimensionless units)
        dim % System dimension (dimensionless units)
        dimRC % RC dimension (dimensionless units)
        beta % Inverse temperature (1/fs^-1)
        gamma % Relaxation rate (fs^-1)
        lambda % Reorganization energy (fs^-1)
        reactDaznis % RC frequency (fs^-1)
        eta % Strength of RC-system interraction (dimensionless units)
        operatoriai % System - RCs (bath) interraction operators
        alfa % Strength of new spectral density (dimensionless units)
        aukstinimo % RC creation operator
        zeminimo % RC annihilation operator
    end
    
    methods
        
        function this = RCMEcoreBundle(parameterStruct)
            this.takeFromParameterStruct(parameterStruct)
        end
        
        function takeFromParameterStruct(this, parameterStruct)
            fieldNames = fieldnames(parameterStruct);
            for fieldNr = 1:length(fieldNames)
                switch fieldNames{fieldNr}
                    case 'dim'
                        this.dim = parameterStruct.dim;
                    case 'J'
                        this.J = parameterStruct.J;
                    case 'T'
                        this.beta = (this.DALIKLIS)/(parameterStruct.T*(this.kB));
                    case 'gamma'
                        this.gamma = parameterStruct.gamma;
                    case 'lambda'
                        this.lambda = parameterStruct.lambda;
                    case 'delta_e'
                        this.delta_e = parameterStruct.delta_e;
                    case 'operatoriai'
                        this.operatoriai = parameterStruct.operatoriai;
                        this.n = length(this.operatoriai);
                    case 'alfa'
                        this.alfa = parameterStruct.alfa;
                    case 'dimRC'
                        this.dimRC = parameterStruct.dimRC;
                    otherwise
                        warning(strcat('Nezinomas kintamasis: ', fieldNames{fieldNr}));
                end
            end
        end
        
        function [H0s] = hamiltonianSystem(this)
            if this.getDim()==1 ||  this.getDim()-3>0.0001
                error('HsError:badDim', 'Sistemos hamiltoniana galima sukurti tik 2 ar 3 dim');
            else
                H0s = zeros(this.getDim(), this.getDim());
                H0s(end-1:end,end-1:end) = [this.delta_e this.J; this.J  0];
            end
        end
        
        function spc = spectral(this, omega)
            spc = this.alfa*omega;
        end
        
        function spc = spectralDivOmega(this, omega)
            spc = this.alfa;
        end
        
        function HRC = hamiltonianRC(this)
            if isempty(this.aukstinimo)
                [this.aukstinimo, this.zeminimo] = RCMEutil.aukstinimoZeminimo(this.dimRC);
            end
            HRC = this.getReactDaznis()*this.aukstinimo*this.zeminimo;
        end
        
        function coord = coordRC(this)
            if isempty(this.aukstinimo)
                [this.aukstinimo, this.zeminimo] = RCMEutil.aukstinimoZeminimo(this.dimRC);
            end
            coord = (this.aukstinimo + this.zeminimo)/sqrt(2);
        end
        
        function I = identitySystem(this)
            I = eye(this.dim);
        end
        
        function I = identityRC(this)
            I = eye(this.dimRC);
        end
        
        function n = getN(this)
            n=this.n;
        end
        
        function dim = getDim(this)
            dim=this.dim;
        end
        
        function dimRC = getDimRC(this)
            dimRC=this.dimRC;
        end
        
        function dimRCs = getDimRCs(this)
            dimRCs=this.dimRC^this.n;
        end
        
        function dimAll = getDimAll(this)
            dimAll = this.dim*this.dimRC^this.n;
        end
        
        function beta = getBeta(this)
            beta = this.beta;
        end
        
        function operatoriai = getOperatoriai(this)
            operatoriai = this.operatoriai;
        end
        
        function reactDaznis = getReactDaznis(this)
            if isempty(this.reactDaznis)
                this.reactDaznis = (this.alfa)*(this.gamma);
            end
            reactDaznis = this.reactDaznis;
        end
        
        function eta = getEta(this)
            if isempty(this.eta)
                this.eta = sqrt( 2*(this.lambda)*(this.getReactDaznis()) );
            end
            eta = this.eta;
        end
        
    end
    
end

