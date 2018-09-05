classdef RCFLfluorescence < handle
    %UNTITLED Fluorescence with RCME
    
    properties(Constant, Access=private)
        DALIKLIS = 5308.837458877; % Conversion between cm^-1 and fs^-1 (cm^-1/fs^-1)
    end
    
    properties
        dataLoader
        paramsRC
        RCfolder
        fi
    end
    
    methods

        function this = RCFLfluorescence(dataLoader, paramsRC, RCfolder, fi)
            this.dataLoader = dataLoader;
            this.paramsRC = paramsRC;
            this.RCfolder = RCfolder;
            this.fi = fi;
        end
        
        function [ freqRC, fluorescenceRC ] = calcFluorescenceSpectrum(this)

            [ density1, density2, t, equil ] = this.loadRCdataForFluorescence( );

            % Finds Green operator/dynamical map
            G = this.getG(density1, density2);
            
            step = this.findStep(t);
            
            % Adds zeros
            [ G ] = this.zeroPadding( G, this.paramsRC.numOfRCZeros );
            
            miu = this.getDipolesFromAngle();

            R = this.getLinearResponse( G, miu, equil );

            [freqRC, fluorescenceRC] = this.fluorescenceFFT(R, step);

            [ freqRC ] = this.shiftFreq( freqRC, this.paramsRC.energyShift );
            
            %[ fluorescenceRC ] = this.multiplyFluorescence( freqRC, fluorescenceRC );

        end
        
        function [ density1, density2, t, equil ] = loadRCdataForFluorescence( this )
            [ density1, t1, equil1 ] = this.dataLoader.loadRCdata( this.RCfolder, ...
                this.paramsRC.alfa, this.paramsRC.dimRC, 1 );
            [ density2, t, equil ] = this.dataLoader.loadRCdata( this.RCfolder, ...
                this.paramsRC.alfa, this.paramsRC.dimRC, 2 );
            
            this.checkT(t1, t);
            this.checkEquil(equil1, equil);
        end
        function checkT(this, t1, t2)
            if ~isequal(t1, t2)
                warning('fluor:diffT', 'RC failuose laikai nesutampa');
            end
        end
        function checkEquil(this, equil1, equil2)
            if ~isequal(equil1, equil2)
                warning('fluor:diffEquil', 'RC failuose pusiausvyros busenos nesutampa');
            end
        end
        
        function step = findStep(this, t)
            step = t(2) - t(1);
        end
        
        function [ G ] = getG( this, density1, density2 )
            
            timePoints = size(density1, 1);
            A = [density1(1,2,1) density1(1,3,1);density2(1,2,1) density2(1,3,1)];
            for i=1:timePoints
                B1 = [density1(i,2,1); density2(i,2,1)];
                B2 = [density1(i,3,1); density2(i,3,1)];
                X1 = linsolve(A, B1);
                X2 = linsolve(A, B2);
                G{3,2}(i,1) = X2(1);
                G{2,3}(i,1) = X1(2);
                G{2,2}(i,1) = X1(1);
                G{3,3}(i,1) = X2(2);
            end
            
            % Or change everything to
            %G{2,2}(:,1) = density1(:,2,1);
            %G{3,2}(:,1) = density1(:,3,1);
            %G{2,3}(:,1) = density2(:,2,1);
            %G{3,3}(:,1) = density2(:,3,1);
        end
        
        function [ G ] = zeroPadding( this, G, numOfZeros )
            for i = 2:3
                for j = 2:3
                    G{i,j} = cat(1, G{i,j}, zeros(numOfZeros, 1));
                end
            end
        end
        
        function miu = getDipolesFromAngle(this)
            miu{2} = [1 0 0];
            miu{3} = [cos(this.fi*pi/180) sin(this.fi*pi/180) 0];
        end

        function [ R ] = getLinearResponse( this, G, miu, equil )
            R = zeros(size(G{3,2}, 1), 1);
            for a=2:3
                for b=2:3
                    for c=2:3
                        R = R + (1./3).*dot(miu{b}, miu{c}).*G{c,a}.*equil(a-1,b-1);
                    end
                end
            end
        end
        
        function [ w, F ] = fluorescenceFFT( this, R, step )

            %Frequency step
            N = size(R, 1);
            stepW = 2.*pi/(N*step);

            % Frequency axist after fftshift
            w = 0:stepW:(N-1)*stepW;
            if mod(N, 2) == 0
                w = w - N/2*stepW;
            else
                w = w - (N-1)*stepW/2;
            end
            w = w*this.DALIKLIS;

            F =  conj(fftshift( fft( conj(R) ))) ;
        end

        function freqRC = shiftFreq( this, freqRC, energyShift )
            freqRC = freqRC + energyShift;
        end
        
        function [ fluorescenceRC ] = multiplyFluorescence( this, freqRC, fluorescenceRC )
            
            fluorescenceRC = fluorescenceRC.*freqRC';

        end

    end
    
end

