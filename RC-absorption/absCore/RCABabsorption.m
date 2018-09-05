classdef RCABabsorption < handle
    %UNTITLED Absorption with RCM
    
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

        function this = RCABabsorption(dataLoader, paramsRC, RCfolder, fi)
            this.dataLoader = dataLoader;
            this.paramsRC = paramsRC;
            this.RCfolder = RCfolder;
            this.fi = fi;
        end
        
        function [ freqRC, absorptionRC ] = calcAbsorptionSpectrum(this)

            [ density21, density31, t ] = this.loadRCdataForAbsorption( );

            % Creates Green operator/dynamical map
            G = this.getG(density21, density31);
            
            step = this.findStep(t);
            
            [ G ] = this.zeroPadding( G, this.paramsRC.numOfRCZeros );
            
            miu = this.getDipolesFromAngle();

            R = this.getLinearResponse( G, miu );

            [freqRC, absorptionRC] = this.absorptionFFT(R, step);

            [ freqRC ] = this.shiftFreq( freqRC, this.paramsRC.energyShift );
            
            %[ absorptionRC ] = this.multiplyAbsorption( freqRC, absorptionRC );

        end
        
        function [ density21, density31, t ] = loadRCdataForAbsorption( this )
            [ density21, t1 ] = this.dataLoader.loadRCdata( this.RCfolder, ...
                this.paramsRC.alfa, this.paramsRC.dimRC, 21 );
            [ density31, t ] = this.dataLoader.loadRCdata( this.RCfolder, ...
                this.paramsRC.alfa, this.paramsRC.dimRC, 31 );
            
            this.checkT(t1, t);
        end
        function checkT(this, t1, t2)
            if ~isequal(t1, t2)
                warning('fluor:diffT', 'RC failuose laikai nesutampa');
            end
        end
        
        function step = findStep(this, t)
            step = t(2) - t(1);
        end
        
        % Green operator/ dynamical map
        function [ G ] = getG( this, density21, density31 )
            
            G{3,2} = squeeze(density21(:,3,1));
            G{2,2} = squeeze(density21(:,2,1));
            G{2,3} = squeeze(density31(:,2,1));
            G{3,3} = squeeze(density31(:,3,1));
        end
        
        % Green operator/ dynamical map 2 (not used)
        function [ G ] = getGgeneral( this, density1, density2 )
            
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

        function [ R ] = getLinearResponse( this, G, miu )
            R = zeros(size(G{3,2}));
            for a=2:3
                for b=2:3
                    R = R + (1./3).*dot(miu{a}, miu{b}).*G{a,b};
                end
            end
        end
        
        function [ w, F ] = absorptionFFT( this, R, step )

            % Frequency step
            N = size(R, 1);
            stepW = 2.*pi/(N*step);

            % Frequency axis after fftshift
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
        
        function [ absorptionRC ] = multiplyAbsorption( this, freqRC, absorptionRC )
            
            absorptionRC = absorptionRC.*freqRC';

        end

    end
    
end

