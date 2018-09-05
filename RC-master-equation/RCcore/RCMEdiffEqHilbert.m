classdef RCMEdiffEqHilbert
    %RCMEDIFFEQHILBERT Master equation, Hilbert space.
    
    properties
        Q
        H0
        Xim
        Xid
    end
    
    methods
        
        function this = RCMEdiffEqHilbert(operatorBuilder)
            [ this.Q, this.H0, this.Xim, this.Xid ] = operatorBuilder.getRCMEoperators();
        end
        
        function drodtVector = nextDensDerivativeVector(this, densVector)
            densMatrix = this.densMatrixFromVector(densVector);
            drodtMatrix = this.DensDerivative(densMatrix);
            drodtVector = this.densVectorFromMatrix(drodtMatrix);
        end
        
        function drodtMatrix = DensDerivative(this, densMatrix)
            drodtMatrix = -1i*RCMEmath.commutator(this.H0, densMatrix);
            for i=1:length(this.Q)
                drodtMatrix = drodtMatrix - RCMEmath.commutator(this.Q{i}, RCMEmath.commutator(this.Xim{i}, densMatrix)) + ...
                    RCMEmath.commutator(this.Q{i}, RCMEmath.antiCommutator(this.Xid{i}, densMatrix));
            end
        end
        
        function densMatrix = densMatrixFromVector(this, densVector)
            dimPilnas = sqrt(length(densVector));
            densMatrix = reshape(densVector, [dimPilnas dimPilnas]);
        end
        
        function densVector = densVectorFromMatrix(this, densMatrix)
            densVector = densMatrix(:);
        end
        
    end
    
end

