classdef RCMEdiffEqLiouville
    %RCMEDIFFEQHILBERT Master equation, Liouville space.
    
    properties
        propogator
    end
    
    methods
        
        function this = RCMEdiffEqLiouville(operatorBuilder)
            [ A, H0, Xim, Xid ] = operatorBuilder.getRCMEoperators( );
            this.propogator = operatorBuilder.operatorPropogator(A, H0, Xim, Xid);
        end
        
        function drodtVector = nextDensDerivativeVector(this, densVector)
            drodtVector = this.propogator*densVector;
        end
        
    end
    
end

