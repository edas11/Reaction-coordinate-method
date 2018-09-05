classdef RCMEsolver < handle
    %RCMESOLVER Solves RCME
    
    properties
    end
    
    methods
        
        function [ t, densityFull ] = solveRCME( this, diffEq, initCondition, evolTimes )
            %SOLVERCME Solves RCME with ode45 iki evolTime
            %options = odeset('MaxStep',1e-2, 'RelTol',1e-10,'AbsTol',1e-10);
            [t, densityFull] = ode45(@(t,y) this.RCME(t, y, diffEq), evolTimes, initCondition);
            
        end
        
        function [ drodtVector ] = RCME( this, t, densVector, diffEq )
            % Status variable
            persistent index;
            if isempty(index)
                index = 0;
            end
            
            % RCME
            drodtVector = diffEq.nextDensDerivativeVector(densVector);
            
            % Prints current calculation time, if status variable allows it
            if mod(index,1000) == 0
                t
            end
            %  Updates status
            index = index + 1;
        end
        
    end
    
end

