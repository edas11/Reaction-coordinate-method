classdef RCMEdynamics < handle
    %RCMEMAINHILBERT Dynamics with RCM.
    
    properties
        evolTime % Evolution length, for ode45
        diffEq % Diff. equation
        RCMEbundle
        operatorBuilder
        initConditionBuilder
        solver
    end
    
    methods
        
        function this = RCMEdynamics( params )
            this.RCMEbundle = RCMEcoreBundle(params);
            this.operatorBuilder = RCMEoperators(this.RCMEbundle);
            this.initConditionBuilder = RCMEinitialCondition(this.operatorBuilder);
            this.diffEq = this.prepareDiffEq(this.operatorBuilder);
            this.solver = RCMEsolver();
        end
        
    end
    
    methods(Abstract)
        [ diffEq ] = prepareDiffEq(this, operatorBuilder)
        [ densityOperator ] = solve(this, initConditionParameter)
    end
    
end

