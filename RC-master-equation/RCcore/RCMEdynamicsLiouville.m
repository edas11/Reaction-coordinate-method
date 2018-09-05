classdef RCMEdynamicsLiouville < RCMEdynamics
    %RCMEMAINHILBERT  Dynamics with RCM. Liouville space.
    
    properties
    end
    
    methods
        
        function this = RCMEdynamicsLiouville( params )
            this = this@RCMEdynamics( params );
        end
        
        function diffEq = prepareDiffEq(this, operatorBuilder)
            diffEq = RCMEdiffEqLiouville(operatorBuilder);
        end
        
        function [ densityOperator ] = solve(this, evolTime, initConditionParameter)
            this.evolTime = evolTime;
            densityOperator = RCMEdensityFull(this.RCMEbundle);

            initCondition = this.initConditionBuilder.get(initConditionParameter);
            [t, densityFull] = this.solver.solveRCME(this.diffEq, initCondition, evolTime );
            densityOperator.addDensity(t, densityFull);
        end
        
    end
end

