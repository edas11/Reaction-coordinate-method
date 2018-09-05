classdef RCMEdynamicsHilbert < RCMEdynamics
    %RCMEMAINHILBERT Dynamics with RCM. Hilbert space.
    
    properties
    end
    
    methods
        
        function this = RCMEdynamicsHilbert( params )
            this = this@RCMEdynamics( params );
        end
        
        function diffEq = prepareDiffEq(this, operatorBuilder)
            diffEq = RCMEdiffEqHilbert(operatorBuilder);
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

