classdef RCMEdynamicsHilbertTraced < RCMEdynamics
    %RCMEMAINHILBERT Dynamics with RCM. Hilbert space with repeated tracing through RCs.
    
    properties
        numberOfBlocks
        lastKnownDensityFull
        initConditionParameter
    end
    
    methods
        
        function this = RCMEdynamicsHilbertTraced( params, numberOfBlocks )
            this = this@RCMEdynamics( params );
            
            this.numberOfBlocks = numberOfBlocks;
        end
        
        function diffEq = prepareDiffEq(this, operatorBuilder)
            diffEq = RCMEdiffEqHilbert(operatorBuilder);
        end
        
        function [ densityOperator ] = solve(this, evolTime, initConditionParameter)
            this.initConditionParameter = initConditionParameter;
            densityOperator = RCMEdensitySystem(this.RCMEbundle);
            evolInd = this.getBlockIndices( evolTime, this.numberOfBlocks );
            
            for cycle = 1:this.numberOfBlocks
                currentEvolTime = evolTime(evolInd{cycle}(1):evolInd{cycle}(2));
                
                initCondition = this.prepareInitCondition();
                [t, densityFull] = this.solver.solveRCME(this.diffEq, initCondition, currentEvolTime );
                densityOperator.addDensity(t, densityFull);
                
                this.lastKnownDensityFull = densityFull();
            end
        end
        
        function initCondition = prepareInitCondition(this)
            if isempty(this.lastKnownDensityFull)
                initCondition = this.initConditionBuilder.get(this.initConditionParameter);
            else
                initCondition = this.lastKnownDensityFull(end,:,:);
            end
            
        end
        
    end
    
    methods(Static)
        
        function [ evolIndices ] = getBlockIndices( evolTimes, numOfBlocks )
            %GETBLOCKINDICES Begining and end times of time blocks
            
            if numOfBlocks<1
                error('blockIndices:WrongBlocksNum', 'Blogas bloku skaicius');
            end
            
            if numOfBlocks>length(evolTimes)
                error('blockIndices:TooManyBlocks', 'Per daug bloku');
            end
            
            totalLength = length(evolTimes);
            dividedLength = floor(totalLength/numOfBlocks);
            dividedMod = mod(totalLength, numOfBlocks);
            
            evolIndices = cell(numOfBlocks, 1);
            for blockIdx = 1:numOfBlocks
                evolIndices{blockIdx} = [dividedLength*(blockIdx-1)+1 dividedLength*blockIdx];
            end
            evolIndices{end}(2) = evolIndices{end}(2) + dividedMod;
            
        end
        
    end
end

