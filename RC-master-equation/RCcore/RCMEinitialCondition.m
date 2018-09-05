classdef RCMEinitialCondition < handle
    %RCMEINITIALCONDITION System + RCs initial condition.
    
    properties (Constant, Access = private)
        % initialConditionStruct.parameter is initial system condition and
        % RCs initial cond. is Boltzmann equilibrium according to their
        % hamiltonian. Full init. cond. is product of both
        DIRECT = 'direct'
        
        % initialConditionStruct.parameter is matrix miu and full initial
        % condition is equilibrium density operator * miu. Equilibrium is
        % Boltzmann equilibrium with respect to system + RCs hamiltonian
        EQUILMIU = 'equilWithMiuShift'
    end
    
    properties (Access = private)
        operatorBuilder
    end
    
    methods
        
        function this = RCMEinitialCondition(operatorBuilder)
            this.operatorBuilder = operatorBuilder;
        end
        
        function [ initCondition ] = get(this, initialConditionStruct)
            initMode = initialConditionStruct.mode;
            initParameter = initialConditionStruct.parameter;
            
            if strcmp(initMode, this.DIRECT)
                [ initCondition ] = this.getDirect( initParameter );
                return
            end
            if strcmp(initMode, this.EQUILMIU)
                [ initCondition ] = this.getEquilMiuShifted( initParameter );
                return
            end
            error('initCond:badMode', 'Negeras pradiniu salygu mode');
        end
        
    end
    
    methods (Access=private)
        
        function [ initCondition ] = getDirect( this, pradSystemDensity )
            pradSystemDensity
            pradRC = this.operatorBuilder.getRCsBoltzmannEquil();
            
            initCondition = this.operatorBuilder.combineSystemRCsOperators(...
                pradSystemDensity, pradRC);
            initCondition = initCondition(:);
        end
        
        function [ initCondition ] = getEquilMiuShifted( this, miuSystemBasis )
            systemAndRCsEquilibrium = this.operatorBuilder.getH0BoltzmannEquil(true);
            miu = this.operatorBuilder.getSystemOperatorFullBasis(miuSystemBasis);
            
            initCondition = systemAndRCsEquilibrium*miu;
            initCondition = initCondition(:);
        end
        
    end
    
end

