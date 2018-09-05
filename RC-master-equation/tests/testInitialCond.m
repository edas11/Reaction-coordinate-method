classdef testInitialCond < matlab.unittest.TestCase
    %TESTRC Testuoja RC pradiniu salygu klase
    
    properties
        params;
    end
    
    methods(TestMethodSetup)
        function prepareParamsStruct(testCase)
            DALIKLIS = 5308.837458877;
            testCase.params.dim = 2;  
            testCase.params.J = 100/DALIKLIS; 
            testCase.params.delta_e = 100/DALIKLIS;
            testCase.params.operatoriai = {[1 0;0 0] [0 0; 0 1]};
            testCase.params.dimRC = 5;
            testCase.params.alfa = 2*pi; 
            testCase.params.T = 300;
            testCase.params.gamma = 1/100;
            testCase.params.lambda = 20/DALIKLIS;
        end
    end
    
    methods (Test)
        
        function testBadMode(testCase)
            bundle = RCMEcoreBundle(testCase.params);
            operatorBuilder = RCMEoperators(bundle);
            iniInstance = RCMEinitialCondition(operatorBuilder);
            initialConditionStruct.mode = 'badMode';
            initialConditionStruct.parameter = 1;
            testCase.verifyError(@() iniInstance.get(initialConditionStruct)...
                ,'initCond:badMode' );
        end
        
        function testDirectMode(testCase)
            bundle = RCMEcoreBundle(testCase.params);
            operatorBuilder = RCMEoperators(bundle);
            iniInstance = RCMEinitialCondition(operatorBuilder);
            initialConditionStruct.mode = 'direct';
            initialConditionStruct.parameter = [1 0;0 0];
            ro = iniInstance.get(initialConditionStruct);
            
            [a_aukstinimo, a_zeminimo] = RCMEutil.aukstinimoZeminimo(bundle.getDimRC() );
            Hrc = bundle.getReactDaznis()*a_aukstinimo*a_zeminimo;
            roRC = expm(-bundle.getBeta()*Hrc);
            roRC = roRC./trace(roRC);
            rotest = kron(roRC,roRC);
            rotest = kron(initialConditionStruct.parameter, rotest);
            rotest = rotest(:);
            testCase.verifyLessThan(sum(sum(abs(ro-rotest))), 0.0001);
        end
        
        function testEquiltMode(testCase)
            testCase.params.operatoriai = {[0 0 0;0 1 0;0 0 0] [0 0 0;0 0 0;0 0 1]};
            testCase.params.dim = 3; 
            bundle = RCMEcoreBundle(testCase.params);
            operatorBuilder = RCMEoperators(bundle);
            iniInstance = RCMEinitialCondition(operatorBuilder);
            initialConditionStruct.mode = 'equilWithMiuShift';
            initialConditionStruct.parameter = [0 1 0;0 0 0;0 0 0];
            ro = iniInstance.get(initialConditionStruct);
            
            testCase.verifyTrue( sum(size(ro) == [5625 1])==2, 'dimensijos');
            testCase.verifyTrue(sum(sum((isnan(ro)))) == 0, 'ne NAN');
        end
        
    end
    
end

