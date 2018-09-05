classdef testRCMEh0 < matlab.unittest.TestCase
    %TESTRC Testuoja RC operatoriaus H0 klase
    
    properties
        params
    end
    
    methods(TestMethodSetup)
        function prepareParamsStruct(testCase)
            testCase.params.dim = 2;  
            testCase.params.J = 100; 
            testCase.params.delta_e = 100;
            testCase.params.operatoriai = {[1 0;0 0] [0 0; 0 1]};
            testCase.params.dimRC = 5;
            %Parenkami taip, kad reactDaznis ir lambda = 20
            testCase.params.alfa = 2*pi; 
            testCase.params.gamma = 10/pi;
            testCase.params.lambda = 20;
        end
    end
    
    methods (Test)
        
        %Testuoja H0 funkcija
        function testOperatorH0Dimer(testCase)
            
            bundle = RCMEcoreBundle(testCase.params);
            operatorBuilder =  RCMEoperators(bundle);
            H0 = operatorBuilder.operatorH0();
            
            [a_aukstinimo, a_zeminimo] = RCMEutil.aukstinimoZeminimo( testCase.params.dimRC );
            operatoriai = testCase.params.operatoriai;
            dimRC = testCase.params.dimRC;
            dim = testCase.params.dim;
            lambda = 20;
            reactDaznis = 20;
            H0s = [100 100;100 0];
            % Pagal formule
            H0test = kron(kron(H0s,eye(dimRC)),eye(dimRC)) -...
                lambda.*kron(kron(operatoriai{1},a_aukstinimo+a_zeminimo),eye(dimRC)) -...
                lambda.*kron(kron(operatoriai{2},eye(dimRC)),a_aukstinimo+a_zeminimo) +...
                reactDaznis.*kron(kron(eye(dim),a_aukstinimo*a_zeminimo),eye(dimRC)) + ...
                reactDaznis.*kron(kron(eye(dim),eye(dimRC)),a_aukstinimo*a_zeminimo);
                  
            testCase.verifyLessThan( abs(H0 - H0test), 0.001.*ones(dim*dimRC*dimRC) );
        end
        
        %Testuoja H0 funkcija
        function testOperatorH0SpinBoson(testCase)
            testCase.params.operatoriai = {[1 0;0 -1]};
            testCase.params.dimRC = 3;
                        
            bundle = RCMEcoreBundle(testCase.params);
            operatorBuilder =  RCMEoperators(bundle);
            H0 = operatorBuilder.operatorH0();
            
            % Ranka suskaiciuota
            H0test = [100 -20 0 100 0 0; -20 120 -20*sqrt(2) 0 100 0; 0 -20*sqrt(2) 140 0 0 100;...
                100 0 0 0 20 0; 0 100 0 20 20 20*sqrt(2); 0 0 100 0 20*sqrt(2) 40];
                  
            testCase.verifyLessThan( abs(H0 - H0test), 0.001.*ones(testCase.params.dim*...
                testCase.params.dimRC) );
        end
        
        function testOperatorH0SystemOnly(testCase)
            
            testCase.params.J = 50;
            testCase.params.dim = 2;
            testCase.params.delta_e = 100;
            bundle = RCMEcoreBundle(testCase.params);
            H0s2 = bundle.hamiltonianSystem();
            testCase.verifyEqual(H0s2, [100 50;50 0], 'dim 2');
            
            testCase.params.J = 50;
            testCase.params.dim = 3;
            testCase.params.delta_e = 100;
            bundle = RCMEcoreBundle(testCase.params);
            testH0s3 = [0 0 0;0 100 50;0 50 0;];
            H0s3 = bundle.hamiltonianSystem();
            testCase.verifyEqual(H0s3, testH0s3, 'dim 3');
        end
        
        function testErrorsWithOneOrMoreThanThreeDim(testCase)
            testCase.params.J = 50;
            testCase.params.dim = 1;
            testCase.params.delta_e = 100;
            bundle = RCMEcoreBundle(testCase.params);
            testCase.verifyError(@() bundle.hamiltonianSystem(), 'HsError:badDim', 'dim 1');
            
            testCase.params.J = 50;
            testCase.params.dim = 4;
            testCase.params.delta_e = 100;
            bundle = RCMEcoreBundle(testCase.params);
            testCase.verifyError(@() bundle.hamiltonianSystem(), 'HsError:badDim', 'dim 4');
        end
        
    end
    
end

