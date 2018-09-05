classdef testRCMEbundle < matlab.unittest.TestCase
    %TESTRC Testuoja RC parametru klase
    
    properties
        paramStruct
    end
    
    methods(TestMethodSetup)
        function prepareParamsStruct(testCase)
            testCase.paramStruct.dim = 1;
            testCase.paramStruct.J = 2;
            testCase.paramStruct.T = 3;
            testCase.paramStruct.gamma = 4;
            testCase.paramStruct.lambda = 5;
            testCase.paramStruct.delta_e = 6;
            testCase.paramStruct.operatoriai = 7;
            testCase.paramStruct.alfa = 8;
            testCase.paramStruct.dimRC = 9;
        end
    end
    
    methods (Test)
        
        function testCreateParamsNoComplete(testCase)
            
            bundle = RCMEcoreBundle( testCase.paramStruct );
            
            testCase.verifyEqual(bundle.getDim(), 1);
            testCase.verifyEqual(bundle.getOperatoriai(), 7);
            testCase.verifyEqual(bundle.getDimRC(), 9);
            
        end
        
        function testCreateParamsWithComplete(testCase)
            
            bundle = RCMEcoreBundle( testCase.paramStruct );
            
            testCase.verifyEqual(bundle.getBeta(), 5308.837458877/(3*0.6946));
            testCase.verifyEqual(bundle.getN(), 1);
            testCase.verifyEqual(bundle.getReactDaznis(), 8*4);
            testCase.verifyEqual(bundle.getEta(), sqrt( 2*5*8*4 ));
            
        end
        
    end
    
end

