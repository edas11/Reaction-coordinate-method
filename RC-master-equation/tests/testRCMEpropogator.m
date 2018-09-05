classdef testRCMEpropogator < matlab.unittest.TestCase
    %TESTRC Testuoja RC propogatoriaus kurima
    
    properties
    end
    
    methods (Test)
        
        %Testuoja propogator kurima
        function testOperatorPropogator(testCase)
            A = {[1 2; 3 4]};
            H0 = [2 4; 1 3];
            Xim = {[1 1; 2 1]};
            Xid = {[3 3; 3 4]};
            [ propogator ] = RCMEoperators.operatorPropogator( A, H0, Xim, Xid );
            
            % Testuojamos vertes suskaiciuotos ranka
            testCase.verifyEqual(propogator(1,1), sparse(-10), '11');
            testCase.verifyEqual(propogator(3,4), sparse(10 - 4*1i), '34');
            testCase.verifyEqual(propogator(2,3), sparse(12), '23');
        end
        
        %Testuoja propogator kurima
        function testOperatorPropogatorDimer(testCase)
            A = {[1 2; 3 4] [1 2; 3 4]};
            H0 = [2 4; 1 3];
            Xim = {[1 1; 2 1] [1 1; 2 1]};
            Xid = {[3 3; 3 4] [3 3; 3 4]};
            [ propogator ] = RCMEoperators.operatorPropogator( A, H0, Xim, Xid );
            
            % Testuojamos vertes suskaiciuotos ranka
            testCase.verifyEqual(propogator(1,1), sparse(-20), '11');
            testCase.verifyEqual(propogator(3,4), sparse(20 - 4*1i), '34');
            testCase.verifyEqual(propogator(2,3), sparse(24), '23');
        end
        
    end
    
end

