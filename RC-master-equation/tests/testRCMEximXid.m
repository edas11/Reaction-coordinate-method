classdef testRCMEximXid < matlab.unittest.TestCase
    %TESTRC Testuoja RC xim xid operatoriu klase
    
    properties
        H0
        A
        params
    end
    
    methods(TestMethodSetup)
        function prepareParamsStruct(testCase)
            testCase.H0 = [1 0 0 0; 0 2 0 0; 0 0 3 0; 0 0 0 4];
            testCase.A = {[1 0 1 0; 0 2 0 2; 3 0 3 0; 0 4 0 4] ...
                [1 0 1 0; 0 2 0 2; 3 0 3 0; 0 4 0 4]};
            testCase.params.alfa = 2*pi; 
            testCase.params.T = 5308.837458877/(0.6946*2);             
            testCase.params.operatoriai = {[1 0;0 0] [0 0; 0 1]};       
            testCase.params.dim = 4;  
            testCase.params.dimRC = 1;
        end
    end
    
    methods (Test)
        
        %Testuoja Xim ir Xid funkcija
        function testOperatorsXidXimDimer(testCase)
            %Is karto tikrines vertes, netestuoja eig ir pan.
            
            bundle = RCMEcoreBundle(testCase.params);
            operatorBuilder =  RCMEoperators(bundle);
            [ Xim, Xid ] = operatorBuilder.operatorsXimXid(testCase.H0, testCase.A);
            
            % Suskaiciuota ranka
            XimTest = [pi 0 2*coth(2)*pi 0; 0 2*pi 0 4*coth(2)*pi; 6*coth(2)*pi 0 3*pi 0;...
                0 8*coth(2)*pi 0 4*pi];
            XimTest = {XimTest XimTest};
            XidTest = [0 0 -2*pi 0; 0 0 0 -4*pi; 6*pi 0 0 0; 0 8*pi 0 0];
            XidTest = {XidTest XidTest};

            dimPilnas = 4;
            testCase.verifyLessThan( abs(Xim{1} - XimTest{1}), 0.001.*ones(dimPilnas), 'Klaida Xim 1' );
            testCase.verifyLessThan( abs(Xid{1} - XidTest{1}), 0.001.*ones(dimPilnas), 'Klaida Xid 1' );
            testCase.verifyLessThan( abs(Xim{2} - XimTest{2}), 0.001.*ones(dimPilnas), 'Klaida Xim 1' );
            testCase.verifyLessThan( abs(Xid{2} - XidTest{2}), 0.001.*ones(dimPilnas), 'Klaida Xid 1' );
        end
        
        %Testuoja Xim ir Xid funkcija
        function testOperatorsXidXim(testCase)
            %Is karto tikrines vertes, netestuoja eig ir pan.
            testCase.H0 = [1 0 0 0; 0 2 0 0; 0 0 3 0; 0 0 0 4];
            testCase.A = {[1 0 1 0; 0 2 0 2; 3 0 3 0; 0 4 0 4]};
                        
            testCase.params.operatoriai = {[1 0;0 -1]};            
            testCase.params.dim = 2;  
            testCase.params.dimRC = 2;
            
            bundle = RCMEcoreBundle(testCase.params);
            operatorBuilder =  RCMEoperators(bundle);
            [ Xim, Xid ] = operatorBuilder.operatorsXimXid(testCase.H0, testCase.A);
            
            % Suskaiciuota ranka
            XimTest = [pi 0 2*coth(2)*pi 0; 0 2*pi 0 4*coth(2)*pi; 6*coth(2)*pi 0 3*pi 0;...
                0 8*coth(2)*pi 0 4*pi];
            XidTest = [0 0 -2*pi 0; 0 0 0 -4*pi; 6*pi 0 0 0; 0 8*pi 0 0];
            
            dimPilnas = testCase.params.dim*testCase.params.dimRC;
            testCase.verifyLessThan( abs(Xim{1} - XimTest), 0.001.*ones(dimPilnas), 'Klaida Xim' );
            testCase.verifyLessThan( abs(Xid{1} - XidTest), 0.001.*ones(dimPilnas), 'Klaida Xid' );
        end
        
    end
    
end

