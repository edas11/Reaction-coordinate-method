classdef testRCfluorescence < matlab.unittest.TestCase
    %TESTRC Testuoja RC duomenu funkcijas
    
    properties
        fluorInstance
    end
    
    methods(TestClassSetup)
        function prepareParamsStruct(testCase)
            params.J = 1;
            params.T = 2;
            params.gamma = 1/3;
            params.lambda = 4;
            params.delta_e = 5;
            fi = 6;
            dataLoader = RCFLdataLoader(params, fi);
            
            paramsRC.alfa = 8;
            paramsRC.dimRC = 20;
            paramsRC.energyShift = 1.5;
            paramsRC.numOfRCZeros = 1000;
            RCfolder = '../../data/test2/';
            
            testCase.fluorInstance = RCFLfluorescence(dataLoader, paramsRC, RCfolder, fi);
        end
    end
    
    methods (Test)
        
        function testLoadRCfilesForFluorescence(testCase)
            warning('off','all')
            [ density1, density2, t, equil ] = testCase.fluorInstance.loadRCdataForFluorescence( );
            
            file1 = open('../../data/test2/RC-(T)2(g-1)3(d)5(l)4(J)1(RCalfa)8(RCdim)20(S)1.mat');
            file2 = open('../../data/test2/RC-(T)2(g-1)3(d)5(l)4(J)1(RCalfa)8(RCdim)20(S)2.mat');
            density1Test = file1.densityReturned;
            density2Test = file2.densityReturned;
            equilTest = file2.equil;
            tTest = file2.t;
            
            testCase.verifyEqual(density1Test, density1, 'density1');
            testCase.verifyEqual(density2Test, density2, 'density2');
            testCase.verifyEqual(equilTest, equil, 'equil');
            testCase.verifyEqual(tTest, t, 't');
            warning('on','all')
        end
        
        function testWarningIfRCfilesDifferInT(testCase)
            t1 = [1 2 3];
            t = [1 2 4];
            testCase.verifyWarning(@() testCase.fluorInstance.checkT(t1, t),...
                'fluor:diffT');
        end
        
        function testWarningIfRCfilesDifferInEquil(testCase)
            equil1 = [1 2;3 4];
            equil2 = [1 2;3 5];
            testCase.verifyWarning(@() testCase.fluorInstance.checkEquil(equil1, equil2),...
                'fluor:diffEquil');
        end
        
        function testFindStep(testCase)
            t = [0.13 0.26 0.39];
            step = testCase.fluorInstance.findStep(t);
            testCase.verifyEqual(step, 0.13);
        end
        
        function testGetG(testCase)
            density1 = zeros(2, 3, 3);
            density2 = zeros(2, 3, 3);
            density1(1, :, :) = [1 2 3;4 5 6; 7 8 9];
            density1(2, :, :) = [11 21 31;41 51 61; 71 81 91];
            density2(1, :, :) = [12 22 32;42 52 62; 72 82 92];
            density2(2, :, :) = [13 23 33;43 53 63; 73 83 93];
            [ G ] = testCase.fluorInstance.getG( density1, density2 );
            
            Gtest{3,2}(1,1) = 0;
            Gtest{2,3}(1,1) = 0;
            Gtest{2,2}(1,1) = 1;
            Gtest{3,3}(1,1) = 1;
            
            A = [4 7; 42 72];
            B1 = [41;43];
            B2 = [71;73];
            X1 = linsolve(A, B1);
            X2 = linsolve(A, B2);
            Gtest{2,3}(2,1) = X1(2);
            Gtest{2,2}(2,1) = X1(1);
            Gtest{3,2}(2,1) = X2(1);
            Gtest{3,3}(2,1) = X2(2);
            
            testCase.verifyEqual(G, Gtest);
        end
        
        function testZeroPadding(testCase)
            G{2,3} = ones(10,1);
            G{2,2} = ones(10,1);
            G{3,2} = ones(10,1);
            G{3,3} = ones(10,1);
            numOfZeros = 90;
            Gpad = testCase.fluorInstance.zeroPadding( G, numOfZeros );
            testVal = [ones(10,1); zeros(90,1)];
            
            testCase.verifyEqual(Gpad{2,3}, testVal);
            testCase.verifyEqual(Gpad{2,2}, testVal);
            testCase.verifyEqual(Gpad{3,2}, testVal);
            testCase.verifyEqual(Gpad{3,3}, testVal);
        end
        
        function testGetDipolesFromAngle(testCase)
            miu = testCase.fluorInstance.getDipolesFromAngle();
            
            miuTest{2} = [1 0 0];
            miuTest{3} = [cos(6*pi/180) sin(6*pi/180) 0];
            
            testCase.verifyEqual(miu, miuTest);
        end
        
        function testLinearResponse(testCase)
            miu{2} = [1 0 0];
            miu{3} = [0.5 0.866 0];
            G{2,3} = [1 2];
            G{2,2} = [3 4];
            G{3,2} = [5 6];
            G{3,3} = [7 8];
            equil = [0.6 -0.1;-0.1 0.4];
            
            [ R ] = testCase.fluorInstance.getLinearResponse( G, miu, equil );

            Rtest = (1./3).*( G{2,2}.*equil(1,1) + 0.5*G{2,2}.*equil(1,2) + ...
                G{2,3}.*equil(2,1) + 0.5*G{2,3}.*equil(2,2) + 0.5*G{3,2}.*equil(1,1) + ...
                G{3,2}.*equil(1,2) + 0.5*G{3,3}.*equil(2,1) + G{3,3}.*equil(2,2));
            
            testCase.verifyLessThan( abs(R-Rtest), [0.0001 0.0001]);
        end
        
        function testMultiplyFluorescence(testCase)
            freqRC = [4; 5; 6];
            fluorescenceRC = [0.1 0.2 0.3];
            [ fluorescenceRC ] = testCase.fluorInstance.multiplyFluorescence( freqRC, fluorescenceRC );
            fluorescenceRCTest = [0.4 1 1.8];
            testCase.verifyLessThan( abs(fluorescenceRC - fluorescenceRCTest), 0.00001);
        end
        
        function testShiftFreq(testCase)
            freqRC = [-0.1 0 0.1];
            energyShift = 0.2;
            freqRC = testCase.fluorInstance.shiftFreq( freqRC, energyShift );
            freqRCTest = [0.1 0.2 0.3];
            testCase.verifyLessThan( abs(freqRC - freqRCTest), 0.00001);
        end
    end
    
end

