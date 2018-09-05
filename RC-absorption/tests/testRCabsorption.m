classdef testRCabsorption < matlab.unittest.TestCase
    %TESTRC Testuoja RC sugerties klase
    
    properties
        absInstance
    end
    
    methods(TestClassSetup)
        function prepareParamsStruct(testCase)
            params.J = 1;
            params.T = 2;
            params.gamma = 1/3;
            params.lambda = 4;
            params.delta_e = 5;
            fi = 6;
            dataLoader = RCABdataLoader(params, fi);
            
            paramsRC.alfa = 8;
            paramsRC.dimRC = 20;
            paramsRC.energyShift = 1.5;
            paramsRC.numOfRCZeros = 1000;
            RCfolder = '../../data/test2/';
            
            testCase.absInstance = RCABabsorption(dataLoader, paramsRC, RCfolder, fi);
        end
    end
    
    methods (Test)
        
        function testLoadRCfilesForFAbsorption(testCase)
            warning('off','all')
            [ density21, density31, t ] = testCase.absInstance.loadRCdataForAbsorption( );
            
            file21 = open('../../data/test2/RC-(T)2(g-1)3(d)5(l)4(J)1(RCalfa)8(RCdim)20(IC)21.mat');
            file31 = open('../../data/test2/RC-(T)2(g-1)3(d)5(l)4(J)1(RCalfa)8(RCdim)20(IC)31.mat');
            density21Test = file21.densityReturned;
            density31Test = file31.densityReturned;
            tTest = file31.t;
            
            testCase.verifyEqual(density21Test, density21, 'density1');
            testCase.verifyEqual(density31Test, density31, 'density2');
            testCase.verifyEqual(tTest, t, 't');
            warning('on','all')
        end
        
        function testWarningIfRCfilesDifferInT(testCase)
            t1 = [1 2 3];
            t = [1 2 4];
            testCase.verifyWarning(@() testCase.absInstance.checkT(t1, t),...
                'fluor:diffT');
        end
        
        function testFindStep(testCase)
            t = [0.13 0.26 0.39];
            step = testCase.absInstance.findStep(t);
            testCase.verifyEqual(step, 0.13);
        end
        
        function testGetG(testCase)
            density21 = zeros(2, 3, 3);
            density31 = zeros(2, 3, 3);
            density21(1, :, :) = [1 2 3;4 5 6; 7 8 9];
            density21(2, :, :) = [11 21 31;41 51 61; 71 81 91];
            density31(1, :, :) = [12 22 32;42 52 62; 72 82 92];
            density31(2, :, :) = [13 23 33;43 53 63; 73 83 93];
            [ G ] = testCase.absInstance.getG( density21, density31 );
            
            Gtest{3,2} = [7; 71];
            Gtest{2,2} = [4; 41];
            Gtest{2,3} = [42; 43];
            Gtest{3,3} = [72; 73];
            
            testCase.verifyEqual(G, Gtest);
        end
        
        function testZeroPadding(testCase)
            G{2,3} = ones(10,1);
            G{2,2} = ones(10,1);
            G{3,2} = ones(10,1);
            G{3,3} = ones(10,1);
            numOfZeros = 90;
            Gpad = testCase.absInstance.zeroPadding( G, numOfZeros );
            testVal = [ones(10,1); zeros(90,1)];
            
            testCase.verifyEqual(Gpad{2,3}, testVal);
            testCase.verifyEqual(Gpad{2,2}, testVal);
            testCase.verifyEqual(Gpad{3,2}, testVal);
            testCase.verifyEqual(Gpad{3,3}, testVal);
        end
        
        function testGetDipolesFromAngle(testCase)
            miu = testCase.absInstance.getDipolesFromAngle();
            
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
            
            [ R ] = testCase.absInstance.getLinearResponse( G, miu );

            Rtest = (1./3).*( G{2,2} + 0.999956*G{3,3} + 0.5.*G{2,3} + 0.5.*G{3,2});
            
            testCase.verifyLessThan( abs(R-Rtest), [0.0001 0.0001]);
        end
        
        function testMultiplyAbsorption(testCase)
            freqRC = [4; 5; 6];
            absorptionRC = [0.1 0.2 0.3];
            [ absorptionRC ] = testCase.absInstance.multiplyAbsorption( freqRC, absorptionRC );
            absorptionRCTest = [0.4 1 1.8];
            testCase.verifyLessThan( abs(absorptionRC - absorptionRCTest), 0.00001);
        end
        
        function testShiftFreq(testCase)
            freqRC = [-0.1 0 0.1];
            energyShift = 0.2;
            freqRC = testCase.absInstance.shiftFreq( freqRC, energyShift );
            freqRCTest = [0.1 0.2 0.3];
            testCase.verifyLessThan( abs(freqRC - freqRCTest), 0.00001);
        end
    end
    
end

