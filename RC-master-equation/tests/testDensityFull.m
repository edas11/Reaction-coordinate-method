classdef testDensityFull < matlab.unittest.TestCase
    %TESTRC Testuoja RC tankio operatoriaus rezultatu klase
    
    properties
        bundle
        densSystem = [0.7 2;2 0.3];
        densRC1 = [0.1 1 2;1 0.5 4;2 4 0.4];
        densRC2 = [0.25 11 21;11 0.3 41;21 41 0.45];
        t = [1; 2];
        densityFull
    end
    
    methods(TestMethodSetup)
        function prepareParamsStruct(testCase)
            params.dim = 2;
            params.operatoriai = {[1 0;0 0] [0 0; 0 1]};
            params.dimRC = 3;
            testCase.bundle = RCMEcoreBundle(params);
            
            density = kron(testCase.densSystem, testCase.densRC1);
            density = kron(density, testCase.densRC2);
            
            testCase.densityFull = zeros([2 size(density)]);
            testCase.densityFull(1,:,:) = density;
            testCase.densityFull(2,:,:) = density;
        end
    end
    
    methods (Test)

        function testTakesFullMatrixGivesFullMatrix(testCase)
            densOp = RCMEdensityFull(testCase.bundle, testCase.t, testCase.densityFull);
            densityTest = densOp.getDensityFull();
            
            testCase.verifyEqual(densityTest, testCase.densityFull);
        end
        
        function testTakesFullMatrixGivesSystem(testCase)
            densOp = RCMEdensityFull(testCase.bundle, testCase.t, testCase.densityFull);
            densityTest = densOp.getDensitySystem();
            
            densityExpexted = zeros(2,2,2);
            densityExpexted(1,:,:) = testCase.densSystem;
            densityExpexted(2,:,:) = testCase.densSystem;
            
            testCase.verifyEqual(densityTest, densityExpexted);
        end
        
        function testTakesFullVectorGivesFullMatrix(testCase)
            densitySample = reshape(testCase.densityFull, [2 324]);
            densOp = RCMEdensityFull(testCase.bundle, testCase.t, densitySample);
            densityTest = densOp.getDensityFull();
            
            testCase.verifyEqual(densityTest, testCase.densityFull);
        end
        
        function testRejectsBadDensityMatrix(testCase)
            densitySample = zeros(2,13,13);
            
            testCase.verifyError(@() RCMEdensityFull(testCase.bundle, testCase.t, densitySample)...
                ,'densError:badDensity' );
        end
        
        function testRejectsBadDensityVecotr(testCase)
            densitySample = zeros(2,16);
            
            testCase.verifyError(@() RCMEdensityFull(testCase.bundle, testCase.t, densitySample)...
                ,'densError:badDensity' );
        end
        
        function testRejectsLargerThan3Ddensity(testCase)
            densitySample = zeros(5,5,5,5);
            
            testCase.verifyError(@() RCMEdensityFull(testCase.bundle, testCase.t, densitySample)...
                ,'densError:badDensityLargeD' );
        end
        
    end
    
end

