classdef testDensitySystem < matlab.unittest.TestCase
    %TESTRC Testuoja RC tankio operatoriaus rezultatu klase
    
    properties
        bundle
        densSystem1 = [0.7 2;2 0.3];
        densSystem2 = [0.11 -1;-1 0.89];
        densSystemTime
        t = [1 2];
    end
    
    methods(TestMethodSetup)
        function prepareParamsStruct(testCase)
            params.dim = 2;
            testCase.bundle = RCMEcoreBundle(params);
            
            testCase.densSystemTime = zeros([2 2 2]);
            testCase.densSystemTime(1,:,:) = testCase.densSystem1;
            testCase.densSystemTime(2,:,:) = testCase.densSystem2;
        end
    end
    
    methods (Test)
        
        function testCanCreateEmptyDensityObject(testCase)
            densOp = RCMEdensitySystem(testCase.bundle);
            testCase.verifyEqual(densOp.getDensitySystem(), zeros(0,2,2), 'density');
            testCase.verifyEqual(densOp.getT(), zeros(1,0), 'time');
        end
        
        function testTakesSystemMatrixGivesSystemMatrix(testCase)
            densOp = RCMEdensitySystem(testCase.bundle, testCase.t, testCase.densSystemTime);
            densityTest = densOp.getDensitySystem();
            
            testCase.verifyEqual(densityTest, testCase.densSystemTime);
        end
        
        function testTakesLargerVecotrOrMatrixGivesSystemMatrix(testCase)
            densityLarger = zeros(2,6,6);
            densRC = [0.8 1 2;1 0.05 3;2 3 0.15];
            densityLarger(1,:,:) = kron(testCase.densSystem1, densRC);
            densityLarger(2,:,:) = kron(testCase.densSystem2, densRC);
            
            densOp = RCMEdensitySystem(testCase.bundle, testCase.t, densityLarger);
            densityTest = densOp.getDensitySystem();
            testCase.verifyLessThan(densityTest-testCase.densSystemTime, 0.001, 'Matrix');
            
            densityLarger = reshape(densityLarger, [2 36]);
            densOp = RCMEdensitySystem(testCase.bundle, testCase.t, densityLarger);
            densityTest = densOp.getDensitySystem();
            testCase.verifyLessThan(densityTest-testCase.densSystemTime, 0.001, 'Vector');
        end
        
        function testErrorsIfVectorCantBeConverterToMatrix(testCase)
            densityInput = ones(2, 11);
            testCase.verifyError(@() RCMEdensitySystem(testCase.bundle, testCase.t, densityInput),'densError:cantCovertVectorMatrix' )
        end
        
        function testTakesSystemVectorGivesSystemMatrix(testCase)
            densityInput = reshape(testCase.densSystemTime, [2,4]);
            
            densOp = RCMEdensitySystem(testCase.bundle, testCase.t, densityInput);
            densityTest = densOp.getDensitySystem();
            
            testCase.verifyEqual(densityTest, testCase.densSystemTime);
        end
        
        function testAddsNewSystemMatrixWithTimeToCurrent(testCase)
            densityAdd = zeros(2,2,2);
            densityAdd(1,:,:) = [0.2 5;5 0.8];
            densityAdd(2,:,:) = [0.2 5;5 0.8];
            tt = [3 4];
            
            densOp = RCMEdensitySystem(testCase.bundle, testCase.t, testCase.densSystemTime);
            densOp.addDensity(tt, densityAdd);
            
            densityExpexted = zeros(4,2,2);
            densityExpexted(1,:,:) = testCase.densSystem1;
            densityExpexted(2,:,:) = testCase.densSystem2;
            densityExpexted(3,:,:) = [0.2 5;5 0.8];
            densityExpexted(4,:,:) = [0.2 5;5 0.8];
            tExpected = [1 2 3 4];
            
            densityActual = densOp.getDensitySystem();
            testCase.verifyEqual(densityActual, densityExpexted, 'density');
            testCase.verifyEqual(densOp.getT(), tExpected, 'time');
        end
        
        function testAddsNewSystemVectorWithTimeToCurrent(testCase)
            densityAdd = zeros(2,4);
            densityAdd(1,:,:) = [0.2; 5;5; 0.8];
            densityAdd(2,:,:) = [0.2; 5;5; 0.8];
            tt = [3 4];
            
            densOp = RCMEdensitySystem(testCase.bundle, testCase.t, testCase.densSystemTime);
            densOp.addDensity(tt, densityAdd);
            
            densityExpexted = zeros(4,2,2);
            densityExpexted(1,:,:) = testCase.densSystem1;
            densityExpexted(2,:,:) = testCase.densSystem2;
            densityExpexted(3,:,:) = [0.2 5;5 0.8];
            densityExpexted(4,:,:) = [0.2 5;5 0.8];
            tExpected = [1 2 3 4];
            
            densityActual = densOp.getDensitySystem();
            testCase.verifyEqual(densityActual, densityExpexted, 'density');
            testCase.verifyEqual(densOp.getT(), tExpected, 'time');
        end
        
        function testRejectsBadDensityMatrix(testCase)
            densitySample = zeros(2,13,13);
            
            testCase.verifyError(@() RCMEdensitySystem(testCase.bundle, testCase.t, densitySample)...
                ,'densError:badDensity' );
        end
        
        function testRejectsBadDensityVector(testCase)
            densitySample = zeros(2,25);
            
            testCase.verifyError(@() RCMEdensitySystem(testCase.bundle, testCase.t, densitySample)...
                ,'densError:badDensity' );
        end
        
        function testRejectsLargerThan3Ddensity(testCase)
            densitySample = zeros(5,5,5,5);
            
            testCase.verifyError(@() RCMEdensitySystem(testCase.bundle, testCase.t, densitySample)...
                ,'densError:badDensityLargeD' );
        end
        
        function testRejectsNot2DTime(testCase)
            tSample = zeros(2,3,3);
            testCase.verifyError(@() RCMEdensitySystem(testCase.bundle, tSample, testCase.densSystemTime)...
                ,'densError:badTimeDim' );
        end
        
        function testRejectsNotVectorTime(testCase)
            tSample = zeros(2,3);
            testCase.verifyError(@() RCMEdensitySystem(testCase.bundle, tSample, testCase.densSystemTime)...
                ,'densError:badTimeDim' );
        end
        
        function testTakesTrowVectorGivesTrowVector(testCase)
            tSample = [1 2 3 4];
            densOp = RCMEdensitySystem(testCase.bundle, tSample, testCase.densSystemTime);
            testCase.verifyEqual(densOp.getT(), tSample);
        end
        
        function testTakesTcolumnVectorGivesTrowVector(testCase)
            tSample = [1; 2; 3; 4];
            densOp = RCMEdensitySystem(testCase.bundle, tSample, testCase.densSystemTime);
            testCase.verifyEqual(densOp.getT(), tSample');
        end
        
    end
    
end

