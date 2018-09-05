classdef testRCMEutil < matlab.unittest.TestCase
    %TESTRC Testuoja RC programa
    
    properties
    end
    
    methods (Test)
        
        %Testuoja trace funkcija
        function testTraceOutRCdimer(testCase)
            dim = 2;
            dimRC = 3;
            n = 2;
            % Sukuria densityFull 2 laikams ir duoda trace funkcijai
            densityFull = zeros(2, dim*dimRC^n, dim*dimRC^n );
            desnistyHilb = kron(kron([1 2; 3 4],[0.5 1 2; 3 0.3 0.1;2 1 0.2]),...
                [0.5 1 2; 3 0.3 0.1;2 1 0.2]);
            densityFull(1,:,:) = desnistyHilb;
            desnistyHilb = kron(kron([5 1; 7 8],[0.1 1i 2i; 3i 0.8 0.1i;2i 1i 0.1]),...
                [0.5 1 2; 3 0.3 0.1;2 1 0.2]);
            densityFull(2,:,:) = desnistyHilb;
            densitySystem = RCMEutil.traceOutRC( densityFull, dim, dimRC, n );
            
            % Kas turetu gautis
            densitySystemTest = zeros(2,dim,dim);
            densitySystemTest(1,:,:) = [1 2;3 4];
            densitySystemTest(2,:,:) = [5 1;7 8];
            
            testCase.verifyLessThan( abs(densitySystem-densitySystemTest), 0.000001*ones(size(densitySystem)) );
        end
        
        %Testuoja trace funkcija
        function testTraceOutRC(testCase)
            dim = 2;
            dimRC = 3;
            n = 1;
            % Sukuria densityFull 2 laikams ir duoda trace funkcijai
            densityFull = zeros(2, dim*dimRC, dim*dimRC );
            desnistyHilb = kron([1 2; 3 4],[0.5 1 2; 3 0.3 0.1;2 1 0.2]);
            densityFull(1,:,:) = desnistyHilb;
            desnistyHilb = kron([5 1; 7 8],[0.1 1i 2i; 3i 0.8 0.1i;2i 1i 0.1]);
            densityFull(2,:,:) = desnistyHilb;
            densitySystem = RCMEutil.traceOutRC( densityFull, dim, dimRC, n );
            
            % Kas turetu gautis
            densitySystemTest = zeros(2,dim,dim);
            densitySystemTest(1,:,:) = [1 2;3 4];
            densitySystemTest(2,:,:) = [5 1;7 8];
            
            testCase.verifyLessThan( abs(densitySystem-densitySystemTest), 0.000001*ones(size(densitySystem)) );
        end
        
        %Testuoja bandyma dalinti i per daug bloku
        function testBlockIndicesTooManyBlocks(testCase)
            
            evolTimes = 0:0.1:0.5;
            numOfBlocks = 10;
            
            testCase.verifyError(@() RCMEdynamicsHilbertTraced.getBlockIndices(evolTimes, numOfBlocks), 'blockIndices:TooManyBlocks');

        end
        
        %Testuoja bloga bloku skaiciu
        function testBlockIndicesBadBlocksNum(testCase)
            
            evolTimes = 0:0.1:0.5;
            numOfBlocks = -1;
            
            testCase.verifyError(@() RCMEdynamicsHilbertTraced.getBlockIndices(evolTimes, numOfBlocks), 'blockIndices:WrongBlocksNum');

        end
        
        %Testuoja su vienu bloku
        function testBlockIndicesOneBlock(testCase)
           
            evolTimes = 0:0.1:0.5;
            numOfBlocks = 1;
            
            shoulEvolIndices = {[1 6]};
            
            testCase.verifyEqual(RCMEdynamicsHilbertTraced.getBlockIndices(evolTimes, numOfBlocks), shoulEvolIndices);
        end
        
        %Testuoja su 2 blokais, dalinasi graziai
        function testBlockIndicesTwoBlocksDivides(testCase)
           
            evolTimes = 0:0.1:0.5;
            numOfBlocks = 2;
            
            shoulEvolIndices = {[1 3]; [4 6]};
            
            testCase.verifyEqual(RCMEdynamicsHilbertTraced.getBlockIndices(evolTimes, numOfBlocks), shoulEvolIndices);
        end
        
        %Testuoja su 2 blokais, dalinasi negraziai
        function testBlockIndicesTwoBlocksNotDivides(testCase)
           
            evolTimes = 0:0.1:0.6;
            numOfBlocks = 2;
            
            shoulEvolIndices = {[1 3]; [4 7]};
            
            testCase.verifyEqual(RCMEdynamicsHilbertTraced.getBlockIndices(evolTimes, numOfBlocks), shoulEvolIndices);
        end
        
        %Testuoja su 5 blokais, dalinasi negraziai
        function testBlockIndicesFiveBlocksNotDivides(testCase)
           
            evolTimes = 1:54;
            numOfBlocks = 5;
            
            shoulEvolIndices = {[1 10]; [11 20]; [21 30]; [31 40]; [41 54]};
            
            testCase.verifyEqual(RCMEdynamicsHilbertTraced.getBlockIndices(evolTimes, numOfBlocks), shoulEvolIndices);
        end
        
         %Testuoja Hilber Liuoville map
        function testCreateHilbertLiouvilleMap(testCase)
           hilbertLiouvilleMap = RCMEutil.createHilbertLiouvilleMap(4);
           hilbertLiouvilleMapTest = {[1 1]; [2 1]; [3 1]; [4 1]; [1 2]; [2 2]; [3 2]; [4 2]; ...
               [1 3]; [2 3]; [3 3]; [4 3]; [1 4]; [2 4]; [3 4]; [4 4] };
           
           testCase.verifyEqual(hilbertLiouvilleMap, hilbertLiouvilleMapTest);
        end
        
        %Testuoja ar transformHilbertLiuvil() duoda vektoriu su mappinimu
        %kaip HilbertLiuoville map
        function testHilbertLiouvilleMapAndT(testCase)
            matrica = [1 2 3; 4 5 6; 7 8 9];
            vektorius = RCMEutil.transformHilbertLiuvil( matrica, false );
            hilbertLiouvilleMap = RCMEutil.createHilbertLiouvilleMap(3);
            for i=1:length(vektorius)
                matricosIdxs = hilbertLiouvilleMap{i};
                testCase.verifyEqual(vektorius(i), matrica(matricosIdxs(1), matricosIdxs(2)));
            end
        end
        
        %Testuoja ar transformHilbertLiuvil() duoda vektoriu su mappinimu
        %kaip HilbertLiuoville map, kai yra laikas.
        function testHilbertLiouvilleMapAndTT(testCase)
            matricos = zeros(2,3,3);
            matricos(1,:,:) = [1 2 3; 4 5 6; 7 8 9];
            matricos(2,:,:) = [10 11 12; 13 14 15; 16 17 18];
            vektoriai = RCMEutil.transformHilbertLiuvil( matricos, true );
            hilbertLiouvilleMap = RCMEutil.createHilbertLiouvilleMap(3);
            for t=1:2
                for i=1:size(vektoriai, 2)
                    matricosIdxs = hilbertLiouvilleMap{i};
                    testCase.verifyEqual(vektoriai(t,i), matricos(t, matricosIdxs(1), matricosIdxs(2)));
                end
            end
        end
         
        %Testuoja ar abi transform funkcijos yra atvirkscios viena
        %kitai,kai yra laikas
        function testTransformInverseT(testCase)
            matricos = zeros(2,3,3);
            matricos(1,:,:) = [1 2 3; 4 5 6; 7 8 9];
            matricos(2,:,:) = [10 11 12; 13 14 15; 16 17 18];
            vektoriai = RCMEutil.transformHilbertLiuvil( matricos, true );
            matricosTest = RCMEutil.transformLiuvilHilbert(vektoriai, true);
            testCase.verifyEqual(matricos, matricosTest);
        end
        
        %Testuoja ar abi transform funkcijos yra atvirkscios viena
        %kitai,kai nera laiko
        function testTransformInverse(testCase)
            matrica = [1 2 3; 4 5 6; 7 8 9];
            vektorius = RCMEutil.transformHilbertLiuvil( matrica, false );
            matricosTest = RCMEutil.transformLiuvilHilbert(vektorius, false);
            testCase.verifyEqual(matrica, matricosTest);
        end
        
        function testCellContainsString(testCase)
            cell = {'asd' 'qwe' 'zxc'};
            
            string = 'asd';
            contains = RCMEutil.cellContainsString(cell, string);
            testCase.verifyTrue(contains, 'asd');
            
            string = 'qwe';
            contains = RCMEutil.cellContainsString(cell, string);
            testCase.verifyTrue(contains, 'qwe');
            
            string = 'zxc';
            contains = RCMEutil.cellContainsString(cell, string);
            testCase.verifyTrue(contains, 'zxc');
            
            string = 'zzz';
            contains = RCMEutil.cellContainsString(cell, string);
            testCase.verifyFalse(contains, 'zzz');
        end
        
    end
    
end

