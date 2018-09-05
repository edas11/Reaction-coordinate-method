classdef testCalculateData < matlab.unittest.TestCase
    %TESTRC Testuoja calculateData funkcija
    
    properties
       calcul 
    end
    
    %Suskaiciuoja duomenis
    methods (TestClassSetup)
        function runCalculateData(testCase)
            cd '..\funcScriptData';
            params = {[300 100 100 100 100 8 4] [150 100 100 100 100 8 4]};
            initCond = { {'21' [0 0 0;1 0 0;0 0 0]} {'31' [0 0 0;0 0 0;1 0 0]} };
            initCondMode = 'direct';
            initCondPrefix = 'IC';
            operatoriai = { [0 0 0;0 1 0;0 0 0] [0 0 0;0 0 0;0 0 1] };
            relFolder = 'test\';
            evolTimes = 0:0.1:1;
            availableModes = {'hilbert', 'hilbert-traced'};
            numberOfBlocks = 2; % Naudoja tik hilbert-traced
            currentMode = availableModes{1};
            need = {'dynamics' 'equilNoGround'};
            testCase.calcul = RCdataCalculator(params, initCond, initCondMode, initCondPrefix, ...
                relFolder, operatoriai);
            [ skipParam ] = testCase.calcul.calculateAllAndSave(need, evolTimes, currentMode, numberOfBlocks);
        end
    end
    
    % Istrina visus failus test aplankale
    methods (TestClassTeardown)
        function deleteAllTestFiles(testCase)
            delete '..\..\data\test\*.mat'
            cd '..\tests';
        end
    end
    
    methods (Test)
        
        %Testuoja, ar tesingas failu kiekis
        function testCoorectNumOfFiles(testCase)
            listing = dir('..\..\data\test');
            
            testCase.verifyEqual(length(listing), 6);
        end
        
        %Testuoja ar sukurti failai turi teisingus pavadinimus
        function testCoorectNamesExist(testCase)
            expectedNames = {'RC-(T)300(g-1)100(d)100(l)100(J)100(RCalfa)8(RCdim)4(IC)21.mat' ...
                'RC-(T)150(g-1)100(d)100(l)100(J)100(RCalfa)8(RCdim)4(IC)21.mat'...
                'RC-(T)300(g-1)100(d)100(l)100(J)100(RCalfa)8(RCdim)4(IC)31.mat'...
                'RC-(T)150(g-1)100(d)100(l)100(J)100(RCalfa)8(RCdim)4(IC)31.mat'};
            folder = '..\..\data\test\';
            for i=1:4
                testCase.verifyEqual( exist(strcat(folder,expectedNames{i}), 'file'), 2 );
            end
        end

        %Testuoja ar failai turi teisinga turini
        function testFileContent(testCase)
            folder = '..\..\data\test\';
            file = open(strcat(folder, 'RC-(T)300(g-1)100(d)100(l)100(J)100(RCalfa)8(RCdim)4(IC)31.mat'));
            t = file.t;
            densitySystem = file.densityReturned;
            equil = file.equil
            equilTest = [0.389167403870193  -0.192598949688765;-0.192598949688765   0.610832596129807];
            testCase.verifyEqual( t, 0:0.1:1, 'laikas' );
            testCase.verifyLessThan( abs(squeeze(densitySystem(9,2,1))+0.000113354374868 + 0.015063411164729i),...
                0.0001, 'densitySystem' );
            testCase.verifyLessThan( abs(equilTest-equil), 0.0001*ones(2));
        end
        
        %Testuoja ar gerai veikia skipParam
        function testSkipParam(testCase)
            params = {[300 100 100 100 100 8 4] [100 100 100 100 100 8 4]};
            initCond = { {'21' [0 0 0;1 0 0;0 0 0]}  };
            initCondMode = 'direct';
            initCondPrefix = 'IC';
            operatoriai = { [0 0 0;0 1 0;0 0 0] [0 0 0;0 0 0;0 0 1] };
            relFolder = 'test\';
            evolTimes = 0:0.1:1;
            availableModes = {'hilbert', 'hilbert-traced'};
            numberOfBlocks = 2; % Naudoja tik hilbert-traced
            currentMode = availableModes{1};
            need = {'dynamics'};
            calcul2 = RCdataCalculator(params,initCond, initCondMode, initCondPrefix, ...
                relFolder, operatoriai);
            [ skipParam ] = calcul2.calculateAllAndSave(need, evolTimes, currentMode, numberOfBlocks);
            testCase.verifyEqual( skipParam, [1 0], 'skipParam' );
            delete '..\..\data\test\RC-(T)100(g-1)100(d)100(l)100(J)100(RCalfa)8(RCdim)4(IC)21.mat'
        end
        
        function testLaunchFromWrongDirectoryGivesError(testCase)
            cd '..\';
            addpath('funcScriptData');
            testCase.verifyError(@() RCdataCalculator(), 'data:wrongDir');
            cd 'funcScriptData';
        end
        
        function testCheckNeed(testCase)
            need = {'asd' 'dynamics' 'equil'};
            isGood = testCase.calcul.checkNeed(need);
            testCase.verifyTrue(isGood, 'geras1');
            
            need = {'asd' 'equilNoGround'};
            isGood = testCase.calcul.checkNeed(need);
            testCase.verifyTrue(isGood, 'geras2');
            
            need = {'asd' 'dynamics'};
            isGood = testCase.calcul.checkNeed(need);
            testCase.verifyTrue(isGood, 'geras3');
            
            need = {'asd'};
            isGood = testCase.calcul.checkNeed(need);
            testCase.verifyFalse(isGood, 'blogas');
        end
        
        function testErrorWithBadNeed(testCase)
            params = {[300 50 50 50 100 8 4]};
            initCond = { {'21' [0 0 0;1 0 0;0 0 0]} };
            initCondPrefix = 'IC';
            initCondMode = 'direct';
            operatoriai = { [0 0 0;0 1 0;0 0 0] [0 0 0;0 0 0;0 0 1] };
            relFolder = 'test\';
            evolTimes = 0:0.1:0.5;
            hilbertMode = 'standard';
            need = {'ghhh'};
            calcul = RCdataCalculator(params, initCond, initCondMode, initCondPrefix, ...
                relFolder, operatoriai);
            testCase.verifyError(@() calcul.calculateAllAndSave(need, ...
                evolTimes, hilbertMode, false), 'data:badNeed');
        end
        
    end
    
end

