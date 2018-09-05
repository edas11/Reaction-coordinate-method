classdef testDataLoader < matlab.unittest.TestCase
    %TESTRC Testuoja RC duomenu funkcijas
    
    properties
        dataLoader
    end
    
    methods(TestClassSetup)
        function prepareParamsStruct(testCase)
            params.J = 1;
            params.T = 2;
            params.gamma = 1/3;
            params.lambda = 4;
            params.delta_e = 5;
            fi = 6;
            testCase.dataLoader = RCFLdataLoader(params, fi);
        end
    end
    
    methods (Test)
        
        function testCreatesRCparamString(testCase)
            alfa = 8;
            dimRC = 20;
            paramString = testCase.dataLoader.createRCparamString( alfa, dimRC );
            paramStringTest = '(T)2(g-1)3(d)5(l)4(J)1(RCalfa)8(RCdim)20';
            testCase.verifyEqual(paramString, paramStringTest);
        end
        
        function testCreatesFullRCfileName(testCase)
            alfa = 8;
            dimRC = 20;
            sample = 1;
            fullName = testCase.dataLoader.createRCfileName( alfa, dimRC, sample );
            fullNameTest = 'RC-(T)2(g-1)3(d)5(l)4(J)1(RCalfa)8(RCdim)20(S)1.mat';
            testCase.verifyEqual(fullName, fullNameTest);
        end
        
        function testLoadsRCdensTequil(testCase)
            RCfolder = '../../data/test2/';
            alfa = 8;
            dimRC = 20;
            sample = 1;
            [ density, t, equil ] = testCase.dataLoader.loadRCdata( RCfolder...
                , alfa, dimRC, sample );
            testCase.verifyLessThan( abs(density(104,2,1) - 0.887793078021228 ...
                - 0.156702475634704*1i), 0.000001, 'density');
            testCase.verifyEqual(t(87), 8.6, 't');
            testCase.verifyEqual(equil, [0.6 -0.1;-0.1 0.4], 'equil');
        end
                
        function testCreatesHEOMparamString(testCase)
            paramString = testCase.dataLoader.createHEOMparamString();
            paramStringTest = '(T)2(g-1)3(d)5(l)4(J)1(fi)6';
            testCase.verifyEqual(paramString, paramStringTest);
        end
        
        function testCreatesFullHEOMfileName(testCase)
            fileName = testCase.dataLoader.createHEOMfileName();
            fileNameTest = 'HEOM-(T)2(g-1)3(d)5(l)4(J)1(fi)6.txt';
            testCase.verifyEqual(fileName, fileNameTest);
        end
        
        function testLoadsHEOMdata(testCase)
            HEOMfolder = '../../data/test2/';
           [ freqHEOM, fluorescenceHEOM ] = ...
               testCase.dataLoader.loadHEOMdata( HEOMfolder ) ;
           testCase.verifyEqual(freqHEOM(11), 12012.33, 'dazniai');
           testCase.verifyEqual(fluorescenceHEOM(11), 4.094970*10^(-3)-1.152198*1i, 'fluorescencija');
        end
        
    end
    
end

