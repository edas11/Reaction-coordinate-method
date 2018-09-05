classdef testDataCode < matlab.unittest.TestCase
    %TESTRC Testuoja RC duomenu funkcijas
    
    properties
    end
    
    methods (Test)
        
        %Testuoja HEOM duomenu ikelimo funkcija
        function testImportHeom(testCase)
            densityHEOM = importHEOM( '..\..\data\HEOM-dynamics-(IC)11\', ...
                '(T)150(g-1)100(d)100(l)100(J)100' );
            densityHEOMLast = squeeze(densityHEOM(end,:,:));
            
            densityTest = zeros(4,2,2);
            % Pirmi 4
            densityTest(1,:,:) = [1 0;0 0];
            densityTest(2,:,:) = [0.999593 0.000202+0.018810*1i; 0.000202-0.018810*1i 0.000407];
            densityTest(3,:,:) = [0.998509 0.000733+0.037536*1i; 0.000733-0.037536*1i 0.001491];
            densityTest(4,:,:) = [0.996735 0.001593+0.056114*1i; 0.001593-0.056114*1i 0.003265];
            % Paskutinis
            densityTestLast = [0.312145 -0.284892; -0.284892 0.687855];
            
            testCase.verifyEqual(size(densityHEOM), [4001 2 2]);
            testCase.verifyLessThan(densityHEOM(1:4,:,:)-densityTest, ones(4,2,2)*0.00000001);
            testCase.verifyLessThan(densityHEOMLast-densityTestLast, ones(2,2)*0.00000001);
        end
        
        %Testuoja createParamString
        function testcreateParamString(testCase)
            params = [300 100 10 25 50];
            names = {'(T)' '(g-1)' '(d)' '(l)' '(J)'};
            paramString = createParamString( names, params );
            paramStringTest = '(T)300(g-1)100(d)10(l)25(J)50';
            testCase.verifyEqual(paramString, paramStringTest);
        end
        
    end
    
end

