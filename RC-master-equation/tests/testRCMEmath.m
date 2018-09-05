classdef testRCMEmath < matlab.unittest.TestCase
    %TESTRC Testuoja RC matematines funkcijas
    
    properties
    end
    
    methods (Test)
        
        % Testuoja aukstinimo ir zeminimo funkcija
        function testAukstinimoZeminimo(testCase)
            [ a_aukstinimo, a_zeminimo ] = RCMEutil.aukstinimoZeminimo( 5 );
            aukstinimo_test = [0 0 0 0 0; 1 0 0 0 0; 0 sqrt(2) 0 0 0; ...
                0 0 sqrt(3) 0 0; 0 0 0 2 0];
            zeminimo_test = aukstinimo_test';
            testCase.verifyEqual(a_aukstinimo, aukstinimo_test);
            testCase.verifyEqual(a_zeminimo, zeminimo_test);
        end
        
        %Testuoja trace begining funkcija
        function testTraceOutBegining(testCase)
            dims.dim = 2;
            dims.dimRC = 3;

            % Sukuria densityFull 2 laikams ir duoda trace funkcijai
            densityFull = zeros(2, dims.dim*dims.dimRC, dims.dim*dims.dimRC);
            desnistyHilb = kron([0.4 2; 3 0.6],[0.5 1 2; 3 0.3 0.1;2 1 0.2]);
            densityFull(1,:,:) = desnistyHilb;
            desnistyHilb = kron([0.88 1; 7 0.12],[0.1 1i 2i; 3i 0.8 0.1i;2i 1i 0.1]);
            densityFull(2,:,:) = desnistyHilb;
            densitySystem = RCMEmath.traceOutBegining( 3, 2, densityFull );

            % Kas turetu gautis
            densitySystemTest = zeros(2,dims.dimRC,dims.dimRC);
            densitySystemTest(1,:,:) = [0.5 1 2; 3 0.3 0.1;2 1 0.2];
            densitySystemTest(2,:,:) = [0.1 1i 2i; 3i 0.8 0.1i;2i 1i 0.1];
            
            testCase.verifyLessThan( abs(densitySystem-densitySystemTest), 0.000001*ones(size(densitySystem)) );
        end
        
        %Testuoja trace begining funkcija
        function testTraceOutBeginingNoTime(testCase)
            dims.dim = 2;
            dims.dimRC = 3;

            % Sukuria densityFull 2 laikams ir duoda trace funkcijai
            densityFull = kron([0.4 2; 3 0.6],[0.5 1 2; 3 0.3 0.1;2 1 0.2]);
            densitySystem = RCMEmath.traceOutBegining( 3, 2, densityFull );

            % Kas turetu gautis
            densitySystemTest = [0.5 1 2; 3 0.3 0.1;2 1 0.2];
            
            testCase.verifyLessThan( abs(densitySystem-densitySystemTest), 0.000001*ones(size(densitySystem)) );
        end
        
        %Testuoja trace end funkcija
        function testTraceOutEnd(testCase)
            dims.dim = 2;
            dims.dimRC = 3;

            % Sukuria densityFull 2 laikams ir duoda trace funkcijai
            densityFull = zeros(2, dims.dim*dims.dimRC, dims.dim*dims.dimRC);
            desnistyHilb = kron([0.4 2; 3 0.6],[0.5 1 2; 3 0.3 0.1;2 1 0.2]);
            densityFull(1,:,:) = desnistyHilb;
            desnistyHilb = kron([0.88 1; 7 0.12],[0.1 1i 2i; 3i 0.8 0.1i;2i 1i 0.1]);
            densityFull(2,:,:) = desnistyHilb;
            densitySystem = RCMEmath.traceOutEnd( 2, 3, densityFull );

            % Kas turetu gautis
            densitySystemTest = zeros(2,dims.dim,dims.dim);
            densitySystemTest(1,:,:) = [0.4 2; 3 0.6];
            densitySystemTest(2,:,:) = [0.88 1; 7 0.12];
            
            testCase.verifyLessThan( abs(densitySystem-densitySystemTest), 0.000001*ones(size(densitySystem)) );
        end
        
        %Testuoja trace end funkcija
        function testTraceOutEndNoTime(testCase)
            dims.dim = 2;
            dims.dimRC = 3;

            % Sukuria densityFull 2 laikams ir duoda trace funkcijai
            densityFull = kron([0.4 2; 3 0.6],[0.5 1 2; 3 0.3 0.1;2 1 0.2]);
            densitySystem = RCMEmath.traceOutEnd( 2, 3, densityFull );

            % Kas turetu gautis
            densitySystemTest = [0.4 2; 3 0.6];
            
            testCase.verifyLessThan( abs(densitySystem-densitySystemTest), 0.000001*ones(size(densitySystem)) );
        end
        
    end
    
end

