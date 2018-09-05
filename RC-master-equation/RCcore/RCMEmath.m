classdef RCMEmath
    %RCMEMATH Some math functions for RC program
    
    properties
    end
    
    methods( Static)
        
        function [ densityTraced ] = traceOutBegining( dimSumazintas, dimToTrace, density )
            %TRACEOUTBEGINING If desnity = A*B, this function returns densityTraced=B.
            %If 3 dimensions, first is not touched.
            if dimToTrace==0
                densityTraced = density;
                return;
            end
            numOfDims = length(size(density));
            if numOfDims==2
                densityTimePoints = 1;
                density = reshape(density, [1 size(density)]);
            elseif numOfDims==3
                densityTimePoints = size(density, 1);
            end
            densityTraced = zeros(densityTimePoints, dimSumazintas, dimSumazintas);
            for row = 1:dimSumazintas
                for column = 1:dimSumazintas
                    % trace RC
                    for indexRC = 0:dimSumazintas:dimSumazintas*(dimToTrace-1)
                        densityTraced(:, row, column) = densityTraced(:, row, column) +...
                            density(:, row + indexRC, column + indexRC);
                    end
                end
            end
            if numOfDims==2
                newDims = size(densityTraced);
                densityTraced = reshape(densityTraced, newDims(2:end));
            end
        end
        
        function [ densityTraced ] = traceOutEnd( dimSumazintas, dimToTrace, density )
            %TRACEOUTBEGINING If desnity = A*B, this function returns densityTraced=A.
            %If 3 dimensions, first is not touched.
            if dimToTrace==0
                densityTraced = density;
                return;
            end
            numOfDims = length(size(density));
            if numOfDims==2
                densityTimePoints = 1;
                density = reshape(density, [1 size(density)]);
            elseif numOfDims==3
                densityTimePoints = size(density, 1);
            end
            densityTraced = zeros(densityTimePoints, dimSumazintas, dimSumazintas);
            for row = 1:dimSumazintas
                for column = 1:dimSumazintas
                    % trace RC
                    for indexRC = 1:dimToTrace
                        densityTraced(:, row, column) = densityTraced(:, row, column) +...
                            density(:, (row-1)*dimToTrace + indexRC, (column-1)*dimToTrace + indexRC);
                    end
                end
            end
            if numOfDims==2
                newDims = size(densityTraced);
                densityTraced = reshape(densityTraced, newDims(2:end));
            end
        end
        
        function [ product ] = kron_AIBI( system, RC, RCpos, lngth )
            %TENSORPRD1SYST1RC Tensor product. First matrix is system, RCpos
            %matrix is RC, other are RC dimension identity matrices. In total lngth
            %terms. A*I*B*I...
            dimRC = length(RC);
            product = system;
            for i=2:lngth
                if i==RCpos
                    product = kron(product, RC);
                else
                    product = kron(product, eye(dimRC));
                end
            end
        end
        
        function [ product ] = kron_ABBB( system, RC, lngth )
            %TENSORPRD1SYSTMANYRC Tensor product. First matric is system.
            %other are RC. In total n matrices. A*B*B*B*
            product = system;
            for i=2:lngth
                product = kron(product, RC);
            end
        end
        
        function [ product ] = kron_BBB( B, lngth )
            %TENSORPRD1SYSTMANYRC Tensor product. All matrices are RC, B*B*B*...
            product = B;
            for i=2:lngth
                product = kron(product, B);
            end
        end
        
        function [ C ] = commutator( A, B )
            C = A*B - B*A;
        end
        
        function [ C ] = antiCommutator( A, B )
            C = A*B + B*A;
        end
        
    end
    
end

