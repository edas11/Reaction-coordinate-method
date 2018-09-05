classdef RCMEutil < handle
    %RCMEUTIL Helper functions
    
    properties
    end
    
    methods(Static)
        
        function [ densitySystem ] = traceOutRC( densityFull, dim, dimRC, n )
            %TRACEOUTRC Traces out RCs from full density operator
            densitySystem = RCMEmath.traceOutEnd(dim, dimRC^n, densityFull );
        end
        
        function [ a_aukstinimo, a_zeminimo ] = aukstinimoZeminimo( dimRC )
            %AUKSTINIMOZEMINIMO Creates dimRC dimension creation and
            %anahilation operators
            a_aukstinimo = zeros( dimRC );
            a_zeminimo = zeros( dimRC );
            for i=1:dimRC-1
                a_aukstinimo( i+1, i ) = sqrt(i);
                a_zeminimo( i, i+1 ) = sqrt(i);
            end
        end
        
        
        function [ hilbertLiouvilleMap ] = createHilbertLiouvilleMap( dimPilnas )
            %HILBERTLIOUVILLEMAP Correspondance between Liuouville and
            %Hilbert dimensions
            % Pvz 1 -> 11, 2 -> 21, 3 -> 12, 4 - >22.
            hilbertLiouvilleMap = cell(dimPilnas^2, 1);
            for i=1:dimPilnas^2
                hilbertLiouvilleMap{i} = [mod(i-1,dimPilnas)+1 floor((i-1)/dimPilnas)+1];
            end
        end
        
        function [ density ] = transformHilbertLiuvil( density, hasTime )
            %TRANSFORMHILBERTLIUVIL Density operator in Hilbert space
            %transformed to Liouville space
            if hasTime
                densityTimePoints = size(density, 1);
                densityDim = size(density, 2);
                density = reshape( density, [densityTimePoints  densityDim*densityDim]);
            else
                density = density(:);
            end
        end
        
        function [ density ] = transformLiuvilHilbert( density, hasTime )
            %TRANSFORMLIUVILHILBERT Density operator in Liouville space
            %transformed to Hilbert space
            densityDims = size(density);
            if hasTime
                densityTimePoints = densityDims(1);
                dimLiouville = densityDims(2);
                density = reshape(density,[densityTimePoints sqrt(dimLiouville) sqrt(dimLiouville)]);
            else
                dimLiouville = densityDims(1);
                density = reshape(density,[sqrt(dimLiouville) sqrt(dimLiouville)]);
            end
        end
        
        function contains = cellContainsString(cell, string)
            for i = 1:length(cell)
                contains = strcmp(cell{i}, string);
                if contains
                    return;
                end
            end
        end
        
    end
    
end

