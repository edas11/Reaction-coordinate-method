classdef RCMEdensity < handle
    %RCMEDENSITYOPERATOR Density operator results (with time)
    
    properties
        t = zeros(1,0);
    end
    
    methods
        
        function addDensity(this, t, density)
            densityMatrix = this.ensureDensityIsMatrix(density);
            tVector = this.ensureTisCorrectVector(t);
            
            this.checkConsistencyOfMatrixDimensions(densityMatrix);
            
            this.addDensityMatrix(densityMatrix);
            this.addT(tVector);
        end
        
        function addT(this,tToAdd)
            this.t = cat(2,this.t, tToAdd);
        end
        
        function densityMatrix = densityVectorToMatrix(this, densityVector)
            timeDim = size(densityVector, 1);
            matrixDim = sqrt(size(densityVector, 2));
            if round(matrixDim)~=matrixDim
                error('densError:cantCovertVectorMatrix', 'Density gali buti vektorius tik jei is jo dim galima istraukti sakni');
            end
            densityMatrix = reshape(densityVector, [timeDim matrixDim matrixDim]);
        end
        
        function densityMatrix = ensureDensityIsMatrix(this, density)
            dims = size(density);
            numOfDims = length(dims);
            if numOfDims == 2
                densityMatrix = this.densityVectorToMatrix(density);
                return
            elseif numOfDims~=3
                error('densError:badDensityLargeD', 'Density turi tureti 2 arba 3 dimensijas');
            end
            densityMatrix = density;
        end
        
        function tVector = ensureTisCorrectVector(this, t)
            dims = size(t);
            numOfDims = length(dims);
            if numOfDims~=2
                error('densError:badTimeDim', 'Laikas turi buti vektorius');
            end
            if dims(1)~=1 && dims(2)~=1
                error('densError:badTimeDim', 'Laikas turi buti vektorius');
            end
            tVector = reshape(t, [1, length(t)]);
        end
        
        function checkConsistencyOfMatrixDimensions(this, density)
            dims = size(density);
            possibleDim1 = dims(2);
            possibleDim2 = dims(3);
            if this.badMatrixDims(possibleDim1, possibleDim2)
                error('densError:badDensity', 'Density turi blogas dimensijas');
            end
        end
        
        function t = getT(this)
            t = this.t;
        end
        
    end
    
    methods(Abstract)
        densitySystem = getDensitySystem(this)
        addDensityMatrix(this, t, density)
        isBadDims = badMatrixDims(this, possibleDim1, possibleDim2)
    end
    
end

