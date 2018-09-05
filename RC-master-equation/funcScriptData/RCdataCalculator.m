classdef RCdataCalculator < handle
    %RCDATACALCULATOR Calculates multiple times and writes to disk.
    %Calculates for every top level cell array element.
    
    properties (Constant, Access = private)
        names = {'(T)' '(g-1)' '(d)' '(l)' '(J)' '(RCalfa)' '(RCdim)'}; % Parameter names
        DALIKLIS = 5308.837458877; % Conversion between cm^-1 and fs^-1 (cm^-1/fs^-1)
    end
    
    properties (Access = private)
        inputParams % cell array with vector elements whose values correspond to
                    % parameters in names
        initCond    % cell array that hold cell arrays which have initial condition
                    % name and initial condition parameter.
        initCondMode
        initCondPrefix  % Initial condition prefix in file name
        folder      % Where to save
        operatoriai % System - RCs (bath) interraction operators
    end
    
    methods
        
        function this = RCdataCalculator(inputParams, initCond, initCondMode, initCondPrefix,...
                relFolder, operatoriai)
            addpath('../RCcore');
            addpath('../');
            isInWrongDirectory = isempty(strfind(pwd,'\RC-master-equation\funcScriptData'));
            if isInWrongDirectory
                error('data:wrongDir','Reikia paleisti is funcScriptData direktorijos');
            end
            
            this.inputParams = inputParams;
            this.initCond = initCond;
            this.initCondMode = initCondMode;
            this.initCondPrefix = initCondPrefix;
            this.folder = strcat('..\..\data\', relFolder);
            this.operatoriai = operatoriai;
        end
        
        % Calculates
        function [ skipParam ] = calculateAllAndSave( this, need, evolTimes, currentMode, numberOfBlocks )

            numOfParams = length(this.inputParams);
            numOfinit = length(this.initCond);

            % Some parameters can be skiped. Here skiped ones will be
            % marked.
            skipParam = zeros(numOfinit, numOfParams);

            for initCondIndex = 1:numOfinit
                currentInit = this.initCond{initCondIndex};
                for paramIndex = 1:numOfParams
                    currentInputParams = this.inputParams{paramIndex};

                    % Prepares
                    params = this.prepareParametersForCalculation( currentInputParams, currentInit{2} );

                    fullPath = this.createFullPathName( currentInputParams, currentInit{1});
                    
                    % For one set of parameters
                    initCond.mode = this.initCondMode;
                    initCond.parameter = currentInit{2};
                    success = this.calculateOneAndSave(need, params, initCond, fullPath,...
                        evolTimes, currentMode, numberOfBlocks);
                    
                    %Mark if failed
                    if ~success
                        skipParam(initCondIndex, paramIndex) = 1;
                    end

                    % Disgnostics
                    this.getCompletePercent(initCondIndex, paramIndex)
                end
            end
        end
        
        % Calculates for one set of parametres and writes to disk,
        function success = calculateOneAndSave(this, need, params, initCond, fullPath, ...
                 evolTimes, currentMode, numberOfBlocks)
            
            if exist(fullPath, 'file')
                success = false;
                return;
            end

            if ~this.checkNeed(need)
                error('data:badNeed','Need turi tureti bent viena is: dynamics, equilNoGround');
            end
            
            % Output'as is writen to struct
            output = struct();
            
         	% Calculates
            if RCMEutil.cellContainsString(need, 'dynamics')
                [densityReturned, t, evolExeTime] = this.calcRCDynamics(params, initCond, ...
                     evolTimes, currentMode, numberOfBlocks);
                output.densityReturned = densityReturned;
                output.t = t;
                output.evolExeTime = evolExeTime;
            end
            
            if RCMEutil.cellContainsString(need, 'equilNoGround')
                [equil] = this.calcRCequilNoGround(params);
                output.equil = equil;
            end
             
            % Saves
            save( fullPath, '-struct', 'output');
            
            success = true;
        end
        
        % Prepares
        function params = prepareParametersForCalculation(this, inputParams, initCond )
            params.dim = length(initCond);
            params.J = inputParams(5)/this.DALIKLIS;
            params.T = inputParams(1);
            params.gamma = 1/inputParams(2);
            params.lambda = inputParams(4)/this.DALIKLIS;
            params.delta_e = inputParams(3)/this.DALIKLIS;
            params.alfa = inputParams(6);
            params.dimRC = inputParams(7);
            params.operatoriai = this.operatoriai;
        end
        
        % RC dynamics
        function  [densityReturned, t, evolExeTime] = calcRCDynamics(this, params, initCond, ...
                 evolTimes, currentMode, numOfBlocks)
            evol = tic;
            if strcmp(currentMode,'hilbert')
                RCMEinstance = RCMEdynamicsHilbert(params);
            elseif strcmp(currentMode,'hilbert-traced')
                RCMEinstance = RCMEdynamicsHilbertTraced(params, numOfBlocks);
            end
            [densityOperator] = RCMEinstance.solve( evolTimes, initCond );
            densityReturned = densityOperator.getDensitySystem();
            t = densityOperator.getT();
            evolExeTime = toc(evol);
        end
        
        % RC equilibrium
        function equil = calcRCequilNoGround(this, params)
            bundle = RCMEcoreBundle(params);
            operatorBuilder = RCMEoperators(bundle);
            equil = operatorBuilder.getH0BoltzmannEquil(true);
            dimToTrace = size(equil,2)/params.dim;
            equil = RCMEmath.traceOutEnd( params.dim, dimToTrace, equil );
            equil = equil(2:end,2:end);
        end
        
        % What percentage of calculations is done
        function completePercent = getCompletePercent(this, initCondIndex, paramIndex)
            numOfParams = length(this.inputParams);
            numOfinit = length(this.initCond);
            completePercent = strcat(num2str( 100*((initCondIndex-1)*numOfParams + paramIndex)/...
                    (numOfParams*numOfinit) ), '%');
        end
        
        function [fullPath] = createFullPathName(this, params, initName)
                paramString = this.createParamString(params);
                fullPath = strcat(this.folder, 'RC-',paramString,...
                    '(', this.initCondPrefix,')', initName, '.mat');
        end
        
        % File name according to parameters
        function [ paramString ] = createParamString( this, params )
            paramString = '';
            for i=1:length(this.names)
                paramString = strcat(paramString, this.names{i}, num2str(params(i)));
            end
        end
        
        function isGood = checkNeed(this, need)
            if RCMEutil.cellContainsString(need, 'dynamics') ...
                    || RCMEutil.cellContainsString(need, 'equilNoGround')
                isGood = true;
            else
                isGood = false;
            end
        end
        
    end
    
end

