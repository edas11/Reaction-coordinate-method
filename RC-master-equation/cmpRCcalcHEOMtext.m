clear all;
addpath('RCcore');
addpath('funcScriptData');
%
%   ** System parameters **
%
DALIKLIS = 5308.837458877; % Conversion between cm^-1 and fs^-1 (cm^-1/fs^-1)
params.dim = 2;   % System dimension (dimensionless units)
params.J = 100/DALIKLIS; % Resonance coupling (fs^-1)
params.T = 300; % Temperature (K)
params.gamma = 1/100; % Relaxation rate (fs^-1)
params.lambda = 100/DALIKLIS; % Reorganization energy (fs^-1)
params.delta_e = 100/DALIKLIS; % Energy difference between sites (fs^-1)
params.operatoriai = { [1 0;0 0] [0 0;0 1]}; % System - RCs (bath) interraction operators
%
%   ** Reaction coordinate (RC) additions **
%
params.alfa = 5; % Strength of new spectral density (dimensionless units)
params.dimRC = 5; % RC dimension (dimensionless units)

% Calculation strategy
% liouville - calculations in liouville space (not efficient)
% hilbert - calculations in hilbert space
% hilbert-traced - calculations in hilbert space but RCs are traced out.
%                  Memory efficient, but information about RCs is lost.
availableModes = {'liouville', 'hilbert', 'hilbert-traced'};
numberOfBlocks = 4; % Only used by hilbert-traced, increasing it should reduce memory usage.
currentMode = availableModes{2};

% Initial conditions
% direct - parameter is initial system condition
% equilWithMiuShift - parameter is matrix miu and initial cond is
%                     equilibrium density operator * miu
initConditionParameter.mode = 'direct';
initConditionParameter.parameter = [1 0;0 0];

%
%   Calculates
%
tic
if strcmp(currentMode,'liouville')
    this.dynamicsInstance = RCMEdynamicsLiouville(params);
elseif strcmp(currentMode,'hilbert')
    this.dynamicsInstance = RCMEdynamicsHilbert(params);
elseif strcmp(currentMode,'hilbert-traced')
    this.dynamicsInstance = RCMEdynamicsHilbertTraced(params, numberOfBlocks);
end
[densityOperator] = this.dynamicsInstance.solve(0:500, initConditionParameter);
exec = toc
densitySystem = densityOperator.getDensitySystem();
t = densityOperator.getT();

% Takes HEOM dynamics
paramString = createParamString({'(T)' '(g-1)' '(d)' '(l)' '(J)'}, ...
    [params.T 1/params.gamma params.delta_e*DALIKLIS params.lambda*DALIKLIS ...
    params.J*DALIKLIS]);
folderHEOM = strcat('..\data\HEOM-dynamics-(IC)11\');
densityHEOM = importHEOM(folderHEOM, paramString);

% Max time(according to RC)
t = t(end);
    
% Figures
fig = figure(3);
for r=1:2
    for j=1:2
        subplot(2,2,(r-1)*2+j);
        grid on;
        plot(1:t, squeeze(real(densitySystem(1:t,r,j))));
        hold on;
        plot(1:10:t, squeeze(real(densityHEOM(1:10:t,r,j))), 'o', ...
            'MarkerSize', 3);
        plot(1:t, squeeze(imag(densitySystem(1:t,r,j))), 'r');
        plot(1:10:t, squeeze(imag(densityHEOM(1:10:t,r,j))), 'ro', ...
            'MarkerSize', 3);
        ylabel( strcat('\rho_{', int2str(10*r+j), '}') )
        xlabel( 'Laikas, fs' );
    end
end
% RC dim and alfa
text(0, 0, strcat('RC dim=', int2str(params.dimRC), ', Exec=', num2str(exec), ...
    ', gamma=', num2str(params.alfa)), 'BackgroundColor', 'w', 'EdgeColor', 'k' );
legend('real RC', 'real HEOM', 'imag RC', 'imag HEOM');
set(fig, 'Units', 'centimeters');
set(fig, 'Position', [5 5 15 10]);
