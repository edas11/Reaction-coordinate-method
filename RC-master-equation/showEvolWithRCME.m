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
currentMode = availableModes{3};

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
[densityOperator] = this.dynamicsInstance.solve(0:1:500, initConditionParameter);
toc

% Figures
densitySystem = densityOperator.getDensitySystem();
t = densityOperator.getT();
for i=1:2
    for j=1:2
        subplot(2,2,(i-1)*2+j)
        plot(t, squeeze(real(densitySystem(:,i,j))),'b');
        grid on;
        hold on;
        plot(t, squeeze(imag(densitySystem(:,i,j))), 'r');
    end
end
legend('real','imag');
