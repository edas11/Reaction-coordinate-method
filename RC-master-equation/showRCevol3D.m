clear all;
addpath('RCcore');
addpath('funcScriptData');
%
%   ** System parameters **
%
DALIKLIS = 5308.837458877; % Conversion between cm^-1 and fs^-1 (cm^-1/fs^-1)
params.dim = 3;   % System dimension (dimensionless units)
params.J = 100/DALIKLIS; % Resonance coupling (fs^-1)
params.T = 300; % Temperature (K)
params.gamma = 1/100; % Relaxation rate (fs^-1)
params.lambda = 100/DALIKLIS; % Reorganization energy (fs^-1)
params.delta_e = 100/DALIKLIS; % Energy difference between sites (fs^-1)
params.operatoriai = { [0 0 0;0 1 0;0 0 0] [0 0 0;0 0 0;0 0 1] }; % System - RCs (bath) interraction operators
%
%   ** Reaction coordinate (RC) additions **
%
params.alfa = 7; % Strength of new spectral density (dimensionless units)
params.dimRC = 3; % RC dimension (dimensionless units)

% Calculation strategy
% liouville - calculations in liouville space (not efficient)
% hilbert - calculations in hilbert space
availableModes = {'liouville', 'hilbert'};
currentMode = availableModes{1};

% Initial conditions
% direct - parameter is initial system condition
% equilWithMiuShift - parameter is matrix miu and initial cond is
%                     equilibrium density operator * miu
initConditionParameter.mode = 'direct';
initConditionParameter.parameter = [0 0 0;0 1 0;0 0 0];

%
%   Calculates
%
tic
if strcmp(currentMode,'liouville')
    this.dynamicsInstance = RCMEdynamicsLiouville(params);
elseif strcmp(currentMode,'hilbert')
    this.dynamicsInstance = RCMEdynamicsHilbert(params);
end
[densityOperator] = this.dynamicsInstance.solve(0:1:500, initConditionParameter);
densRC = densityOperator.getDensityRC(1);
t = densityOperator.getT();
toc

% Figures
for i=1:size(densRC, 2)
        figure(i);
        plot(t, squeeze(real(densRC(:,i,i))), '-b');
        grid on;
        hold on;
        plot(t, squeeze(imag(densRC(:,i,i))), '-r');
end
legend('real','imag');

