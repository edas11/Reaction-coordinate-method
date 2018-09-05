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
params.alfa = 2*pi; % Strength of new spectral density (dimensionless units)
params.dimRC = 4; % RC dimension (dimensionless units)

% Calculation strategy
% liouville - calculations in liouville space (not efficient)
% hilbert - calculations in hilbert space
availableModes = {'liouville', 'hilbert'};
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
end
[densityOperator] = this.dynamicsInstance.solve(0:1:1200, initConditionParameter);
t = densityOperator.getT();
toc

% Max time(according to RC)
t = t(end);

% Coordinate operators (dimer)
bundle = RCMEcoreBundle(params);
q = bundle.coordRC();
% Average coordinates
for rcId = 1:length(params.operatoriai)
    densRC = densityOperator.getDensityRC(rcId);
    for laikas=1:t
        densityRCcurrent = squeeze(densRC(laikas,:,:));
        qAvg{rcId}(laikas) = trace(q*densityRCcurrent);
    end
end


% Figures
fig = figure(5);
%Nubreziama
plot(0:t-1, real(qAvg{1}), 'LineWidth', 1);
hold on;
plot(0:t-1, real(qAvg{2}), '-r', 'LineWidth', 1);
ylabel( strcat('<q>') )
xlabel( 'Laikas, fs' );
grid on;

legend('RC saveikaujanti su [1 0;0 0]', 'RC saveikaujanti su [0 0;0 1]');

% RC dim and alfa
text(0, 0, strcat('RC dim=', int2str(params.dimRC), ...
    ', alfa=', num2str(params.alfa)), 'BackgroundColor', 'w', 'EdgeColor', 'k' );
set(fig, 'Units', 'centimeters');
set(fig, 'Position', [5 5 15 10]);
