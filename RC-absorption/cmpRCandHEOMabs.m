clear all;
addpath('absCore');

% Angle between dipole moments
fi = 45;
miu{2,1} = [1 0 0];
miu{3,1} = [cos(fi*pi/180) sin(fi*pi/180) 0];

% RC data folder
RCfolder = '../data/RCM-absorp-example/';

% HEOM data folder
HEOMfolder = '..\data\HEOM-absorption\';

% System parameters
params.J = 100; % Resonance coupling (fs^-1)
params.T = 300; % Temperature (K)
params.gamma = 1/100; % Relaxation rate (fs^-1)
params.lambda = 100; % Reorganization energy (fs^-1)
params.delta_e = 100; % Energy difference between sites (fs^-1)

% Reaction coordinate (RC) additions
paramsRC.alfa = 6.28; % Strength of new spectral density (dimensionless units)
paramsRC.dimRC = 8; % RC dimension (dimensionless units)
% Number if zeros added to density operator end.
paramsRC.numOfRCZeros = 50000;
% RC fluoresence spectra shift
paramsRC.energyShift = 15000;

% Calculates
[RCABinstance] = RCABmain(params, fi);
[ freqHEOM, absorptionHEOM, freqRC, absorptionRC ] = RCABinstance.getBothAbsorption(paramsRC, RCfolder, HEOMfolder);
absorptionHEOM = real(absorptionHEOM);
%
% Figures
%

Etikrines = paramsRC.energyShift + eig([params.delta_e params.J; params.J 0]);

fig = figure(1);
a = plot( freqRC, real(absorptionRC)./max(real(absorptionRC)), 'LineWidth', 2);
hold on;
b = plot( [Etikrines(1) Etikrines(1)], [-0.1 1], '-r');
plot( [Etikrines(2) Etikrines(2)], [-0.1 1], '-r');
grid on;

% Plots
maxFreq = 4000;
color = [1 0.8 0.2];
[ freqHEOMreduced, absorptionHEOMreduced ] =...
    eqPointDistance( freqHEOM./maxFreq, absorptionHEOM./max(absorptionHEOM), 0.04 );
d = plot( freqHEOMreduced.*maxFreq, absorptionHEOMreduced, 'ko', 'MarkerSize', 3,...
    'MarkerFaceColor', color);
grid on;

set(fig, 'Units', 'centimeters');
set(fig, 'Position', [10 4 13 9]);
axis([14000 17500 0 1.1]);
xlabel('\omega, cm^{-1}');
ylabel('Normuota sugertis');
legend([a d], 'RC sugertis', 'HEOM sugertis' );