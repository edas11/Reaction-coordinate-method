clear all;
addpath('fluorCore');

% Angle between dipole moments
fi = 360;
miu{2,1} = [1 0 0];
miu{3,1} = [cos(fi*pi/180) sin(fi*pi/180) 0];

% RC data folder
RCfolder = '../data/RCM-fluor-example/'

% HEOM data folder
HEOMfolder = '..\data\HEOM-fluorescence\';

% System parameters
params.J = 100; % Resonance coupling (fs^-1)
params.T = 300; % Temperature (K)
params.gamma = 1/100; % Relaxation rate (fs^-1)
params.lambda = 100; % Reorganization energy (fs^-1)
params.delta_e = 100; % Energy difference between sites (fs^-1)

% Reaction coordinate (RC) additions
paramsRC.alfa = 8; % Strength of new spectral density (dimensionless units)
paramsRC.dimRC = 8; % RC dimension (dimensionless units)
% Number if zeris added to density operator end.
paramsRC.numOfRCZeros = 50000;
% RC fluoresence spectra shift
paramsRC.energyShift = 15000;

% Calculates
[RCFLinstance] = RCFLmain(params, fi);
[ freqRC, fluorescenceRC ] = RCFLinstance.getRCfluorescence( paramsRC, RCfolder );

%
% Figures
%

% Energy eigenvalues and shifted
Etikrines = paramsRC.energyShift + eig([params.delta_e params.J; params.J 0]);
Epaslinktos = Etikrines - 2*params.lambda;

% Plots eigenvalues, shifted and results
fig = figure(1);
a = plot( freqRC, real(fluorescenceRC)./max(real(fluorescenceRC)), 'LineWidth', 2);
hold on;
b = plot( [Etikrines(1) Etikrines(1)], [-0.1 1], '-r');
plot( [Etikrines(2) Etikrines(2)], [-0.1 1], '-r');
c = plot( [Epaslinktos(1) Epaslinktos(1)], [-0.1 1], '-m');
plot( [Epaslinktos(2) Epaslinktos(2)], [-0.1 1], '-m');
grid on;

set(fig, 'Units', 'centimeters');
set(fig, 'Position', [10 4 13 9]);
axis([12000 18000 0 1.1]);
xlabel('\omega, fs^{-1}');
ylabel('Normuota fluorescencija');
legend([a], 'RC fluorescencija' );