clear all;
addpath('fluorCore');

% Angle between dipole moments
fi = 45;
miu{2,1} = [1 0 0];
miu{3,1} = [cos(fi*pi/180) sin(fi*pi/180) 0];

% RC data folder
RCfolder = '../data/RCM-fluor-example/';

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

% File name
names = {'(T)' '(g-1)' '(d)' '(l)' '(J)' '(fi)' '(RCalfa)' '(RCdim)'};
parameters = [params.T 1/params.gamma params.delta_e...
                params.lambda params.J fi paramsRC.alfa paramsRC.dimRC];
paramString = RCFLdataLoader.createParamString(names,parameters);
fileName = strcat('RC-', paramString, '.txt');

% Removes some points
indMin = find(freqRC>10000, 1);
indMax = find(freqRC>20000, 1);
freqRC = freqRC(indMin:indMax);
fluorescenceRC = fluorescenceRC(indMin:indMax);

% Saves
id = fopen(fileName, 'w');
fprintf(id, '%s\t%s\n', 'freqRC', 'fluorescenceRC');
for i=1:length(fluorescenceRC)
    fprintf(id, '%.4f\t%.4f\n', freqRC(i), fluorescenceRC(i));
end
fclose(id);
