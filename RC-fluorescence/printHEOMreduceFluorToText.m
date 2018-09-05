clear all;
addpath('fluorCore');

% Angle between dipole moments
fi = 45;
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

% Gets fluorescence
[RCFLinstance] = RCFLmain(params, fi);
[freqHEOM, fluorescenceHEOM] = RCFLinstance.getHEOMfluorescence( HEOMfolder );
fluorescenceHEOM = real(fluorescenceHEOM);

% File name
names = {'(T)' '(g-1)' '(d)' '(l)' '(J)' '(fi)'};
parameters = [params.T 1/params.gamma params.delta_e...
                params.lambda params.J fi];
paramString = RCFLdataLoader.createParamString(names,parameters);
fileName = strcat('HEOM-', paramString, '.txt');

% Removes some points
maxFreq = 2000;
[ freqHEOM, fluorescenceHEOM ] =...
    eqPointDistance( freqHEOM./maxFreq, fluorescenceHEOM./max(fluorescenceHEOM), 0.08 );
freqHEOM = freqHEOM.*maxFreq;

% Saves
id = fopen(fileName, 'w');
fprintf(id, '%s\t%s\n', 'freqHEOM', 'fluorescenceHEOM');
for i=1:length(fluorescenceHEOM)
    fprintf(id, '%.4f\t%.4f\n', freqHEOM(i), fluorescenceHEOM(i));
end
fclose(id);
