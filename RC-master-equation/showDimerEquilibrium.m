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
params.operatoriai = { [1 0;0 0] [0 0;0 1] }; % System - RCs (bath) interraction operators
%
%   ** Reaction coordinate (RC) additions **
%
params.alfa = 11; % Strength of new spectral density (dimensionless units)
params.dimRC = 2; % RC dimension (dimensionless units)

%
%   Calculates
%
bundle = RCMEcoreBundle(params);
operatorBuilder = RCMEoperators(bundle);
equil = operatorBuilder.getH0BoltzmannEquil(false);
dimToTrace = size(equil,2)/params.dim;
equil = RCMEmath.traceOutEnd( params.dim, dimToTrace, equil )

       