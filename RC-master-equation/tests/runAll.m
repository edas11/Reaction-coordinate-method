clear all;
clc;
addpath('../RCcore');
addpath('../funcScriptData');
addpath('../');

results = [];

test = testCalculateData();
results = [results run(test)];

test = testDataCode();
results = [results run(test)];

test = testRCMEh0();
results = [results run(test)];

test = testRCMEmath();
results = [results run(test)];

test = testRCMEpropogator();
results = [results run(test)];

test = testRCMEutil();
results = [results run(test)];

test = testRCMEximXid();
results = [results run(test)];

test = testRCMEbundle();
results = [results run(test)];

test = testInitialCond();
results = [results run(test)];

test = testDensitySystem();
results = [results run(test)];

test = testDensityFull();
results = [results run(test)];

results
