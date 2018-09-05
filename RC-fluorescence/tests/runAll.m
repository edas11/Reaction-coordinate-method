clear all;
clc;
addpath('../fluorCore');
addpath('../');

test = testDataLoader();
run(test)

test = testRCfluorescence();
run(test)
