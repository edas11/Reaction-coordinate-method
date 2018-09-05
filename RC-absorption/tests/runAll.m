clear all;
clc;
addpath('../absCore');
addpath('../');

test = testDataLoader();
run(test)

test = testRCabsorption();
run(test)
