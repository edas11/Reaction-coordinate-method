# Reaction coordinate method

## About
This is an implementation of the reaction coordinate method which is used for simulations of open quantum systems. This program was used to calculate results presented in an article "Simulations of absorption and fluorescence lineshapes using the reaction coordinate method"(see https://www.sciencedirect.com/science/article/pii/S0301010418304099).

## How to use
The program is designed to be launched through scripts which call the main code. There are many scripts throughout the project directories. Works with MATLAB R2013a.

## Directories
**RC-master-equation** - holds the main program that allows calculations of open quantum system dynamics using the RC master equation.

**RC-master-equation/RCcore** - the main program. Calculates dynamics with some set of parameters.

**RC-master-equation/funcScriptData** - basically script that calls the main RCM program with many parameter sets and saves results to **data** directory.

**data** - holds data files needed for RCM programs. Data file names have format (x)y where x is some parameter's name and y is value.

**data/test**, **data/test2** - needed for tests.

**RC-fluorescence** - fluorescence using RCM. It expects to finds its input data in **data** directory.

**RC-absorption** - absorption using RCM. It expects to finds its input data in **data** directory.
