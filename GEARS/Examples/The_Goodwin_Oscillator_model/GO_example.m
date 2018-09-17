% Copyright 2018 Jake Alan Pitt and Julio R. Banga
 
% This file is part of GEARS.
 
% GEARS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
 
% GEARS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
 
% You should have received a copy of the GNU General Public License
% along with GEARS.  If not, see <http://www.gnu.org/licenses/>.
 
% File author: Jake Alan Pitt (jp00191.su@gmail.com)


%% Model information

Model_information.Model_name = 'GO'; % The name of the model. (string)

Model_information.Diff_eqns = {'dx1 = (k1*Ki^n)/(Ki^n + x3^n) - k2*x1', 'dx2 = k3*x1 - k4*x2', 'dx3 = k5*x2 - k6*x3'}; % The model equations. (cell array)
 
Model_information.Var_names = {'x1','x2', 'x3'}; % The variable names. (cell array)

Model_information.Fitting_params = {'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'Ki', 'n'}; % The parameters that should be estimated. (cell array)

% Model_information.Fixed_params = []; % The parameters that have known values. (cell array) [Optional] 
% All parameters are fit here.

Model_information.Param_bounds_lower = [10^-3*ones(1, 7) 1]; % The lower bounds for the parameter estimation. (vector)

Model_information.Param_bounds_upper = [10^3*ones(1, 7) 12]; % The upper bounds for the parameter estimation. (vector)


%% Fitting data

Fitting_exp_name = 'Fitting_experiment'; % The name of the fitting experiment. (string)

Fitting_data.(char(Fitting_exp_name)).Observables = {'x1', 'x3'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Fitting_data.(char(Fitting_exp_name)).Time_points = [1.0000
   27.5600
   54.1100
   80.6700
  107.2000
  133.8000
  160.3000
  186.9000
  213.4000
  240.0000]; % The time point for the experiment. (vector)

Fitting_data.(char(Fitting_exp_name)).Measurements = [0.0821    2.5180
    0.0041    2.9610
    0.0351    2.4740
    0.0062    1.7450
    0.0125    2.5940
    0.0463    1.7840
    0.0061    2.2890
    0.0103    2.8900
    0.0489    1.4620
    0.0046    1.2890]; % The measurements, should be of the size # time points vs # observables. (matrix)

Fitting_data.(char(Fitting_exp_name)).Standard_deviation = [0.0095    0.2682
    0.0031    0.3100
    0.0048    0.2809
    0.0031    0.2106
    0.0032    0.3031
    0.0058    0.2174
    0.0031    0.2248
    0.0033    0.2937
    0.0059    0.2054
    0.0031    0.2278]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Fitting_data.(char(Fitting_exp_name)).Initial_conditions = [0.1000    0.2000    2.5000]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Cross-validation data

% exp1
Validation_exp_name = 'exp1'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {'x1', 'x2', 'x3'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Validation_data.(char(Validation_exp_name)).Time_points = [36 72 108 144 180 216 252 288 324 360]; % The time points for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = [0.0064    0.0972    1.7050
    0.0075    0.1045    2.2380
    0.0077    0.1666    2.3530
    0.0089    0.2427    2.4400
    0.0134    0.3088    2.7620
    0.0218    0.3326    2.3340
    0.0351    0.3297    2.4310
    0.0371    0.3379    1.9600
    0.0567    0.2844    1.6850
    0.0359    0.1629    1.5910]; % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = [  0.0032    0.0146    0.1972
    0.0031    0.0161    0.2276
    0.0031    0.0193    0.2561
    0.0031    0.0232    0.2779
    0.0033    0.0275    0.2890
    0.0036    0.0314    0.2855
    0.0042    0.0340    0.2653
    0.0050    0.0332    0.2318
    0.0059    0.0271    0.1989
    0.0047    0.0185    0.1861]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Validation_data.(char(Validation_exp_name)).Initial_conditions = [0.0847 0.0633 3.0830]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)

% exp2
Validation_exp_name = 'exp2'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {'x1', 'x2', 'x3'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Validation_data.(char(Validation_exp_name)).Time_points = [36 72 108 144 180 216 252 288 324 360]; % The time points for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = [0.0030    0.0854    1.6200
    0.0060    0.1711    2.4850
    0.0126    0.2062    2.7360
    0.0114    0.2267    2.6990
    0.0246    0.3090    2.4410
    0.0294    0.3186    2.4360
    0.0618    0.2328    1.8020
    0.0386    0.1703    1.5750
    0.0214    0.1063    1.5570
    0.0081    0.1366    2.2060]; % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = [0.0030    0.0143    0.2199
    0.0031    0.0193    0.2686
    0.0032    0.0250    0.2945
    0.0034    0.0303    0.2977
    0.0039    0.0340    0.2789
    0.0048    0.0343    0.2429
    0.0058    0.0289    0.2045
    0.0051    0.0196    0.1862
    0.0034    0.0153    0.1942
    0.0031    0.0156    0.2149]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Validation_data.(char(Validation_exp_name)).Initial_conditions = [0.0965 0.1346 3.9680]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Results folder

Results_folder = 'GO_example_results'; % The results folder. (string)


%% Run GEARS 

% [Results, Processed_fitting_data, Processed_validation_data] =
 Run_GEARS(Results_folder, Model_information, Fitting_data, Validation_data, 'GO_example_options'); % Run GEARS

