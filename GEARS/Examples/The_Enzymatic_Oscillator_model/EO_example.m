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

Model_information.Model_name = 'EO'; % The name of the model. (string)

Model_information.Diff_eqns = {'dalpha1 = v_Km1_r1 - (alpha1*sigma1_r2*(alpha1 + 1)*(beta1 + 1)^2)/(L1_r2*10^6 + (alpha1 + 1)^2*(beta1 + 1)^2)', ...
'dbeta1 = (50*alpha1*sigma1_r2*(alpha1 + 1)*(beta1 + 1)^2)/(L1_r2*10^6 + (alpha1 + 1)^2*(beta1 + 1)^2) - (sigma2_r3*(gamma1 + 1)^2*(d_r3*beta1/100 + 1)*beta1)/(L2_r3 + (gamma1 + 1)^2*(d_r3*beta1/100 + 1)^2)', ...
'dgamma1 = (sigma2_r3*(gamma1 + 1)^2*(d_r3*beta1/100 + 1)*beta1)/(50*(L2_r3 + (gamma1 + 1)^2*(d_r3*beta1/100 + 1)^2)) - ks_r4*gamma1'}; % The model equations. (cell array)
 
Model_information.Var_names = {'alpha1', 'beta1', 'gamma1'}; % The variable names. (cell array)
% alpha, beta and gamma have a 1 as the names are matlab functions.

Model_information.Fitting_params = {'v_Km1_r1' 'L1_r2', 'sigma1_r2', 'L2_r3', 'd_r3', 'sigma2_r3', 'ks_r4'}; % The parameters that should be estimated. (cell array)

% Model_information.Fixed_params = []; % The parameters that have known values. (cell array) [Optional] 
% All parameters are fit here.

Model_information.Param_bounds_lower = 10^-3*ones(1, length(Model_information.Fitting_params)); % The lower bounds for the parameter estimation. (vector)

Model_information.Param_bounds_upper = 10^3*ones(1, length(Model_information.Fitting_params)); % The upper bounds for the parameter estimation. (vector)


%% Fitting data

% Exp1
Fitting_exp_name = 'Exp1'; % The name of the fitting experiment. (string)

Fitting_data.(char(Fitting_exp_name)).Observables = {'alpha1', 'beta1'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Fitting_data.(char(Fitting_exp_name)).Time_points = [500 1000 1600 2200 2550 3400 3900]; % The time point for the experiment. (vector)

Fitting_data.(char(Fitting_exp_name)).Measurements = [   2.4697e+02   4.1249e+00
   2.1911e+02   8.0437e+01
   4.0314e+02   4.3513e+00
   3.5842e+02   5.5829e+00
   1.7770e+02   2.9944e+00
   2.0113e+02   2.7953e+00
   1.8040e+02   4.7403e+01]; % The measurements, should be of the size # time points vs # observables. (matrix)

Fitting_data.(char(Fitting_exp_name)).Standard_deviation = [   2.3421e+01   5.0000e+00
   2.6276e+01   8.5907e+00
   3.8509e+01   5.0003e+00
   3.4032e+01   5.0001e+00
   1.9801e+01   5.0000e+00
   2.5071e+01   5.0000e+00
   1.8031e+01   6.8743e+00]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Fitting_data.(char(Fitting_exp_name)).Initial_conditions = [2.9200e+01 1.8880e+02 3.3670e-01]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Cross-validation data

% Exp1
Validation_exp_name = 'Exp1'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {'alpha1', 'beta1'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Validation_data.(char(Validation_exp_name)).Time_points = [500 1000 1600 2200 2550 3400 3900]; % The time points for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = [247.1967    5.6488
  173.3258    4.2068
  232.9399   72.5647
  407.7735    4.0536
  196.6142    3.3016
  295.8747    1.0921
  183.8546    1.0826]; % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = [27.9752    5.0000
   19.6230    5.0000
   25.5517    8.5005
   38.5925    5.0003
   24.3298    5.0000
   29.6348    5.0000
   21.2656    5.0000]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Validation_data.(char(Validation_exp_name)).Initial_conditions = [75.6878 60.1984 85.7170]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)

% Exp2
Validation_exp_name = 'Exp2'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {'alpha1', 'beta1'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Validation_data.(char(Validation_exp_name)).Time_points = [500 1000 1600 2200 2550 3400 3900]; % The time point for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = [   2.9532e+02   6.5773e+00
   1.9193e+02   4.9232e+00
   1.6829e+02   3.3548e+01
   3.2801e+02   1.1785e+00
   2.6610e+02   3.7272e+00
   2.9889e+02   6.1728e+00
   2.6689e+02   5.2554e+00]; % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = [   2.9933e+01   5.0000e+00
   2.1559e+01   5.0000e+00
   1.7468e+01   6.2071e+00
   4.0532e+01   5.0007e+00
   2.6279e+01   5.0000e+00
   3.1590e+01   5.0001e+00
   2.3202e+01   5.0000e+00]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Validation_data.(char(Validation_exp_name)).Initial_conditions = [9.5754e+01   8.9283e+01   3.5651e+01]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Results folder

Results_folder = 'EO_example_results'; % The results folder. (string)


%% Run GEARS 

[Results, Processed_fitting_data, Processed_validation_data] = Run_GEARS(Results_folder, Model_information, Fitting_data, Validation_data, 'EO_example_options'); % Run GEARS

