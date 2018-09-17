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

Model_information.Model_name = 'RP'; % The name of the model. (string)

Model_information.Diff_eqns = {'dp1 = beta1*(m1-p1)', 'dp2 = beta1*(m2-p2)', 'dp3 = beta1*(m3-p3)', ...
'dm1 = alpha0 + alpha1/(1+p3^n) - m1', 'dm2 = alpha0 + alpha1/(1+p1^n) - m2', 'dm3 = alpha0 + alpha1/(1+p2^n) - m3'}; % The model equations. (cell array)
 
Model_information.Var_names = {'p1', 'p2', 'p3', 'm1', 'm2', 'm3'}; % The variable names. (cell array)

Model_information.Fitting_params = {'alpha0', 'alpha1', 'n', 'beta1'}; % The parameters that should be estimated. (cell array)
% alpha1 and beta1 are used as alpha and beta are matlab function names.

% Model_information.Fixed_params = []; % The parameters that have known values. (cell array) [Optional]
% All parameters are fit here.

Model_information.Param_bounds_lower = [10^-3 10^-3 1 10^-3]; % The lower bounds for the parameter estimation. (vector)

Model_information.Param_bounds_upper = [500 500 10 500]; % The upper bounds for the parameter estimation. (vector)


%% Fitting data

% exp1
Fitting_exp_name = 'Experiment_1'; % The name of the fitting experiment. (string)

Fitting_data.(char(Fitting_exp_name)).Observables = {'m3'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Fitting_data.(char(Fitting_exp_name)).Time_points = [10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200]; % The time points for the experiment. (vector)

Fitting_data.(char(Fitting_exp_name)).Measurements = [34.6327
    4.6454
    1.2468
    6.0185
  282.9133
  268.1860
    2.0021
    4.5692
    0.5448
    8.4656
  326.9984
  261.1481
    5.3017
    0.8844
    3.9266
    6.8485
  303.8762
  368.5677
    7.6193
   10.2852]; % The measurements, should be of the size # time points vs # observables. (matrix)

Fitting_data.(char(Fitting_exp_name)).Standard_deviation = [5.8942
    5.0000
    5.0000
    5.0000
   27.9127
   30.2214
    5.0152
    5.0000
    5.0000
    5.0000
   29.8212
   30.2215
    5.0005
    5.0000
    5.0000
    5.0000
   30.1528
   30.2215
    5.0000
    5.0000]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Fitting_data.(char(Fitting_exp_name)).Initial_conditions = [10.0000 0.0100 1.0000  1.0000 0.0100 10.0000]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Validation data

% Val_exp_1
Validation_exp_name = 'Val_exp_1'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {'p3', 'm3'};

Validation_data.(char(Validation_exp_name)).Time_points = [10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200]; % The time points for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = [62.3283    0.4730
    4.6043    1.6908
    9.1423    3.8183
    1.1345    5.7067
  254.4451  348.4076
  303.9485  274.0937
   38.7122    0.5844
    6.3964    0.8646
    3.3875    3.4589
    5.4487    4.5275
  248.8146  318.7132
  289.5650  358.5250
   23.7054    7.0345
    0.9688    1.5784
    8.5539    2.6369
   14.9342   44.1457
  230.6689  286.6180
  221.6045   88.3544
    9.9416    2.5002
    0.7643    7.9792];  % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = [8.3194    5.0105
    5.0117    5.0000
    5.0000    5.0000
    5.0000    5.0000
   23.9004   30.1531
   29.9042   30.2215
    6.3906    5.0000
    5.0041    5.0000
    5.0000    5.0000
    5.0000    5.0009
   26.4802   30.2098
   30.0346   30.2189
    5.5238    5.0000
    5.0015    5.0000
    5.0000    5.0000
    5.0381    6.9275
   28.0140   30.2195
   25.1235    8.8145
    5.1883    5.0000
    5.0005    5.0000]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Validation_data.(char(Validation_exp_name)).Initial_conditions = [5.7196 0.9749 5.3793 3.5145 3.0534 12.2069]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


% Val_exp_2
Validation_exp_name = 'Val_exp_2'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {'p3', 'm3'};  % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

Validation_data.(char(Validation_exp_name)).Time_points = [10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200]; % The time points for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = [2.9125    4.9417
    1.1111    4.5559
  198.3451  271.4090
  292.8764  265.0093
   31.3105    6.4622
    2.6392    6.5832
    3.1300    4.9686
   13.0680   11.8916
  259.2459  319.7330
  302.9411  119.7499
   16.8771    1.2202
    0.6217   10.5540
    3.3700    4.3391
   35.0610  187.9838
  302.9514  300.4717
  180.2383   17.2577
    3.5605    3.1002
    5.0200    3.5959
    1.3242    1.4523
  135.9557  251.5079]; % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = [5.0029    5.0002
    5.0000    5.0000
   25.5685   30.1971
   29.9887   30.2215
    5.7021    5.0000
    5.0020    5.0000
    5.0000    5.0000
    5.0041    5.2874
   27.6441   30.2181
   27.6289   13.1697
    5.2550    5.0000
    5.0007    5.0000
    5.0000    5.0000
    6.4876   18.9951
   28.7019   30.2209
   18.8611    5.4195
    5.0903    5.0000
    5.0003    5.0000
    5.0000    5.0000
   13.4639   28.0752]; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

Validation_data.(char(Validation_exp_name)).Initial_conditions = [11.7447 8.4997 12.1700 8.6521 14.1605 13.0719]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Results folder

Results_folder = 'RP_example_results'; % The results folder. (string)


%% Run GEARS 

[Results, Processed_fitting_data, Processed_validation_data] = Run_GEARS(Results_folder, Model_information, Fitting_data, Validation_data, 'RP_example_options'); % Run GEARS

