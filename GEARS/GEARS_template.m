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


%% File description
% A template that once filled in will run the GEARS procedure. 
% Details for each input can be found in the documentation.
% All fields are required unless stated otherwise.


%% Model infomation

Model_infomation.Model_name = ''; % The name of the model. (string)

Model_infomation.Diff_eqns = {''}; % The model equations. (cell array)
 
Model_infomation.Var_names = {''}; % The variable names. (cell array)

Model_infomation.Fitting_params = {''}; % The parameters that should be estimated. (cell array)

Model_infomation.Fixed_params = {}; % The parameters that have known values. (cell array) [Optional] 

Model_infomation.Param_bounds_lower = []; % The lower bounds for the parameter estimation. (vector)

Model_infomation.Param_bounds_upper = []; % The upper bounds for the parameter estimation. (vector)


%% Fitting data
% For multiple fitting experiments copy and paste this section changing the Fitting_exp_name.

Fitting_exp_name = 'Fitting_experiment'; % The name of the fitting experiment. (string)

Fitting_data.(char(Fitting_exp_name)).Observables = {''}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

    if isequal({'Custom_function'}, Fitting_data.(char(Fitting_exp_name)).Observables) % Only for observation functions.
        
    Fitting_data.(char(Fitting_exp_name)).Custom_function = {''}; % A custom observation function. (cell array) [Optional] 
        
    end

Fitting_data.(char(Fitting_exp_name)).Time_points = []; % The time points for the experiment. (vector)

Fitting_data.(char(Fitting_exp_name)).Measurements = []; % The measurements, should be of the size # time points vs # observables. (matrix)

Fitting_data.(char(Fitting_exp_name)).Standard_deviation = []; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

    if ~isempty(Model_infomation.Fixed_params) % If there are any Fixed parameters
        
    Fitting_data.(char(Fitting_exp_name)).Fixed_params = []; % The values of the fixed parameters. (vector) [Optional]
        
    end

Fitting_data.(char(Fitting_exp_name)).Initial_conditions = []; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Cross-validation data
% For multiple validation experiments copy and paste this section changing the Validation_exp_name.

Validation_exp_name = 'Validation_experiment'; % The name of the fitting experiment. (string)

Validation_data.(char(Validation_exp_name)).Observables = {''}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)

    if isequal({'Custom_function'}, Validation_data.(char(Validation_exp_name)).Observables) % Only for observation functions.
        
    Validation_data.(char(Validation_exp_name)).Custom_function = {''}; % A custom observation function. (cell array) [Optional] 
        
    end

Validation_data.(char(Validation_exp_name)).Time_points = []; % The time points for the experiment. (vector)

Validation_data.(char(Validation_exp_name)).Measurements = []; % The measurements, should be of the size # time points vs # observables. (matrix)

Validation_data.(char(Validation_exp_name)).Standard_deviation = []; % The standard deviation for the measurements, should be same size as the measurements. (matrix)

    if ~isempty(Model_infomation.Fixed_params) % If there are any Fixed parameters
        
    Validation_data.(char(Validation_exp_name)).Fixed_params = []; % The values of the fixed parameters. (vector) [Optional]   
        
    end

Validation_data.(char(Validation_exp_name)).Initial_conditions = []; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)


%% Results folder

Results_folder = ''; % The results folder. (string)


%% Run GEARS 

[Results, Processed_fitting_data, Processed_validation_data] = Run_GEARS(Results_folder, Model_infomation, Fitting_data, Validation_data); % Run GEARS

