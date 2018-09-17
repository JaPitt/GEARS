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

Model_information.Model_name = 'AP'; % The name of the model. (string)

Model_information.Diff_eqns = {'dx1=-(p1+p2)*x1'; 'dx2= p1*x1'; 'dx3= p2*x1-(p3+p4)*x3+p5*x5'; 'dx4= p3*x3'; 'dx5= p4*x3-p5*x5'}; % The model equations. (cell array)
 
Model_information.Var_names = {'x1', 'x2', 'x3', 'x4', 'x5'}; % The variable names. (cell array)

Model_information.Fitting_params = {'p1', 'p2', 'p3', 'p4', 'p5'}; % The parameters that should be estimated. (cell array)
% p1=5.93e-5;  p2=2.96e-5;  p3=2.05e-5;  p4=27.5e-5;  p5=4e-5;

% Model_information.Fixed_params = []; % The parameters that have known values. (cell array) [Optional] 
% All parameters are fit here.

Model_information.Param_bounds_lower = [10^-6*ones(1, length(Model_information.Fitting_params))]; % The lower bounds for the parameter estimation. (vector)

Model_information.Param_bounds_upper = [100*ones(1, length(Model_information.Fitting_params))]; % The upper bounds for the parameter estimation. (vector)


%% Fitting data 

Fitting_exp_name = 'exp1'; % The name of the fitting experiment. (string)

Fitting_data.(char(Fitting_exp_name)).Time_points = [1230 3060 4920 7800 10680 15030 22620 36420]; % The time points for the experiment. (vector)

Fitting_data.(char(Fitting_exp_name)).Measurements = [88.35 7.3   2.3 0.4 1.75
76.4  15.6  4.5 0.7  2.8
65.1  23.1  5.3  1.1 5.8
50.4  32.9  6.0  1.5 9.3
37.5  42.7  6.0  1.9 12.0
25.9  49.1  5.9  2.2 17.0
14.0  57.4  5.1  2.6 21.0
4.5   63.1  3.8  2.9 25.7]; % The measurements, should be of the size # time points vs # observables. (matrix)

Fitting_data.(char(Fitting_exp_name)).Standard_deviation = ones(size(Fitting_data.(char(Fitting_exp_name)).Measurements)); 

Fitting_data.(char(Fitting_exp_name)).Initial_conditions = [100 0 0 0 0]; % The initial conditions for the states. Use a cell string for initial conditions that depend on parameters. (vector or cell array)

Fitting_data.(char(Fitting_exp_name)).Observables = {'x1', 'x2', 'x3', 'x4', 'x5'}; % The observables, should be either a subsection of the variable names or {'Custom_function'} for observed functions. (cell array)


%% Results folder

Results_folder = 'AP_example_results'; % The results folder. (string) 


%% Run GEARS 

[Results, Processed_fitting_data, Processed_validation_data] = Run_GEARS(Results_folder, Model_information, Fitting_data, [], 'AP_example_options'); % Run GEARS

