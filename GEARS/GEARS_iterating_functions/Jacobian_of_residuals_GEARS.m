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


function [Jobj, Jg, Jres] = Jacobian_of_residuals_GEARS(Param_values, Data, Simulate, Int_opts, ~)
% The jacobian function used for the initial parameter estimation in GEARS. This function will be iterated over.
% Param_values - The parameter vector in log scale. (vector)
% Data         - The set of data for which the cost will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% Simulate     - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Int_opts     - The intergration options for AMICI (structure) Optional


%% Format parameter vector

Param_values = (reshape(Param_values, [length(Param_values), 1])); % Make sure the parameter vector is a column vector and convert from log form.


%% Find experiment names

Exp_names = fieldnames(Data); % The experiment names.


%% Add sensitivity setting for AMICI

Int_opts.sensi = 1; % Turn on the sensitivity calculations.

Int_opts.sensi_meth = 'forward'; % Turn on the sensitivity calculations.


%% Perform simulation for each experiment and calculate residual

Residual = [];

Jres = [];

    for i = 1:length(Exp_names)
       
    Simulation = Simulate(Param_values, Data.(Exp_names{i}), Int_opts, Data.(Exp_names{i}).Time_points); % Simulate the model.    

%% Residual for each experiment

    Residual = [Residual; Data.(Exp_names{i}).Scaling_vector.*Data.(Exp_names{i}).Residual_function(Simulation, Data.(Exp_names{i}), Param_values)]; % Calculate the residuals for this experiment.
    
    
%% Sensitivity for each experiment

        if Data.(char(Exp_names(i))).Param_initials % Include the sensitivity from any initial conditions that are dependant on parameters. 

        Sens_initials = Data.(Exp_names{i}).Initial_jac_format_function(Data.(Exp_names{i}).Initial_jac_scaling_function(Param_values), Simulation.sx(:, :, length(Param_values) + 1 - sum(Data.(Exp_names{i}).Initial_params_only):end), Data.(Exp_names{i}));
        
        Simulation.sx(:,:,1:length(Param_values)) = Simulation.sx(:, :, 1:length(Param_values)) + Sens_initials; % Add the sensitivity from initial conditions.            
        
        end

    Sens = Data.(Exp_names{i}).Jacobian_function(Simulation, Data.(Exp_names{i}), Param_values); % Format sensitivities
    
    Sens = repmat(Data.(Exp_names{i}).Scaling_vector(Data.(Exp_names{i}).Finite_data), [1, length(Param_values)]).*Sens./repmat(Data.(Exp_names{i}).Standard_deviation(Data.(Exp_names{i}).Finite_data), [1, length(Param_values)]); % Include the standard deviation and scaling vector.
    
    Jres = [Jres; Sens];
    
    end
   
   
%% Combine residuals and calculate jacobian of the objective function

Jobj = 2*Residual'*Jres;

Jg = 0;

end

