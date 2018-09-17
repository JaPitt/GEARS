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


function [Cost, g, Residual] = Cost_function_regularised_GEARS(Param_values, Data, Simulate, Int_opts, alpha, P_ref)
% The cost function used for the regularised parameter estimation in GEARS. This function will be iterated over.
% Param_values - The parameter vector in log scale. (vector)
% Data         - The set of data for which the cost will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% Simulate     - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Int_opts     - The intergration options for AMICI (structure) Optional
% alpha        - A regularisation parameter (scaler)
% P_ref        - A regularisation parameter. This should not be in log scale. (vector) 


%% Format parameter vector

Param_values = (reshape(Param_values, [length(Param_values), 1])); % Make sure the Param_values are a column vector and convert from log form.


%% Find experiment names

Exp_names = fieldnames(Data); % The experiment names.


%% Perform simulation for each experiment and calculate residual

Residual = [];

    for i = 1:length(Exp_names)
       
    Simulation = Simulate(Param_values, Data.(Exp_names{i}), Int_opts, Data.(Exp_names{i}).Time_points); % Simulate the model.
    

%% Residual for each experiment

    Residual = [Residual; Data.(Exp_names{i}).Scaling_vector.*Data.(Exp_names{i}).Residual_function(Simulation, Data.(Exp_names{i}), Param_values)]; % Calculate the residuals for this experiment.

    end
      
%% Organise outputs

Residual = [Residual; sqrt(alpha)*diag(1./(P_ref))*(Param_values - P_ref)]; % Add the regularisation term
    
Cost = Residual'*Residual;

g = 0; % We do not handle problems with constraints

end

