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


function [Cost, g, Residual, Sampling_output] = Cost_function_GEARS(Param_values, Data, Simulate, Int_opts, Sample)
% The cost function used for the initial parameter estimation in GEARS. This function will be iterated over.
% Param_values       - The parameter vector in log scale. (vector)
% Data               - The set of data for which the cost will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% Simulate           - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Int_opts           - The intergration options for AMICI (structure) Optional
% Sample             - A flag indicating if the cost function should use the sampling section. (logical)


%% Format parameter vector

Param_values = (reshape(Param_values, [length(Param_values), 1])); % Make sure the parameter values are in a column vector and convert from log form.


%% Find experiment names

Exp_names = fieldnames(Data); % The experiment names.


%% Perform simulation for each experiment and calculate residual

Residual = [];

    for i = 1:length(Exp_names)

    Simulation = Simulate(Param_values, Data.(Exp_names{i}), Int_opts, Data.(Exp_names{i}).Time_points); % Simulate the model.

    
%% Residual for each experiment

    Residual = [Residual; Data.(Exp_names{i}).Scaling_vector.*Data.(Exp_names{i}).Residual_function(Simulation, Data.(Exp_names{i}), Param_values)]; % Calulate residuals for this experiment.

    end
      
%% Organise outputs
    
Cost = Residual'*Residual; 

g = 0; % We do not handle problems with constraints


%% Save sample points 

    if Sample % If we are sampling save the infomation for this iteration.

    persistent Num_sample_points Param_sample_points Param_sample_costs % Static variables

        if isempty(Num_sample_points) % Initialise the variables

        Num_sample_points = 0;

        Param_sample_costs = inf*ones(10^5, 1);

        Param_sample_points = zeros(10^5, length(Param_values));
          
        end

    Num_sample_points = Num_sample_points + 1;

    Param_sample_costs(Num_sample_points) = Cost;

    Param_sample_points(Num_sample_points, :) = Param_values';

    Sampling_output.Num_sample_points = Num_sample_points; 

    Sampling_output.Param_sample_costs = Param_sample_costs;

    Sampling_output.Param_sample_points = Param_sample_points;
    
    else

    Sampling_output = []; % To define the output if no sampling is performed
    
    end
   
end

