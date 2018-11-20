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


function [NRMSE, NRMSE_experiment_wise] = Calculate_NRMSE_GEARS(Param_values, Data, Simulate, Int_opts)
% A function that calculates the normalised root mean square error (NRMSE) for a given parameter vector and data set.
% Param_values - A vector of the parameter values at which the NRMSE will be calculated. (vector)
% Data         - The set of data for which the NRMSE will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% Simulate     - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Int_opts     - The intergration options for AMICI (structure) [Optional]


%% Check inputs 

Ins = who;

    if sum(ismember({'Param_values', 'Data', 'Simulate'}, Ins)) ~= 3
        
    error('The first 3 inputs are required')       

    end

    if ~ismember('Int_opts', Ins)
     
    Int_opts = [];    
        
    end
        
            
%% Reshape parameter vector

Param_values = reshape(Param_values, [length(Param_values), 1]);


%% Find how many experiments are being considered

Exp_names = fieldnames(Data);


%% Calculate normalised residuals for each experiment

Residual = [];

NRMSE_experiment_wise = zeros(length(Exp_names), 1);

    for i = 1:length(Exp_names)
       
    Simulation = Simulate(Param_values, Data.(char(Exp_names(i))), Int_opts, Data.(char(Exp_names(i))).Time_points);    

%% Residual for each experiment

    Norm = repmat(max(Data.(char(Exp_names(i))).Measurements) - min(Data.(char(Exp_names(i))).Measurements), [size(Data.(char(Exp_names(i))).Measurements, 1), 1]);
    
    Norm = Norm(Data.(char(Exp_names(i))).Finite_data);
    
    Residual_this_exp = Data.(char(Exp_names(i))).Residual_function(Simulation, Data.(char(Exp_names(i))), Param_values)./Norm;
    
    Residual_this_exp = Residual_this_exp.*Data.(char(Exp_names(i))).Standard_deviation(Data.(char(Exp_names(i))).Finite_data);

    NRMSE_experiment_wise(i) = sqrt(sum(Residual_this_exp.^2)/numel(Residual_this_exp)); % Do not use the standard deviation in the NRMSE.
    
    Residual = [Residual; Residual_this_exp];

    end

    
%% Calculate NRMSE  
    
NRMSE = sqrt(sum(Residual.^2)/numel(Residual));

    
end

