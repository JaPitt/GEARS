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


function [R2_all_exp, R2_experiment_wise] = Calculate_R2_GEARS(Param_values, Data, Simulate, Int_opts)
% A function that calculates the R2 goodness of fit for a given parameter vector and data set.
% Param_values - A vector of the parameter values at which the NRMSE will be calculated. (vector)
% Data         - The set of data for which the NRMSE will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% Simulate     - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Int_opts     - The intergration options for AMICI (structure) Optional


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

R2_experiment_wise = zeros(length(Exp_names), 1);

Total_data = [];

    for i = 1:length(Exp_names)
       
    Simulation = Simulate(Param_values, Data.(char(Exp_names(i))), Int_opts, Data.(char(Exp_names(i))).Time_points);    

%% R2 for each experiment
    
    Residual_this_exp = Data.(char(Exp_names(i))).Residual_function(Simulation, Data.(char(Exp_names(i))), Param_values);
    
    D = Data.(char(Exp_names(i))).Measurements(Data.(char(Exp_names(i))).Finite_data);    
    
    SS_res = sum(Residual_this_exp.^2);
    
    SS_tot = sum(((D) - mean((D))).^2);
         
    R2_experiment_wise(i) = 1 - SS_res/SS_tot;
        
    Residual = [Residual; Residual_this_exp];
    
    Total_data = [Total_data; D];

    end
    
R2_all_exp = 1 - (sum(Residual.^2)/(sum(Total_data.^2)))/sum((abs(D) - mean(abs(D))).^2);

end
    
    