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


function [Conclusion, p] = Calculate_chi2_GEARS(Param_values, Data, Simulate, Int_opts)
% A function that tests the chi2 hypothesis for a given parameter vector and data set.
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

    for i = 1:length(Exp_names)
       
    Simulation = Simulate(Param_values, Data.(char(Exp_names(i))), Int_opts, Data.(char(Exp_names(i))).Time_points);    

%% R2 for each experiment
    
    Residual_this_exp = Data.(char(Exp_names(i))).Residual_function(Simulation, Data.(char(Exp_names(i))), Param_values);         
    
    Residual = [Residual; Residual_this_exp];
    
    end
    
DOF = numel(Residual) - length(Param_values); 

chi2 = sum(Residual.^2);
  
p = 1 - gammainc(chi2/2, DOF/2); % Requires stats package

% Significance level for test:
alpha = 0.01;


%% Conclusion

    if p < alpha
        
    Conclusion = 'We reject the chi2 hypothesis';        
        
    else 
        
    Conclusion = 'We do not reject the chi2 hypothesis';  
    
    end
    
