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


function [Parameter_confidence_intervals] = Calculate_parameter_confidence_GEARS(Param_values, Data, Simulate, Int_opts, Is_regularised, alpha, P_ref)
% A fucntion that calculates parameter confidence using the FIM.
% Param_values   - The parameter values at which parameter confidence should be calculated.
% Data           - The set of data for which the NRMSE will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% Simulate       - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Int_opts       - The intergration options for AMICI (structure) [Optional]
% Is_regularised - If this set of parameter values are from the regularised estimation (logical) [Optional]
% alpha          - The alpha regularisation parameter. Only required if Is_regularised (scaler) [Optional]
% P_ref          - The P_ref regularisation parameter. Only required if Is_regularised (vector) [Optional]


%% Check inputs

Ins = who;
        
    if sum(ismember({'Param_values', 'Data', 'Simulate'}, Ins)) ~= 3
        
    error('The first 3 inputs are required to use this function')    
        
    end       
        
    if ~ismember('Int_opts', Ins)
     
    Int_opts = [];    
        
    end

    if ~ismember('Is_regularised', Ins)
     
    Is_regularised = 0;    
        
    end
    
    
%% Reshape parameter vector

Param_values = reshape(Param_values, [length(Param_values), 1]);    
    
    
%% Add sensitivity setting for AMICI

Int_opts.sensi = 1;

Int_opts.sensi_meth = 'forward';


%% Perform simulation for each experiment and calculate residual

Parameter_confidence_intervals = inf*ones(length(Param_values), 1);

    if Is_regularised

    [~, ~, Jres] = Jacobian_of_residuals_regularised_GEARS(Param_values, Data, Simulate, Int_opts, alpha, P_ref);

    else
    
    [~, ~, Jres] = Jacobian_of_residuals_GEARS(Param_values, Data, Simulate, Int_opts);

    end

Exp_names = fieldnames(Data); 
 
Scaling_tester = [];

    for i = 1:length(Exp_names)

    Scaling_tester = [Scaling_tester; Data.(Exp_names{i}).Scaling_vector];
    
    end

    if ~isequal(Scaling_tester, ones(size(Scaling_tester)))

    [Cost, ~, Residual] = Cost_function_GEARS(Param_values, Data, Simulate, Int_opts, 0);

    Scaling = Cost/(length(Residual) - 1);

    else

    Scaling = 1;

    end


%% FIM

FIM = Jres'*Jres;

var_cov_matrix = zeros(size(FIM));

zero_cols = ~any(FIM);

    if sum(zero_cols) == length(Param_values)

    warning('The FIM is singular')

    return

    end
    
red_FIM = FIM(~zero_cols, ~zero_cols); 

red_FIM = 0.5*(red_FIM + red_FIM'); % Regularise the FIM

r_cond_FIM = rcond(red_FIM);

    if r_cond_FIM <= 1e-50 || isnan(r_cond_FIM)

    warning('The Fisher Information Matrix is singular')
 
    return

    end


%% Covariance matrix

var_cov_matrix(~zero_cols, ~zero_cols) = Scaling*(eye(length(red_FIM))/red_FIM); % = Scaling*Inv(red_FIM)  


%% Parameter confidence

    if sum(zero_cols) ~= length(zero_cols)
    
    Parameter_confidence_intervals(~zero_cols) = 1.96*sqrt(diag(var_cov_matrix(~zero_cols, ~zero_cols)));

    end

end

