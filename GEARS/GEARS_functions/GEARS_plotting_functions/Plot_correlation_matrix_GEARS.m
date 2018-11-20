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


function [Figure_handle, Corr_matrix] = Plot_correlation_matrix_GEARS(Param_values, Data, Simulate, Fitting_params, Save_loc, Int_opts, Is_regularised, alpha, P_ref)
% A function that plots the parameter correlation matrix for a given parameter value.
% Param_values   - A vector of the parameter values at which the NRMSE will be calculated. (vector)
% Data           - The set of data for which the NRMSE will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% Simulate       - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Fitting_params - The names of the parameters that have been fit. (cell array)
% Save_loc       - The location in which the figures should be saved. (string)
% Int_opts       - The intergration options for AMICI (structure) Optional
% Is_regularised - If this set of parameter values are from the regularised estimation (logical) [Optional]
% alpha          - The alpha regularisation parameter. Only required if Is_regularised (scaler) [Optional]
% P_ref          - The P_ref regularisation parameter. Only required if Is_regularised (vector) [Optional]


%% Check inputs

Ins = who;

    if sum(ismember({'Param_values', 'Data', 'Simulate', 'Fitting_params','Save_loc'}, Ins)) ~= 5
        
    error('The first 5 inputs are required to use this function')   
        
    end   

    if ~ischar(Save_loc)
        
    error('The save location should be given as a string')
        
    end
    
    if ~ismember('Int_opts', Ins)
     
    Int_opts = [];    
        
    end    

    if ~ismember('Is_regularised', Ins)
     
    Is_regularised = 0;    
        
    end

    
%% Reshape parameter vector

Param_values = reshape(Param_values, [length(Param_values), 1]);


%% Calculate Jacobian

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

    [Cost, ~, Residual] = Cost_function_GEARS(Param_values, Data, Simulate, Int_opts, 0, 0);

    Scaling = Cost/(length(Residual) - 1);

    else

    Scaling = 1;

    end

    
%% FIM

Correlation_plot = figure('color', 'w');

Correlation_plot.Position = [0 0 1250 750];

FIM = Jres'*Jres;

var_cov_matrix = zeros(size(FIM));

zero_cols = ~any(FIM);

red_FIM = FIM(~zero_cols, ~zero_cols); 

red_FIM = 0.5*(red_FIM + red_FIM'); % Regularise the FIM

r_cond_FIM = rcond(red_FIM);

Corr_matrix = zeros(size(red_FIM));

    if sum(zero_cols) == length(Param_values)

    warning('The Fisher Information Matrix is singular, it is not possible to calculate the correlation matrix. This will be skipped.')
    
    Figure_handle = Correlation_plot;

    return

    end

    if r_cond_FIM <= 1e-50 || isnan(r_cond_FIM)

    warning('The Fisher Information Matrix is singular, it is not possible to calculate the correlation matrix. This will be skipped.')
    
    Figure_handle = Correlation_plot;
 
    return

    end


%% Calculate correlation matrix

var_cov_matrix(~zero_cols, ~zero_cols) = Scaling*(eye(length(red_FIM))/red_FIM); % = Scaling*Inv(red_FIM)

var_cov_matrix_red = var_cov_matrix(~zero_cols, ~zero_cols);
 
    for i = 1:size(var_cov_matrix_red, 1)

        for j = 1:size(var_cov_matrix_red, 1)

        Corr_matrix(i, j) = var_cov_matrix_red(i, j)/(abs((var_cov_matrix_red(i, i)*var_cov_matrix_red(j, j)))^0.5);

        end

    end


%% Plot Correlation matrix

Param_names_red = Fitting_params(~zero_cols);

M = zeros(length(Param_names_red) + 1);

M(length(Param_names_red):-1:1, 1:length(Param_names_red)) = Corr_matrix;

pcolor(real(M))

Colour_map = [3, 17, 81
3, 34, 178
52, 100, 249
140, 182, 255
206, 224, 255
255 255 255
255, 221, 211
226 179 179
226 111 111
204 40 40
122 0 0]/255; % /255 to go from rbg to matlab rgb.

colormap(gca, flipud(Colour_map))

colorbar

caxis([-1 1])

title('Parameter correlations')

set(gca, 'Xtick', (1:length(Param_names_red)) + 0.5,'Ytick', (1:length(Param_names_red)) + 0.5)  

xlabel('Parameter name')

ylabel('Parameter name')

set(gca, 'XtickLabel', Param_names_red, 'YtickLabel', Param_names_red(length(Param_names_red):-1:1), 'Fontsize', 16)


%% Save figure

savefig(Correlation_plot, Save_loc)

Figure_handle = Correlation_plot;

end

