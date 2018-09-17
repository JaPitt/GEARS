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


function [Figure_handle] = Plot_bounds_param_confidence_GEARS(Param_values_non_regularised, Param_values_regularised, Fitting_params, problem, Param_bounds_upper, Param_bounds_lower, Parameter_confidence_intervals_non_regularised,Parameter_confidence_intervals_regularised, Save_loc)
% A function that plots a figure visualising the reduction in parameter bounds.
% Param_values_non_regularised                   - The solution vector of the estimation without regularisation. (vector)
% Param_values_regularised                       - The solution vector of the estimation with regularisation. (vector)
% Fitting_params                                 - The names of the parameters that have been fit. (cell array)
% Data                                           - The set of data for which the confidence will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% problem                                        - The MEIGO formulation of the problem. Please note that this problem should contain the updated parameter bounds (structure)
    % .x_L - The reduced lower bounds. (vector)
    % .x_U - The reduced upper bounds. (vector)
    % Nothing else from this structure is used.
% Param_bounds_upper                             - The orginal upper parameters bounds. (vector)
% Param_bounds_lower                             - The orginal lower parameters bounds. (vector)
% Parameter_confidence_intervals_non_regularised - The confidence intervals calculated for the non-regularised estimation. Can be calculated using "Calculate_parameter_confidence_GEARS". (vector)
% Parameter_confidence_intervals_regularised     - The confidence intervals calculated for the regularised estimation.  Can be calculated using "Calculate_parameter_confidence_GEARS". (vector)
% Save_loc                                       - The location in which the figure should be saved. (string)


%% Check inputs

Ins = who;
    
    if length(Param_bounds_upper) ~= length(Param_bounds_upper)
        
    error('The number of upper bounds provided must be equal to the number of lower bounds provided')
    
    end
    
    if length(Fitting_params) ~= length(Param_bounds_upper)
        
    error('The number of Fitting parameters must be equal to the number of bounds provided')
    
    end
    
    if ~ischar(Save_loc)
        
    error('The save location should be given as a string')
        
    end
    
    if sum(ismember({'Param_values_non_regularised', 'Param_values_regularised', 'Fitting_params', 'problem', 'Param_bounds_upper', 'Param_bounds_lower', ...
    'Parameter_confidence_intervals_non_regularised', 'Parameter_confidence_intervals_regularised', 'Save_loc'}, Ins)) ~= 9

    error('All inputs are required to run this function')

    end
    

%% Restrict the number of bars per figure

Max_number_params_each_fig = 20; % If more than 20 parameters the figure will be split across multiple figures to maintain visability

Ind = 1:Max_number_params_each_fig:length(Fitting_params);

    if Ind(end) ~= length(Fitting_params)

    Ind = [Ind length(Fitting_params)];

    end

    for i = 2:length(Ind)

    Param_values_non_regularised_s{i - 1} = Param_values_non_regularised(Ind(i-1):Ind(i));
    
    Param_values_regularised_s{i - 1} = Param_values_regularised(Ind(i-1):Ind(i));

    Fitting_params_s{i - 1} = Fitting_params(Ind(i-1):Ind(i));

    problem_s{i - 1}.x_L = problem.x_L(Ind(i - 1):Ind(i));

    problem_s{i - 1}.x_U = problem.x_U(Ind(i - 1):Ind(i));
    
    Param_bounds_upper_s{i - 1} = Param_bounds_upper(Ind(i - 1):Ind(i));

    Param_bounds_lower_s{i - 1} = Param_bounds_lower(Ind(i - 1):Ind(i));

    Parameter_confidence_intervals_non_regularised_s{i - 1} = Parameter_confidence_intervals_non_regularised(Ind(i - 1):Ind(i));

    Parameter_confidence_intervals_regularised_s{i - 1} = Parameter_confidence_intervals_regularised(Ind(i - 1):Ind(i));

    end


%% Plot bounds and parameter confidence

Bounds_plot = gobjects([length(Ind) - 1, 1]);

    for i = 1:length(Ind) - 1

    Bounds_plot(i) = figure('color', 'w');

    Bounds_plot(i).Position = [0 0 1250 750];

    hold all    

    x = [(1:length(Fitting_params_s{i}))' - 0.25 (1:length(Fitting_params_s{i}))' (1:length(Fitting_params_s{i}))' (1:length(Fitting_params_s{i}))' - 0.25];

    p1 = [problem_s{i}.x_U' problem_s{i}.x_U' problem_s{i}.x_L' problem_s{i}.x_L'];

    p2 = [Param_bounds_upper_s{i}' Param_bounds_upper_s{i}' Param_bounds_lower_s{i}' Param_bounds_lower_s{i}'];

    patch(x', p1', [.5 .8 .5])

    x = x + 0.25;

    patch(x', p2', [.5 .5 0.8])

    h = errorbar((1:length(Fitting_params_s{i})) - 0.125, Param_values_regularised_s{i}, Parameter_confidence_intervals_regularised_s{i}, 'ok', 'Linewidth', 2); 

    h.LData = h.YData - max(problem_s{i}.x_L, h.YData - h.LData); % Restrict confidence to the bounds

    h.UData = min(problem_s{i}.x_U, h.YData + h.UData) - h.YData; % Restrict confidence to the bounds

    h = errorbar((1:length(Fitting_params_s{i})) + 0.125, Param_values_non_regularised_s{i}, Parameter_confidence_intervals_non_regularised_s{i}, 'or', 'Linewidth', 2);   

    h.LData = h.YData - max(Param_bounds_lower_s{i}, h.YData - h.LData); % Restrict confidence to the bounds

    h.UData = min(problem_s{i}.x_U, h.YData + h.UData) - h.YData; % Restrict confidence to the bounds

    legend('Reduced bounds', 'Original bounds', 'Regularised estimate with uncertainty', 'Non-regularised estimate with uncertainty')                  

    ylabel('Parameter bounds value')

    xlabel('Parameter name')

    grid on

    box on

    xlim([0.5 length(Fitting_params_s{i}) + 0.5])

    ylim([min(Param_bounds_lower_s{i})*10^-1, max(Param_bounds_upper_s{i})*10])

    Axes = gca;

    Axes.XTick = 1:length(Fitting_params_s{i});

    Axes.XTickLabel = strrep(Fitting_params_s{i}, '_', '\_');

    Axes.YScale = 'log';

    Axes.FontSize = 16;

    Axes.YMinorGrid = 'off';

    end

            
%% Save figure

savefig(Bounds_plot, Save_loc)

Figure_handle = Bounds_plot;

end

