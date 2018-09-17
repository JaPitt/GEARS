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


function [Figure_handle] = Plot_samples_GEARS(Param_sample_points, Param_sample_costs, Cost_cut_offs, Fitting_params, problem, Param_bounds_upper, Param_bounds_lower, Save_loc, Max_number_subplots)
% A function that plots the samples found in the intial optimisation within GEARS.
% Param_sample_points    - A Matrix containing all of the parameter points sampled. Should be size # sample points # parameters. (matrix)
% Param_sample_costs     - A vector containing all the costs of the parameter points sampled. (vector)
% Cost_cut_off           - The value at which parameter points below contribute to the automated regularisation tuning and parameter bounding. (scaler)
% Fitting_params         - The names of the parameters that have been fit. (cell array)
% problem                - The MEIGO formulation of the problem. Please note that this problem should contain the updated parameter bounds (structure)
    % .x_L - The reduced lower bounds. (vector)
    % .x_U - The reduced upper bounds. (vector)
    % Nothing else from this structure is used.
% Param_bounds_upper     - The orginal upper parameters bounds. (vector)
% Param_bounds_lower     - The orginal lower parameters bounds. (vector)
% Save_loc               - The location in which the figure should be saved. (string)
% Max_number_subplots    - The max number of subplots that should be in a single figure. (scaler) Optional


%% Check inputs

Ins = who;

    if sum(ismember({'Param_sample_points', 'Param_sample_costs', 'Cost_cut_offs', 'Fitting_params', 'problem', 'Param_bounds_lower', 'Param_bounds_upper', 'Save_loc'}, Ins)) ~= 8
        
    error('The first 8 inputs are needed to run this function.')   
        
    end
    
    if length(Param_bounds_upper) ~= length(Param_bounds_upper)
        
    error('The number of upper bounds provided must be equal to the number of lower bounds provided')
    
    end
    
    if length(Fitting_params) ~= length(Param_bounds_upper)
        
    error('The number of Fitting parameters must be equal to the number of bounds provided')
    
    end
    
    if ~ischar(Save_loc)
        
    error('The save location should be given as a string')
        
    end
    
    if ~ismember('Max_number_subplots', Ins)
        
    Max_number_subplots = 20; % A feasible maximum number of subplots
    
    end

    if isempty(Max_number_subplots)

    Max_number_subplots = 20; % A feasible maximum number of subplots  

    end
    

%% Plot samples

Num_figs_needed = length(Fitting_params)/Max_number_subplots;

        if round(Num_figs_needed) == Num_figs_needed % Is integer

        Final_fig_num_subfigs = Max_number_subplots;

        else 

        Num_figs_needed = ceil(Num_figs_needed);

        Final_fig_num_subfigs = mod(length(Fitting_params), Max_number_subplots);

        end  

Sample_plot = gobjects([Num_figs_needed, 1]);

    if Num_figs_needed > 1
       
    Num_subplots_per_fig = [repmat(Max_number_subplots, [Num_figs_needed - 1, 1]); Final_fig_num_subfigs];
        
    else
                
    Num_subplots_per_fig = Final_fig_num_subfigs;
        
    end

Counter = 0;

    for i = 1:Num_figs_needed
        
    Sample_plot(i) = figure('color', 'w');  
    
    Sample_plot(i).Position = [0 0 1250 750];
        
        for j = 1:Num_subplots_per_fig(i)
            
        Counter = Counter + 1;     
        
            if Num_subplots_per_fig(i) ~= 1

            subplotrows = floor(sqrt(Num_subplots_per_fig(i))); 

            subplotcolumns = ceil(Num_subplots_per_fig(i)/subplotrows);   

            subplot(subplotrows, subplotcolumns, j)    

            end
        
        hold all 
        
        title(strrep(char(Fitting_params(Counter)), '_', '\_'))

        patch([problem.x_L(Counter) problem.x_U(Counter) problem.x_U(Counter) problem.x_L(Counter)], [min(Param_sample_costs) min(Param_sample_costs) Cost_cut_offs(Counter) Cost_cut_offs(Counter)], 'g')

        plot(Param_sample_points(:, Counter), Param_sample_costs, 'k.')

        xlim([Param_bounds_lower(Counter) Param_bounds_upper(Counter)])

        ylim([min(Param_sample_costs) max(Param_sample_costs)])

        grid on 

        box on
        
        xlabel('Parameter values')
    
        ylabel('Cost')
    
            if j == 1

            legend('Parameter bounds box', 'Sample points')

            end
    
        set(gca, 'Xscale', 'log', 'Yscale', 'log', 'FontSize', 16)  
        
        set(gca, 'YMinorGrid', 'off', 'XMinorGrid', 'off')
                      
        end
    
    end
    
    
%% Save figure
    
savefig(Sample_plot, Save_loc) 

Figure_handle = Sample_plot;

end
    
