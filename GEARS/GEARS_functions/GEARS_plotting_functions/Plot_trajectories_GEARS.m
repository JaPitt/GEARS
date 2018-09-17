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


function [Figure_handle] = Plot_trajectories_GEARS(Data, Var_names, Param_values_for_plot, Param_set_names, Simulate, Save_loc, Max_number_subplots, Int_opts)
% A function that plots the trajectories for multiple paramter vectors for a given data set. 
% Data                  - The set of data for which the plots will be made. Must be in the format created in "Initialise_GEARS_data" (structure)
% Var_names             - The names of the states in the differential equations. (cell array)
% Param_values_for_plot - A matrix containing parameter vectors that should be plotted should be of format # Parameter sets to be plotted # Parameters (matrix)
% Simulate              - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Save_loc              - The location in which the figures will be saved. (string)
% Int_opts              - The intergration options for AMICI (structure) Optional
% Max_number_subplots - The max number of subplots that should be in a single figure. (scaler) Optional


%% Check inputs 

Ins = who;

    if sum(ismember({'Data', 'Var_names', 'Param_values_for_plot', 'Param_set_names', 'Simulate', 'Save_loc'}, Ins)) ~= 6
        
    error('The first 6 inputs are required')       

    end

    if ~ismember('Int_opts', Ins)
     
    Int_opts = [];    
        
    end   
    
    if ~ismember('Max_number_subplots', Ins) 
        
    Max_number_subplots = 20; % A feasible maximum number of subplots  
    
    end

    if isempty(Max_number_subplots)

    Max_number_subplots = 20; % A feasible maximum number of subplots  

    end

    
%% Find how many experiments are being considered

Exp_names = fieldnames(Data);

Simulation = cell([size(Param_values_for_plot, 2), 1]);


%% Simulate and plot figures

Counter1 = 0;

	for i = 1:length(Exp_names)

        if strcmp('Custom_function', Data.(char(Exp_names(i))).Observables)
            
        Num_sub_plots = length(Data.(char(Exp_names(i))).Custom_function);    
                
        Cust = 1;                
            
        else
            
        For_plotting = find(sum(Data.(char(Exp_names(i))).Res_positions) > 0);    
            
        Num_sub_plots = length(Data.(char(Exp_names(i))).Observables);
    
        Cust = 0;    
            
        end
            
        for k = 1:size(Param_values_for_plot, 2)
        
        Simulation{k} = Simulate(Param_values_for_plot(:, k), Data.(char(Exp_names(i))), Int_opts,  0:(Data.(char(Exp_names(i))).Time_points(end))/10000:Data.(char(Exp_names(i))).Time_points(end));   
 
            if Cust 

            Observations{k} = Data.(char(Exp_names(i))).Observation(Simulation{k}, Data.(char(Exp_names(i))), Param_values_for_plot(:, k));
            
            end
        
        end       

    Num_figs_needed = Num_sub_plots/Max_number_subplots;

        if round(Num_figs_needed) == Num_figs_needed % Is an integer

        Final_fig_num_subfigs = Max_number_subplots;

        else 

        Num_figs_needed = ceil(Num_figs_needed);

        Final_fig_num_subfigs = mod(Num_sub_plots, Max_number_subplots);

        end        

        if Num_figs_needed > 1
       
        Num_subplots_per_fig = [repmat(Max_number_subplots, [Num_figs_needed - 1, 1]); Final_fig_num_subfigs];
        
        else
                
        Num_subplots_per_fig = Final_fig_num_subfigs;
        
        end

    Counter2 = 0;

        for j = 1:Num_figs_needed

        Counter1 = Counter1 + 1;

        Fig(Counter1) = figure('color' ,'w'); 
    
        Fig(Counter1).Position = [0 0 1250 750];

            for z = 1:Num_subplots_per_fig(j)

            Counter2 = Counter2 + 1;

                if Num_subplots_per_fig(j) ~= 1

                subplotrows = floor(sqrt(Num_subplots_per_fig(j))); 

                subplotcolumns = ceil(Num_subplots_per_fig(j)/subplotrows);   

                subplot(subplotrows, subplotcolumns, z)    

                end

            hold all
            
                if Cust
                
                    for k = 1:size(Param_values_for_plot, 2) 

                    plot(Simulation{k}.t, Observations{k}(:, Counter2), 'Linewidth', 2)

                    end                                    

                else

                    for k = 1:size(Param_values_for_plot, 2) 

                    plot(Simulation{k}.t, Simulation{k}.x(:, For_plotting(Counter2)), 'Linewidth', 2)

                    end    
            
                end 
   
                if length(Data.(char(Exp_names(i))).Time_points) == size(Data.(char(Exp_names(i))).Measurements, 1)

                    if sum(Data.(char(Exp_names(i))).Standard_deviation(:, Counter2) == 1) ~= length(Data.(char(Exp_names(i))).Standard_deviation(:, Counter2))
            
                    errorbar(Data.(char(Exp_names(i))).Time_points, Data.(char(Exp_names(i))).Measurements(:, Counter2), Data.(char(Exp_names(i))).Standard_deviation(:, Counter2), 'k.', 'Linewidth', 2, 'Markersize', 20)    

                    else 

                    plot(Data.(char(Exp_names(i))).Time_points, Data.(char(Exp_names(i))).Measurements(:, Counter2), 'k.', 'Linewidth', 2, 'Markersize', 20)    

                    end

                else

                    if sum(Data.(char(Exp_names(i))).Standard_deviation(:, Counter2) == 1) ~= length(Data.(char(Exp_names(i))).Standard_deviation(:, Counter2))

                    errorbar(Data.(char(Exp_names(i))).Time_points(2:end), Data.(char(Exp_names(i))).Measurements(:, Counter2), Data.(char(Exp_names(i))).Standard_deviation(:, Counter2), 'k.','Linewidth', 2, 'Markersize', 20)
    
                    else

                    plot(Data.(char(Exp_names(i))).Time_points(2:end), Data.(char(Exp_names(i))).Measurements(:, Counter2), 'k.','Linewidth', 2, 'Markersize', 20)

                    end

                end

                if z == 1

                legend(Param_set_names)

                end

                if Cust

                title(strrep(['Obs ' num2str(Counter2) ' experiment: ' (char(Exp_names(i)))], '_', '\_'))    

                else 

                title(strrep([char(Var_names(For_plotting(Counter2))) ' experiment: ' (char(Exp_names(i)))], '_', '\_'))

                end
        
            xlim([0 Data.(char(Exp_names(i))).Time_points(end)])

            xlabel('Time')

            set(gca, 'FontSize', 16)

            grid on

            box on
            
            end
        
        end

	end               


%% Save figure 
    
savefig(Fig, Save_loc) 

Figure_handle = Fig;
                     
end

