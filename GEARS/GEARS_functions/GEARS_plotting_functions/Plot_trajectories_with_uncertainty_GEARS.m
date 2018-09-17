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


function [Figure_handle] = Plot_trajectories_with_uncertainty_GEARS(Data, Param_values, Var_names, Simulate, Save_loc, Max_number_subplots, Int_opts, Is_regularised, alpha, P_ref)
% A function that plots the trajectories for multiple paramter vectors for a given data set. 
% Data         - The set of data for which the plots will be made. Must be in the format created in "Initialise_GEARS_data" (structure)
% Var_names    - The names of the states in the differential equations. (cell array)
% Param_values - The parameter vector that the model will be simulated at. (vector)
% Simulate     - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Save_loc     - The location in which the figures will be saved. (string)
% Int_opts     - The intergration options for AMICI (structure) Optional
% Is_regularised - If this set of parameter values are from the regularised estimation (logical) [Optional]
% alpha          - The alpha regularisation parameter. Only required if Is_regularised (scaler) [Optional]
% P_ref          - The P_ref regularisation parameter. Only required if Is_regularised (vector) [Optional]


%% Check inputs 

Ins = who;

    if sum(ismember({'Data', 'Var_names', 'Param_values', 'Simulate', 'Save_loc'}, Ins)) ~= 5
        
    error('The first 5 inputs are required')       

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

    if ~ismember('Is_regularised', Ins)
     
    Is_regularised = 0;    
        
    end

    
%% Reshape parameter vector

Param_values = reshape(Param_values, [length(Param_values), 1]);    
      

%% Find how many experiments are being considered

Exp_names = fieldnames(Data);

Fig = gobjects([length(Exp_names), 1]);


%% Caculate confidence intervals

Counter1 = 0;

	for i = 1:length(Exp_names)
        
        if Is_regularised

        [~, ~, Jres] = Jacobian_of_residuals_regularised_GEARS(Param_values, Data, Simulate, Int_opts, alpha, P_ref);

        else

        [~, ~, Jres] = Jacobian_of_residuals_GEARS(Param_values, Data, Simulate, Int_opts);

        end

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

    [Cost, ~, Residual] = Cost_function_GEARS(Param_values, Data, Simulate, Int_opts, 0);

    var_cov_matrix(~zero_cols, ~zero_cols) = Cost*(eye(length(red_FIM))/red_FIM)/length(Residual);

        if r_cond_FIM <= 1e-50 || isnan(r_cond_FIM)

        warning('The Fisher Information Matrix is singular')

        Figure_handle = [];

        return

        end       
        
        if strcmp('Custom_function', Data.(char(Exp_names(i))).Observables)
            
        Num_sub_plots = length(Data.(char(Exp_names(i))).Custom_function);    
                
        Cust = 1;                
            
        else
            
        For_plotting = find(sum(Data.(char(Exp_names(i))).Res_positions) > 0);    
            
        Num_sub_plots = length(Data.(char(Exp_names(i))).Observables);
    
        Cust = 0;    
            
        end    
    
    Int_opts.sensi = 1;

    Int_opts.sensi_meth = 'forward';

    Simulation = Simulate(Param_values, Data.(char(Exp_names(i))), Int_opts,  0:(Data.(char(Exp_names(i))).Time_points(end))/1000:Data.(char(Exp_names(i))).Time_points(end));

        if Data.(char(Exp_names(i))).Param_initials % Include the sensitivity from any initial conditions that are dependant on parameters. 

        Sens_initials = Data.(char(Exp_names(i))).Initial_jac_format_function(Data.(char(Exp_names(i))).Initial_jac_scaling_function(Param_values), Simulation.sx(:, :, length(Param_values) + 1 - sum(Data.(char(Exp_names(i))).Initial_params_only):end), Data.(char(Exp_names(i))));

        Simulation.sx(:,:,1:length(Param_values)) = Simulation.sx(:, :, 1:length(Param_values)) + Sens_initials; % Add the sensitivity from initial conditions.            

        end
        
        if Cust

        Sens = zeros(length(Simulation.t), size(Data.(char(Exp_names(i))).Measurements, 2), length(Param_values));   

            for k = 1:size(Data.(char(Exp_names(i))).Jacobian_all_time_str, 3)

                for j = 1:size(Data.(char(Exp_names(i))).Jacobian_all_time_str, 2)     

                Sens(:, j, k) = eval(char(Data.(char(Exp_names(i))).Jacobian_all_time_str(1, j, k)));

                end

            end

        J = permute(Sens, [2, 3, 1]);
         
        else
            
        J = permute(Simulation.sx(:, :, 1:length(Param_values)), [2, 3, 1]);                  
            
        end
        
    Confidence = zeros(length(Simulation.t), length(Var_names));

        for k = 1:length(Simulation.t)

        C = J(:, :, k)*var_cov_matrix*J(:, :, k)';    

        Confidence(k, :) = sqrt((diag(C)))';

        end     
                
        
%% Plot figures        
                    
	Simulation = Simulate(Param_values, Data.(char(Exp_names(i))), Int_opts,  0:(Data.(char(Exp_names(i))).Time_points(end))/1000:Data.(char(Exp_names(i))).Time_points(end));   
    
        if Cust 

        Observations = Data.(char(Exp_names(i))).Observation(Simulation, Data.(char(Exp_names(i))), Param_values);
            
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
                
                y = [Observations(:, Counter2) - sqrt(Confidence(:, Counter2)) Observations(:, Counter2) + sqrt(Confidence(:, Counter2))];
                
                fill([Simulation.t' fliplr(Simulation.t')] , [(y(:, 1))' fliplr(y(:, 2)') ], [0.8, 0.8, 0.8], 'EdgeColor', 'None');

                plot(Simulation.t, Observations(:, Counter2), 'r', 'Linewidth', 2)
                
                else
                    
                x = Simulation.t;

                y1 = Simulation.x(:, For_plotting(Counter2)) - sqrt(Confidence(:, For_plotting(Counter2)));

                y2 = Simulation.x(:, For_plotting(Counter2)) + sqrt(Confidence(:, For_plotting(Counter2)));

                y = [y1 y2];
        
                fill([x' fliplr(x')] , [(y(:, 1))' fliplr(y(:, 2)') ], [0.8, 0.8, 0.8], 'EdgeColor', 'None');

                plot(Simulation.t, Simulation.x(:, For_plotting(Counter2)), 'r', 'Linewidth', 2)
                
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
    
%% Save the figure.    
    
savefig(Fig, Save_loc) 

Figure_handle = Fig;
                     
end

