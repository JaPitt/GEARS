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


function [Figure_handle] = Plot_pred_vs_meas_GEARS(Param_values, Data, Simulate, Save_loc, Int_opts)
% A function that plots the predictions vs the measurements, both experiment wise and for all experiments combined.
% Param_values - A vector of the parameter values at which the NRMSE will be calculated. (vector)
% Data         - The set of data for which the NRMSE will be calculated. Must be in the format created in "Initialise_GEARS_data" (structure)
% Simulate     - The simulation handle of the model. Should be in the format created in "Generate_simulation_handle_GEARS". (anonymous function)
% Save_loc     - The location in which the figures should be saved. (string)
% Int_opts     - The intergration options for AMICI (structure) Optional
    

%% Check inputs

Ins = who;

    if sum(ismember({'Param_values', 'Data', 'Simulate', 'Save_loc'}, Ins)) ~= 4
        
    error('The first 4 inputs are required to use this function')   
        
    end   

    if ~ischar(Save_loc)
        
    error('The save location should be given as a string')
        
    end
    
    if ~ismember('Int_opts', Ins)
     
    Int_opts = [];    
        
    end 
    

%% Reshape parameter vector

Param_values = reshape(Param_values, [length(Param_values), 1]);


%% Find how many experiments are being considered

Exp_names = fieldnames(Data);


%% Calculate normalised residuals for each experiment

Predictions_all_exp = [];

Measured_data = [];

    for i = 1:length(Exp_names)
       
    Simulation = Simulate(Param_values, Data.(char(Exp_names(i))), Int_opts, Data.(char(Exp_names(i))).Time_points);    

%% Residual for each experiment

    Norm = repmat(max(Data.(char(Exp_names(i))).Measurements) - min(Data.(char(Exp_names(i))).Measurements), [size(Data.(char(Exp_names(i))).Measurements, 1), 1]);
    
    Norm = Norm(Data.(char(Exp_names(i))).Finite_data);
    
        if strcmp('Custom_function', Data.(char(Exp_names(i))).Observables)
            
        Predictions{i} = Data.(char(Exp_names(i))).Observation(Simulation,Data.(char(Exp_names(i))),Param_values);   
        
            if size(Predictions{i}, 2) ~= size(Data.(char(Exp_names(i))).Finite_data, 2)
                
            Predictions{i} = Predictions{i}(2:end, :);    
                
            end
        
        Predictions{i} = Predictions{i}(Data.(char(Exp_names(i))).Finite_data);
            
        else 
            
        Predictions{i} = Simulation.x(Data.(char(Exp_names(i))).Res_positions);    
            
        end
         
    Predictions_all_exp = [Predictions_all_exp; Predictions{i}];
                   
    Measured_data = [Measured_data; Data.(char(Exp_names(i))).Measurements(Data.(char(Exp_names(i))).Finite_data)];  
            
    end
    
    
%% Plot predictions vs measurements for all experiments   
    
    if length(Exp_names) > 1
    
    Pred_vs_meas_plot = gobjects([length(Exp_names) + 1, 1]);
        
    end

Pred_vs_meas_plot(1) = figure('color', 'w');

Pred_vs_meas_plot(1).Position = [0 0 1250 750];

hold all

title('Predictions vs Measurements for all experiments')

plot(Predictions_all_exp, Measured_data, 'or', 'Linewidth', 2)

xlims = get(gca, 'xlim');

ylims = get(gca, 'ylim');

plot([min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)])], [min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)])], 'k--', 'Linewidth', 2)

xlabel('Predictions')

ylabel('Measurements')

box on 

grid on

set(gca, 'FontSize', 16)    


%% Plot predictions vs measurements

    if length(Exp_names) > 1

        for i = 1:length(Exp_names)

        Pred_vs_meas_plot(i + 1) = figure('color', 'w');

        Pred_vs_meas_plot(i + 1).Position = [0 0 1250 750];

        hold all
        
        title(['Predictions vs Measurements for experiment: ' strrep(char(Exp_names(i)), '_', '\_')])
        
        plot(Predictions{i}, Data.(char(Exp_names(i))).Measurements(Data.(char(Exp_names(i))).Finite_data), 'or', 'Linewidth', 2)
            
        xlims = get(gca, 'xlim');

        ylims = get(gca, 'ylim');

        plot([min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)])], [min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)])], 'k--', 'Linewidth', 2)

        xlabel('Predictions')

        ylabel('Measurements')

        box on 

        grid on

        set(gca, 'FontSize', 16)    

        end

    end

    
%% Save figures

savefig(Pred_vs_meas_plot, Save_loc)

Figure_handle = Pred_vs_meas_plot;

end

