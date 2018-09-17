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


function [Figure_handle] = Plot_residuals_GEARS(Param_values, Data, Simulate, Save_loc, Int_opts)
% A function that plots the normalised residuals, both experiment wise and for all experiments combined.
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

Residual_all_exp = [];

Times = [];

    for i = 1:length(Exp_names)
       
    Simulation = Simulate(Param_values, Data.(char(Exp_names(i))), Int_opts, Data.(char(Exp_names(i))).Time_points);    

%% Residual for each experiment

    Norm = repmat(max(Data.(char(Exp_names(i))).Measurements) - min(Data.(char(Exp_names(i))).Measurements), [size(Data.(char(Exp_names(i))).Measurements, 1), 1]);
    
    Norm = Norm(Data.(char(Exp_names(i))).Finite_data);
    
    Residual{i} = Data.(char(Exp_names(i))).Residual_function(Simulation, Data.(char(Exp_names(i))), Param_values)./Norm;
        
    Residual_all_exp = [Residual_all_exp; Residual{i}];
    
        if sum(Data.(char(Exp_names(i))).Res_positions(1, :)) ~= 0
            
            if strcmp('Custom_function', Data.(char(Exp_names(i))).Observables)
                
            T = repmat(Data.(char(Exp_names(i))).Time_points, 1, length(Data.(char(Exp_names(i))).Custom_function));                   
                
            else
                                
            T = repmat(Data.(char(Exp_names(i))).Time_points, 1, length(Data.(char(Exp_names(i))).Observables));    
                
            end
            
        else
            
            if strcmp('Custom_function', Data.(char(Exp_names(i))).Observables)
                
            T = repmat(Data.(char(Exp_names(i))).Time_points(2:end, :), 1, length(Data.(char(Exp_names(i))).Custom_function));                   
                
            else
                                
            T = repmat(Data.(char(Exp_names(i))).Time_points(2:end, :), 1, length(Data.(char(Exp_names(i))).Observables));    
                
            end
                            
        end
     
    Times = [Times; T(Data.(char(Exp_names(i))).Finite_data)];    

    end
    
    
%% Plot residuals for all experiments   
    
    if length(Exp_names) > 1
    
    Residual_plot = gobjects([length(Exp_names) + 1, 1]);
        
    end

Residual_plot(1) = figure('color', 'w');

Residual_plot(1).Position = [0 0 1250 750];

hold all

title('Residuals for all experiments')

plot([0 max(Times)], [0 0], 'k--', 'Linewidth', 2)

[Times, idx] = sort(Times);

Residual_all_exp = Residual_all_exp(idx);

plot(Times, Residual_all_exp, 'or', 'Linewidth', 2)

xlim([0 max(Times)])

xlabel('Time')

ylabel('Normalised residuals')

box on 

grid on

set(gca, 'FontSize', 16)


%% Plot residuals for individual experiments

    if length(Exp_names) > 1

        for i = 1:length(Exp_names)

        Residual_plot(i + 1) = figure('color', 'w');

        Residual_plot(i + 1).Position = [0 0 1250 750];

        hold all
        
        title(['Residuals for experiment: ' strrep(char(Exp_names(i)), '_', '\_')])
        
        Plot_times = repmat(Data.(char(Exp_names(i))).Time_points, [1 size(Data.(char(Exp_names(i))).Measurements, 2)]);

            if length(Data.(char(Exp_names(i))).Finite_data) == length(Residual{i})

            Plot_times = Plot_times(Data.(char(Exp_names(i))).Finite_data);

            else

            Plot_times = Plot_times(2:end, :);

            Plot_times = Plot_times(Data.(char(Exp_names(i))).Finite_data);

            end
            
        plot(Plot_times, Residual{i}, 'or', 'Linewidth', 2)   
            
        plot([0 max(Plot_times)], [0 0], 'k--', 'Linewidth', 2)

        xlim([0 max(Plot_times)])

        xlabel('Time')

        ylabel('Normalised residuals')

        box on 

        grid on   
        
        set(gca, 'FontSize', 16)

        end

    end

    
%% Save figures

savefig(Residual_plot, Save_loc)

Figure_handle = Residual_plot;

end

