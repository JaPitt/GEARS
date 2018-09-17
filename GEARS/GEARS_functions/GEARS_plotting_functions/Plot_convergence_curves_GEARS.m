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


function [Figure_handle] = Plot_convergence_curves_GEARS(MEIGO_results, Save_loc, Max_number_subplots)
% A function that plots convergence curves from MEIGO results format.
% MEIGO_results       - Results of a MEIGO parameter estimation. Should be in standard MEIGO format. (structure)
% Save_loc            - The location in which the figure will be saved. 
% Max_number_subplots - The max number of subplots that should be on a single figure. This will only effect this function if = 1. (scaler) Optional


%% Check inputs

Ins = who;

    if sum(ismember({'MEIGO_results', 'Save_loc'}, Ins)) ~= 2
        
    error('The first 2 inputs are needed to run this function.')   
        
    end

    if ~ismember('Max_number_subplots', Ins)
        
    Max_number_subplots = 20; % Plots both curves in subfigures
    
    end


%% Plot convergance curves

Convergence_plot = figure('color', 'w');

Convergence_plot.Position = [0 0 1250 750];

    if Max_number_subplots ~= 1

    subplot(1, 2, 1)

    end

hold all

title('Convergence curve cpu time')

stairs(MEIGO_results.time, MEIGO_results.f, 'k', 'Linewidth', 2)

set(gca, 'XScale', 'log', 'Yscale', 'log')

xlabel('cpu time (seconds)')

ylabel('Cost')

xlim([min(MEIGO_results.time) max(MEIGO_results.time)])

grid on 

box on

Axes = gca;

Axes.YMinorGrid = 'off';

Axes.XMinorGrid = 'off';

Axes.FontSize = 16;

    if Max_number_subplots ~= 1

    subplot(1, 2, 2)

    else

    Convergence_plot(2) = figure('color', 'w');
    
    Convergence_plot(2).Position = [0 0 1250 750];

    end

hold all

title('Convergence curve function evaulations')

stairs(MEIGO_results.neval, MEIGO_results.f, 'k', 'Linewidth', 2)

xlim([min(MEIGO_results.neval) max(MEIGO_results.neval)])

xlabel('Function evaluations')

ylabel('Cost')

Axes = gca;

set(gca, 'XScale', 'log', 'Yscale', 'log')

grid on 

box on

Axes.YMinorGrid = 'off';

Axes.XMinorGrid = 'off';

Axes.FontSize = 16;


%% Save figure

savefig(Convergence_plot, Save_loc)

Figure_handle = Convergence_plot;

end

