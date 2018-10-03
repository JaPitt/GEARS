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


function [Figure_handle] = Plot_regularised_parameter_summary_GEARS(Results, Save_loc)
% A function that plots a uitable of the final parameters results summary.
% Results  - The results structure from GEARS. (structure)
% Save_loc - Where the figure should be saved.


%% Check inputs

Ins = who;

    if sum(ismember({'Results', 'Save_loc'}, Ins)) ~= 2
        
    error('Both the inputs are required to run this function')           
        
    end
    

%% Create_table

Table_figure = figure('color', 'w');

Table_figure.Position = [0 0 1250 750];

hold all

Txt_title = uicontrol('Style', 'text', 'Position', [500 675 300 50], 'String', ['GEARS results summary']);

Txt_title.BackgroundColor = 'w';

Txt_title.FontSize = 17.6; % Same as the other titles ie 16*Matlabs title factor

Txt_title.FontName = 'Helvetica'; % The normal matlab figure font.

Txt_title.FontWeight = 'bold'; % Like all other matlab titles

Column_names = {'Parameter', 'Value', 'Confidence (95%)', 'Coeff of var (%)', 'Bounds status'};

Table_to_plot = Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary;

    for i = 1:size(Table_to_plot, 1)

        for j = [2, 4]

        Table_to_plot{i, j} = num2str(Table_to_plot{i, j});

        end

    end

Table_handle = uitable('Data', Table_to_plot, 'ColumnName', Column_names);


%% Format table

Table_handle.Position = [20 55 1205 640];

Table_handle.FontSize = 11;

set(gca, 'Visible', 'off')

htmlstart = '<html><h3>'; %html start

htmlend = '</h3></html>'; %html end
    
Column_names = cellfun(@(x)[htmlstart x htmlend], Column_names, 'uni', false); % with html

set(Table_handle,'ColumnName', Column_names)

set(Table_handle,'ColumnWidth', {150, 150, 150, 150, 150});

My_green = [153	204	0]/255;

Table_handle.BackgroundColor(1, :) = My_green;

Table_handle.BackgroundColor(2, :) = [1 1 1];


%% Add condition number text

    if Results.Settings_used.GEARS_options.Calculate_parameter_confidence

        if sum(strcmp('xbest', fieldnames(Results.Global_parameter_estimation.Regularised_estimation))) ~= 0 % If we have regularised results

        Cond = Results.Global_parameter_estimation.Regularised_estimation.FIM_cond_num;

        else

        Cond = Results.Global_parameter_estimation.Non_regularised_estimation.FIM_cond_num;
        
        end

        if Cond > 10^5 % 5 = condition number

        Colour = [1 0 0]; % Red

        Cond_num_message = '. It is highly likely that a lack of identifiability exists here. Metrics calulated using the FIM are probably artefacts.';

        else 

        Colour = [0, 0.8, 0]; % Green  

        Cond_num_message = '. The FIM is not singular.';

        end

    Cond_num_text = annotation('textbox', [0.01, 0.025, 0.8, 0.04], 'String', ['The FIM''s condition number is: ' num2str(Cond) Cond_num_message], 'FitBoxToText', 'off');

    Cond_num_text.LineStyle = 'none';

    Cond_num_text.Color = Colour;

    end


%% Save table

savefig(Table_figure, Save_loc);

Figure_handle = Table_figure;

end
