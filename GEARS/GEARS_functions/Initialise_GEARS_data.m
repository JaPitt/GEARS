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


function [Data, Param_initials] = Initialise_GEARS_data(Data, Var_names, Fitting_params, Parameters_in_eqns, Fixed_params)
% A function that checks the data format and acounts for missing points.
% Data      - A structure containing the infomation regarding each experiment. Data."experiment name" holds the infomation for each experiment. Where the user names the experiment eg Data.exp1.
    % .Time_points        - The time point of each measurement (vector)
    % .Measurements       - Contains the measurements. Should be of dimensions # Time points by # Observables. Missing points can be denoted as nan. (matrix)
    % .Standard_deviation - Contains the standard deviation of each point. Should be equal size to .Measurements. (matrix)
    % .Observables        - Contains the names of the observed states. If the observable is not solely a variable set as {'Custom_function'} and use .Custom_function (cell array).
    % .Custom_function    - For observable functions that are not solely states {'x1^2 - 1'; 'x2 - 69'} for example. (cell array)
    % .Initial_conditions - The initial conditions of the experiment for each state. Please ensure the vector is in the same order as Var_names. (vector)
    % .Fixed_params       - The values of Fixed parameters. Should be in the same order as any Fixed_params (cell array) [Optional]
% Var_names - The variable names of the states. Should be the same order as Data."experiment name".Initial_conditions. (cell array)
% Fitting_params     - The names of the parameters that will be fit. (cell array) 
% Parameters_in_eqns - The logical indication as to if each parametr exists in the model equations (not including initial conditions) (logical vector)
% Fixed_params       - The names of the parameters that are fixed. This input is only required if any initial conditions or observation functions depend on parameters and some parameters are fixed. (cell array) [Optional] 


%% Data checks 

Ins = who;

Exp_names = fieldnames(Data);

    for i = 1:length(Exp_names)

        if sum(abs(size(Data.(char(Exp_names(i))).Measurements) - size(Data.(char(Exp_names(i))).Standard_deviation))) ~= 0
            
        error('The data and standard deiviation matrices must be of equal size. If a measurement is missing use nan to denote this.') 
            
        end
        
        if strcmp('Custom_function', Data.(char(Exp_names(i))).Observables)
            
        	if ~isequal([length(Data.(char(Exp_names(i))).Time_points), length(Data.(char(Exp_names(i))).Custom_function)], size(Data.(char(Exp_names(i))).Measurements))

            error('The data matrix must be # Time points by # Observables in terms of dimensions')    

        	end            
            
        else
        
            if ~isequal([length(Data.(char(Exp_names(i))).Time_points), length(Data.(char(Exp_names(i))).Observables)], size(Data.(char(Exp_names(i))).Measurements))

            error('The data matrix must be # Time points by # Observables in terms of dimensions')    

            end 
        
        end
        
        if strcmp(Data.(char(Exp_names(i))).Observables, 'Custom_function') == 0
           
            if sum(ismember(Data.(char(Exp_names(i))).Observables, Var_names)) ~= length(Data.(char(Exp_names(i))).Observables)
               
            error('Not all the observables were recongnised in the variable names. If you wish to use a function that is not solely a variable please use the {''Custom function''} option')    
                
            end
            
        end
        
        if length(Data.(char(Exp_names(i))).Initial_conditions) ~= length(Var_names)
            
        error('The inital conditions are not of equal size to the variable names. Please also ensure that they are in the correct order.')              
            
        end       
        
    end
    
    if ~ismember('Fitting_params', Ins)
        
    Fitting_params = [];           
        
    end
    
    if ~ismember('Fixed_params', Ins)
        
    Fixed_params = [];           
        
    end
    
Param_initials = 0;

    
%% Data_processing

    for i = 1:length(Exp_names)
        
        if ismember('Fixed_params', fieldnames(Data.(char(Exp_names(i)))))
            
        Data.(char(Exp_names(i))).Fixed_params = reshape(Data.(char(Exp_names(i))).Fixed_params, [length(Data.(char(Exp_names(i))).Fixed_params), 1]);       
            
        end    
        
    Data.(char(Exp_names(i))).Time_points = reshape(Data.(char(Exp_names(i))).Time_points, [length(Data.(char(Exp_names(i))).Time_points), 1]);
    
    Data.(char(Exp_names(i))).Initial_conditions = reshape(Data.(char(Exp_names(i))).Initial_conditions, [1, length(Data.(char(Exp_names(i))).Initial_conditions)]);
                
    Finite_data_pos = isfinite(Data.(char(Exp_names(i))).Measurements);
    
    Data.(char(Exp_names(i))).Finite_data = logical(Finite_data_pos);  
    
        if sum(Data.(char(Exp_names(i))).Standard_deviation(Data.(char(Exp_names(i))).Finite_data) <= 0) ~= 0 
            
        error('All standard deviation data must be larger than zero')
            
        end
    
        if Data.(char(Exp_names(i))).Time_points(1) ~= 0

        Data.(char(Exp_names(i))).Time_points = [0; Data.(char(Exp_names(i))).Time_points];
        
            if strcmp('Custom_function', Data.(char(Exp_names(i))).Observables) == 0  

            Finite_data_pos = [zeros(1, length(Data.(char(Exp_names(i))).Observables)); Finite_data_pos];         

            else

            Finite_data_pos = [zeros(1, length(Data.(char(Exp_names(i))).Custom_function)); Finite_data_pos];
                
            end

        end
        
        if sum(strcmp('Custom_function', Data.(char(Exp_names(i))).Observables)) == 0            
                    
        Data.(char(Exp_names(i))).Res_positions = zeros(length(Data.(char(Exp_names(i))).Time_points), length(Var_names));    
    
            for j = 1:length(Data.(char(Exp_names(i))).Observables)

            Pos = strcmp(Var_names, Data.(char(Exp_names(i))).Observables(j));

            Data.(char(Exp_names(i))).Res_positions(:, Pos) = Finite_data_pos(:, j);            

            end
                        
        Obs_pos = find(ismember(Var_names, Data.(char(Exp_names(i))).Observables));

        Jac_matrix = cell(length(Data.(char(Exp_names(i))).Observables), length(Fitting_params));

        Jac_in = @(Obs_pos) ['Simulation.sx(Data.Res_positions(:,' num2str(Obs_pos) '),' num2str(Obs_pos) ',']; 

        Jacobian_function = [];

            for j = 1:length(Obs_pos)

            Jac_matrix(j, :) = cellstr(strcat(Jac_in(Obs_pos(j)), num2str((1:length(Fitting_params))'), ')'))';

                if j ~= length(Obs_pos)

                Jacobian_function = [Jacobian_function strcat(strjoin(Jac_matrix(j, :)), ';')];    

                else 

                Jacobian_function = [Jacobian_function strjoin(Jac_matrix(j, :))];    

                end

            end

        Jacobian_function = eval(['@(Simulation, Data, Param_values) ' strcat('[', Jacobian_function, ']')]); 
        
        Data.(char(Exp_names(i))).Res_positions = logical(Data.(char(Exp_names(i))).Res_positions);     
            
        Data.(char(Exp_names(i))).Residual_function = @(Simulation, Data, Param_values) (Simulation.x(Data.Res_positions) - Data.Measurements(Data.Finite_data))./Data.Standard_deviation(Data.Finite_data);
        
        Data.(char(Exp_names(i))).Jacobian_function = Jacobian_function;
        
        else
            
        Data.(char(Exp_names(i))) = Custom_function_formater(Data.(char(Exp_names(i))), Finite_data_pos, Var_names, Fitting_params, Fixed_params);  
            
        end
        
        if iscell(Data.(char(Exp_names(i))).Initial_conditions) % For initial conditions dependant on parameters            
        
        Param_initials = Param_initials + 1;

        Data.(char(Exp_names(i))).Param_initials = true;    

        Params_in_intials = sum(~isAlways(jacobian(sym(Data.(char(Exp_names(i))).Initial_conditions), sym(Fitting_params)) == 0, 'Unknown','false')) ~= 0;

            if sum(Params_in_intials | Parameters_in_eqns) ~= length(Parameters_in_eqns)

            error('Parameters were found that do not exist in either the model equations or the initial conditions')

            else

            Data.(char(Exp_names(i))).Initial_params_only = Params_in_intials & ~Parameters_in_eqns;

            end
        
        Data.(char(Exp_names(i))) = Parameter_dep_intials_formater(Data.(char(Exp_names(i))), Fitting_params, Fixed_params);
                  
        else

%             if sum(Parameters_in_eqns) ~= length(Parameters_in_eqns)
%  
%             error('Parameters were found that do not exist in either the model equations or the initial conditions')
%             
%             end
            
        Data.(char(Exp_names(i))).Param_initials = false;    

        Data.(char(Exp_names(i))).Initial_params_only = false(size(Fitting_params));
            
        end

    dummy = ones(size(Data.(char(Exp_names(i))).Measurements));

    pos = sum(Data.(char(Exp_names(i))).Standard_deviation == 1) == size(Data.(char(Exp_names(i))).Standard_deviation, 1);

    dummy(:, pos) = repmat(max(abs(Data.(char(Exp_names(i))).Measurements(:, pos))), [size(Data.(char(Exp_names(i))).Measurements, 1), 1]);
    
    Data.(char(Exp_names(i))).Scaling_vector = 1./dummy(Data.(char(Exp_names(i))).Finite_data);
        
    end
       
Param_initials = Param_initials > 0;
    
end

function [Data] = Custom_function_formater(Data, Finite_data_pos, Var_names, Fitting_params, Fixed_params)
% A function that calculates the residual function and jacobian function for custom observation functions. This function adds the following fields to the data structure: Data.Observation Data.Residual_function Data.Jacobian_function
% Data            - The data for a specific experiment. (structure)
% Finite_data_pos - The position of all finite data values. (matrix)
% Var_names       - The variable names. (cell array)
% Fitting_params  - The names of the parameters to be fit. (cell array)
% Fixed_params    - The names of any fixed parameters. (cell array) [Optional]


%% Check for fixed params
Ins = who;

    if ~ismember('Fixed_params', Ins)
        
    Fixed_params = [];    
        
    end  

    
%% Calculate residual function

Data.Res_positions = zeros(size(Finite_data_pos, 1), length(Var_names), size(Finite_data_pos, 2));

	for j = 1:size(Data.Res_positions, 3)
        
	Data.Res_positions(:, :, j) = repmat(Finite_data_pos(:, j), [1 length(Var_names)]);

	end
    
Data.Res_positions = logical(Data.Res_positions);    

Dummy = subs(Data.Custom_function, reshape(Var_names, [1, numel(Var_names)]), strcat('Dummy(', strrep(cellstr(num2str((1:length(Var_names))')), ' ', ''), ')')');
    
	if ~isempty(Fitting_params)

	Dummy = subs(Dummy, reshape(Fitting_params, [1, numel(Fitting_params)]), strcat('DUMMY(', strrep(cellstr(num2str((1:length(Fitting_params))')), ' ', ''), ')')');

    end
    
	if ~isempty(Fixed_params)

	Dummy = subs(Dummy, reshape(Fixed_params, [1, numel(Fixed_params)]), reshape(Data.Fixed_params, [1, numel(Data.Fixed_params)]));

    end
    
Observation = cell(size(Dummy));    
    
    for j = 1:length(Observation)

	Observation{j} = char(Dummy(j));    
        
	end
            
Observation2 = Observation;
        
Observation2 = strrep(Observation2, 'Dummy(', 'Simulation.x(:, ');
    
	for k = 1:length(Observation)
        
        for j = 1:length(Var_names)

        Observation(k) = strrep(Observation(k), ['Dummy(' num2str(j) ')'], ['Simulation.x(Data.Res_positions(:, ' num2str(j) ', ' num2str(k) '), ' num2str(j) ')']);   

        end
    
    end
            
	if ~isempty(Fitting_params)

	Observation = strrep(Observation, 'DUMMY', 'Param_values');
            
	Observation2 = strrep(Observation2, 'DUMMY', 'Param_values');

    end   

Observation = reshape(Observation, [1, length(Observation)]);           
        
Observation2(1:end-1) = strcat(Observation2(1:end-1), ',');

Observation(1:end-1) = strcat(Observation(1:end-1), ';');

Observation = strcat('[', strjoin(Observation), ']');
        
Observation2 = strcat('[', strjoin(Observation2), ']');

Observation = strrep(Observation, '*', '.*');  

Observation = strrep(Observation, '^', '.^');  

Observation = strrep(Observation, ' ', ''); 

Observation2 = strrep(Observation2, '*', '.*');  

Observation2 = strrep(Observation2, '^', '.^');  

Observation2 = strrep(Observation2, ' ', ''); 


%% Calculate Jacobian function

Custom_function = sym(Data.Custom_function);

Fitting_params_dep = Fitting_params;

Fitting_params_dep(1:end-1) = strcat(Fitting_params_dep(1:end-1), ',');

To_diff = subs(Custom_function, Var_names', strcat(Var_names', ['(' strjoin(Fitting_params_dep) ')']));

Jacobian_eqns = sym(zeros([size(To_diff), length(Fitting_params)]));

    for j = 1:length(Fitting_params)
    
    Jacobian_eqns(:, :, j) = diff(To_diff, char(Fitting_params(j)));    
        
    end
    
Jacobian_eqns = squeeze(Jacobian_eqns);

Sens_vars = cell(length(Fitting_params), length(Var_names));

Dummy_sens_vars = cell(size(Sens_vars));

    for j = 1:length(Var_names)
       
    Sens_vars(:, j) = strcat('diff(', Var_names(j), ['(' strjoin(Fitting_params_dep) '),'], Fitting_params, ')');
    
        for k = 1:length(Fitting_params)
    
        Dummy_sens_vars(k, j) = {['Dummy(' num2str(k) ',' num2str(j) ')']};
        
        end
        
    end

Sens_vars = reshape(Sens_vars, [numel(Sens_vars), 1]);

Dummy_sens_vars = reshape(Dummy_sens_vars, [numel(Dummy_sens_vars), 1]);

Jacobian_eqns = subs(Jacobian_eqns, Sens_vars, Dummy_sens_vars);

dummy_vars = strcat('dummy(', strrep(cellstr(num2str((1:length(Var_names))')), ' ', ''), ')')';

Jacobian_eqns = subs(Jacobian_eqns, reshape(Var_names, [1, numel(Var_names)]), dummy_vars);

DUMMY_params = strcat('DUMMY(', strrep(cellstr(num2str((1:length(Fitting_params))')), ' ', ''), ')')';

Jacobian_eqns = subs(Jacobian_eqns, reshape(Fitting_params, [1, numel(Fitting_params)]), DUMMY_params);

DUMMY_params2 = DUMMY_params;

DUMMY_params2(1:end-1) = strcat(DUMMY_params2(1:end-1), ',');

DUMMY_params2 = strcat('(', strjoin(DUMMY_params2), ')');

Jacobian_eqns = subs(Jacobian_eqns, strcat(dummy_vars, DUMMY_params2), dummy_vars);

    if ~isempty(Fixed_params)
        
    Jacobian_eqns = subs(Jacobian_eqns, Fixed_params', Data.Fixed_params);    
        
    end

Jacobian_str = cell(size(Jacobian_eqns));    
        
	for k = 1:size(Jacobian_eqns, 1)
            
        for z = 1:size(Jacobian_eqns, 2)
            
        Jacobian_str{k, z} = char(Jacobian_eqns(k, z));
            
        end
    
    end
    
Jacobian_str = strrep(Jacobian_str, '*', '.*');  

Jacobian_str = strrep(Jacobian_str, '^', '.^');  

Jacobian_str = strrep(Jacobian_str, ' ', '');  

Jacobian_str = strrep(Jacobian_str, 'DUMMY', 'Param_values');  

Param_idx = repmat((1:length(Fitting_params))', [length(Var_names), 1]);

State_idx = reshape(repmat(1:length(Var_names), [length(Fitting_params), 1]), length(Var_names)*length(Fitting_params), 1);

Simulation_pos_vector = @(obs_num) cellstr(strcat('Simulation.sx(Data.Res_positions(:,', num2str(State_idx), ',', num2str(obs_num), ')', ',', num2str(State_idx), ',', num2str(Param_idx), ')'));

    for j = 1:size(Jacobian_str, 1)        
    
        for z = 1:size(Jacobian_str, 2)
            
        Simulation_pos = Simulation_pos_vector(z);    
    
            for k = 1:length(Simulation_pos)

            Jacobian_str(j, :) =  strrep(Jacobian_str(j, :), Dummy_sens_vars(k), Simulation_pos(k));    

            end 
        
        end
        
        for k = 1:length(Var_names)
            
            for z = 1:size(Jacobian_str, 2)
        
            Jacobian_str(j, :) =  strrep(Jacobian_str(j, :), dummy_vars(k), ['Simulation.x(Data.Res_positions(:,' num2str(k) ',' num2str(z) '),' num2str(k) ')']);
        
            end
       
        end                
        
    end
    
Jacobian_str = Jacobian_str';

Jacobian_function = [];

Jacobian_str_all_t = Jacobian_str;

    for i = 1:size(Jacobian_str, 1)
                     
        if i ~= size(Jacobian_str, 1)

        Jacobian_function = [Jacobian_function strcat(strjoin(Jacobian_str(i, :)), ';')];    

        else 

        Jacobian_function = [Jacobian_function strjoin(Jacobian_str(i, :))];    

        end    
        
        for k = 1:length(Var_names)
        
        Jacobian_str_all_t(i, :) =  strrep(Jacobian_str_all_t(i, :), ['Data.Res_positions(:,' num2str(k) ',' num2str(i) ')'], ':');    
        
        end
                   
    end
    
Jacobian_function = eval(['@(Simulation, Data, Param_values) ' strcat('[', Jacobian_function, ']')]);    

Jacobian_str_all_t = reshape(Jacobian_str_all_t, [1, size(Jacobian_str_all_t, 2), size(Jacobian_str_all_t, 1)]);

Jacobian_str_all_t = permute(Jacobian_str_all_t, [3 1 2]);


%% Update Data structure
        
Data.Observation = eval(['@(Simulation, Data, Param_values) (' Observation2 ')']);

Data.Residual_function = eval(['@(Simulation, Data, Param_values) (' Observation '- Data.Measurements(Data.Finite_data))./Data.Standard_deviation(Data.Finite_data);']);

Data.Jacobian_function = Jacobian_function;

Data.Jacobian_all_time_str = Jacobian_str_all_t;

end

function [Data] = Parameter_dep_intials_formater(Data, Fitting_params, Fixed_params)
% A function that formats the data structure for when the initial conditions depend on parameter values. Adds the following fields to the data structure: Data.Initial_conditions Data.Initial_jac_scaling_function Data.Initial_jac_format_function.
% Data            - The data for a specific experiment. (structure)
% Fitting_params  - The names of the parameters to be fit. (cell array)
% Fixed_params    - The names of any fixed parameters. (cell array) [Optional]


%% Format initials for simulation

Fitting_param_vector = strcat('Param_values(', strrep(cellstr(num2str((1:length(Fitting_params))')), ' ', ''), ')')';

	if ~isempty(Fixed_params)
        
	Call_fun_sym = subs(Data.Initial_conditions, [Fitting_params; Fixed_params]', [Fitting_param_vector, Data.Fixed_params]);
            
    else
                
	Call_fun_sym = subs(Data.Initial_conditions, Fitting_params', Fitting_param_vector);    
            
	end
        
Call_fun = cell(size(Call_fun_sym));

    for j = 1:length(Call_fun)
        
	Call_fun{j} = char(Call_fun_sym(j));

	end

Call_fun(1:end-1) = strcat(Call_fun(1:end-1), ', ');

Call_fun = strcat('[' , strjoin(Call_fun), ']');
                

%% Calculate Jacobian shift

Initials = sym(Data.Initial_conditions);

Initials_eqns = sym(zeros(length(Initials), length(Fitting_params)));

    for j = 1:length(Fitting_params)
        
    Initials_eqns(:, j) = diff(Initials, char(Fitting_params(j)));    
        
    end
    
Fitting_param_vector = strcat('Param_values(', strrep(cellstr(num2str((1:length(Fitting_params))')), ' ', ''), ')')';

Initials_eqns = subs(Initials_eqns, Fitting_params, Fitting_param_vector');

	if ~isempty(Fixed_params)

	Initials_eqns = subs(Initials_eqns, Fixed_params', Data.Fixed_params);

    end
    
Str_initial_eqns = cell(size(Initials_eqns));    

Initials_function = [];    

    for j = 1:size(Initials_eqns, 1)
        
        for k = 1:size(Initials_eqns,2)
                        
        Str_initial_eqns{j, k} = char(Initials_eqns(j, k));    
            
        end
        
        if j~= size(Initials_eqns, 1)
        
        Initials_function = [Initials_function strcat(strjoin(Str_initial_eqns(j, :)), ';')];    
    
        else
            
        Initials_function = [Initials_function strjoin(Str_initial_eqns(j, :))];    
            
        end
        
    end
    
Initial_scaling_function = eval(['@(Param_values) ' strcat('[', Initials_function, ']')]);

Initial_format_function = @(Initials_scaling_matrix, Initials_sens, Data) reshape(sum(repmat(permute(Initials_scaling_matrix, [3 2 1]), [length(Data.Time_points), size(Initials_scaling_matrix, 1), 1]).*repmat(Initials_sens, [1, size(Initials_scaling_matrix, 2), 1]), 3), [length(Data.Time_points), size(Initials_sens, 2), size(Initials_scaling_matrix, 2)]); 


%% Update data structure

Data.Initial_conditions = Call_fun;

Data.Initial_jac_scaling_function = Initial_scaling_function;

Data.Initial_jac_format_function = Initial_format_function;

end

