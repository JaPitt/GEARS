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


function [Parameters_in_eqns] = Create_AMICI_files_GEARS(Model_name, Diff_eqns, Var_names, Fitting_params, Save_loc, Fixed_params)
% A function that generates all the model files required for the use of AMICI. Please ensure that the order of the Var_names is the same as any intial conditions that you give within the GEARS package.
% Model_name     - The name of the model (string)
% Diff_eqns      - The model equations. Should be of the format {'dx = y^2'; 'dy = x'}. Equations not begining with d will be considered non-differential. (cell array) 
% Var_names      - The names of the states in the differential equations. (cell array)
% Save_loc       - The name of the folder in which the files will be saved. (string)
% Fitting_params - The names of any fixed parameters in the system. (cell array) Optional


%% Matlab paths

	if isempty(which('amiwrap.m'))

	error('AMICI needs to be the the Matlab path.');

	end
    
    
%% Check inputs

Ins = who;

    if sum(ismember(Ins, {'Model_name', 'Diff_eqns', 'Var_names', 'Fitting_params', 'Save_loc'})) < 5
      
    error('The first 5 inputs are required to run this function')    
        
    end
    
    if ~ismember({'Fixed_params'}, Ins)
        
    Fixed_params = [];
    
    end
                
    if isempty(Fixed_params)
        
    Fit_all = true;       
        
    else
        
    Fit_all = false;  
    
    Fixed_params = reshape(Fixed_params, [length(Fixed_params), 1]);
        
    end
    

%% Reshape vectors

Diff_eqns = reshape(Diff_eqns, [length(Diff_eqns), 1]);

Fitting_params = reshape(Fitting_params, [length(Fitting_params), 1]);

Var_names = reshape(Var_names, [length(Var_names), 1]);


%% Extract infomation from equations
   
Diff_eqns = strrep(Diff_eqns, ' ', '');    
    
Eqn_split = cellfun(@(x) strsplit(x, '='), Diff_eqns, 'uni', 0);    
    
Eqn_split = cat(1, Eqn_split{:});

LHSs = Eqn_split(:, 1);

RHSs = Eqn_split(:, 2);

Diff_pos = char(LHSs);

Diff_pos = strcmp(cellstr(Diff_pos(:, 1)), 'd');

Alg_eqns = Diff_eqns(~Diff_pos);

LHSs = LHSs(Diff_pos);

RHSs = RHSs(Diff_pos);

Var_names_found = char(LHSs);

Var_names_found = cellstr(Var_names_found(:, 2:end));

Var_pos = zeros(size(Var_names_found));

    for i = 1:length(Var_names_found)
        
    Var_pos(i) = find(strcmp(Var_names_found, Var_names(i)));    
        
    end

LHSs = LHSs(Var_pos);

RHSs = RHSs(Var_pos);


%% Check which parameters appear in the equations.
% These will be checked against the data to see if they appear in the initial conditions.

Parameters_in_eqns = sum(~isAlways(jacobian(sym(RHSs), sym(Fitting_params)) == 0, 'Unknown', 'false')) > 0;


%% Create model file

Model_file_lines = cell(10^6, 1);

Model_file_lines{1} = ['function [model] = ' Model_name '_model()\n'];

Model_file_lines{2} = ['%% The model file for the ' Model_name ' model as created by GEARS for use with AMICI.\n\n\n'];

Model_file_lines{3} = '%%%% Variables\n\n';

Model_file_lines{4} = ['syms t ' strjoin(Var_names') '\n\n'];

Model_file_lines{5} = ['model.sym.x = [' strjoin(Var_names') '];\n\n\n'];      

Model_file_lines{6} = '%%%% Initial conditions\n\n';

Model_file_lines{7} = 'model.sym.x0 = sym(''IC'', (size(model.sym.x)));\n\n\n';

Model_file_lines{8} = '%%%% Parameters and Inputs\n\n';

	if Fit_all
        
	Model_file_lines{9} = ['syms ' strjoin(Fitting_params(Parameters_in_eqns)') '\n\n'];          
        
    else
            
	Model_file_lines{9} = ['syms ' strjoin(Fitting_params(Parameters_in_eqns)') ' ' strjoin(Fixed_params) '\n\n']; 
            
	end
        
Model_file_lines{10} = ['Fitting_params = [' strjoin(Fitting_params(Parameters_in_eqns)') '];\n\n'];
        
	if Fit_all
            
	Model_file_lines{11} = 'model.sym.p = [Fitting_params model.sym.x0];\n\n\n';    
        
	else 
        
	Model_file_lines{11} = ['model.sym.p = [Fitting_params model.sym.x0 ' strjoin(Fixed_params) '];\n\n\n'];    
    
	end
        
    if ~isempty(Alg_eqns)
        
    Model_file_lines{12} = '%%%% Non-differential equations\n\n';
    
        for i = 1:length(Alg_eqns)
    
        Model_file_lines{12 + i} = [char(Alg_eqns(i)) ';\n\n'];
        
        end
        
    Model_file_lines{13 + i}  ='\n';
    
    Line_num = 13 + i + 1;
     
    else
        
    Line_num = 12;
        
    end

Model_file_lines{Line_num} = '%%%% Differential equations\n\n';

Model_file_lines{Line_num + 1} = 'model.sym.xdot = sym(zeros(size(model.sym.x)));\n\n';

Counter = 0;

    for i = 1:length(LHSs)
        
    Counter = Counter + 1;    
       
    Model_file_lines{Line_num + 1 + Counter} = ['%% ' char(LHSs(i)) '\n'];
            
    Counter = Counter + 1;    
    
    Model_file_lines{Line_num + 1 + Counter} = ['model.sym.xdot(' num2str(i) ') = ' char(RHSs(i)) ';\n\n'];    
                          
    end
    
Line_num = Line_num + Counter + 2;
        
Model_file_lines{Line_num}  ='\n';
        
Model_file_lines{Line_num + 1} = '%%%% Observables\n\n';
        
Model_file_lines{Line_num + 2} = 'model.sym.y = [];\n\n';
   
Model_file_lines{Line_num + 3} = 'end\n\n';

Model_file_lines = Model_file_lines(1:Line_num + 3);

    if exist(Save_loc, 'dir') == 0
        
    mkdir(Save_loc)    
        
    end

Z  = fopen([Save_loc filesep Model_name '_model.m'], 'wt');

    for i = 1:length(Model_file_lines)
        
    fprintf(Z, char(Model_file_lines(i)));
                
    end

fclose(Z);


%% Create mex files

addpath(genpath(Save_loc))

    if exist([Save_loc filesep 'simulate_' Model_name '_model'], 'file') == 0 
    
    disp('Creating AMICI files ...')   

        try

        amiwrap([Model_name '_model'], [Model_name '_model'], Save_loc) 

        catch

        error('Something went wrong with the creation of the AMICI files, please ensure you have correctly installed (and tested) AMICI')

        end
        
    end        

end

