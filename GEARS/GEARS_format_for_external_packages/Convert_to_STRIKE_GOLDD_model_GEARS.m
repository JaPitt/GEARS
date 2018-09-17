function [] = Convert_to_STRIKE_GOLDD_model_GEARS(Model_information, Data, Save_loc)
% A function that converts the model from GEARS format into STRIKE-GOLDD format.
% Please not that STRIKE-GOLDD should not be used with models containing discontinuities (i.e. heaviside functions) and is not capable of dealing with unknown initial conditions.
% Model_infomation - The model infomation in standard GEARS input format. (structure)
% Data             - The data to analyse the model for, does not require processing first. (structure)
% Save_loc         - Where the STRIKE-GOLDD model should be saved. (string)


%% Inputs check

Ins = who; % The variables in the workspace.

    if sum(ismember({'Model_information', 'Data', 'Save_loc'}, Ins)) ~= 3
        
    error('All inputs are required to run this function')   
        
    end
    

%% States

x = sym(reshape(Model_information.Var_names, length(Model_information.Var_names), 1)); % The variables.


%% Parameters

p = sym(reshape(Model_information.Fitting_params, length(Model_information.Fitting_params), 1)); % The parameters.


%% Inputs

Fieldn = fieldnames(Model_information); % The fields saved in Model_infomation.

    if ismember('Fixed_params', Fieldn)
        
    u = sym(reshape(Model_information.Fixed_params, length(Model_information.Var_names), 1)); % Any inputs.
        
    else 
        
    u = [];   
        
    end
    

%% Outputs

Exp_names = fieldnames(Data); % The experiment names.

Dummy = [];

    for i = 1:length(Exp_names) % All observables over all experiments.
        
        if sum(strcmp({'Custom_function'}, Data.(char(Exp_names(i))).Observables)) ~= 0
        
        Dummy = [Dummy; reshape(Data.(char(Exp_names(i))).Custom_function, [length(Data.(char(Exp_names(i))).Custom_function), 1])];
                        
        else
                       
        Dummy = [Dummy; reshape(Data.(char(Exp_names(i))).Observables, [length(Data.(char(Exp_names(i))).Observables), 1])];    
            
        end
        
    end
    
h = unique(Dummy); % The unique observables from all the experiments.
       

%% Dynamic equations

f = cell(length(Model_information.Diff_eqns), 1);

    for i = 1:length(Model_information.Diff_eqns)
       
    A = strsplit(char(Model_information.Diff_eqns(i)), '='); 
    
    f(i) = A(2); % Take the RHS from the equations.
        
    end

f = sym(f);


%% Initial conditions 

ics = [];

known_ics = zeros(1,numel(x)); 


%% Save

save(Save_loc, 'x', 'h', 'p', 'f', 'u', 'ics', 'known_ics')

end




