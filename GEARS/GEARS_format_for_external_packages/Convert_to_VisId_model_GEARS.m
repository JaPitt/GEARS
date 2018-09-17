function [Case_study, Save_loc] = Convert_to_VisId_model_GEARS(Results, Data, Save_loc)
% A function that converts the model from GEARS format into VisId format. When using VisId leave the AMICI files created in GEARS in the path.
% Results  - Results from GEARS. (structure)
% Data     - The data to analyse the model for, does not require processing first. (structure)
% Save_loc - Where the VisId model should be saved. (string)


%% Inputs check

Ins = who;

    if sum(ismember({'Results', 'Data', 'Save_loc'}, Ins)) ~= 3
        
    error('All inputs are required to run this function')   
        
    end


%% Convert model 

Case_study = [Results.Model_information.Model_name '_model'];

    if ismember('Fixed_params', fieldnames(Results.Model_information))

    inputs.model.par_names = char([Results.Model_information.Fitting_params Results.Model_information.Fixed_params]');

    else

    inputs.model.par_names = char(Results.Model_information.Fitting_params');

    end

inputs.model.n_stimulus = 0;

inputs.model.eqns = char(Results.Model_information.Diff_eqns);

inputs.model.st_names = char(Results.Model_information.Var_names');

inputs.model.n_st = length(Results.Model_information.Var_names);

results.fit.Rjac = Results.Global_parameter_estimation.Regularised_estimation.Jres;

inputs.PEsol.id_global_theta = char(Results.Model_information.Fitting_params');

Exp_names = fieldnames(Data);

inputs.exps.n_exp = 1;

inputs.model.stimulus_names = [];

A = [];

    for i = 1:length(Exp_names)
    
        if ismember('Custom_function', Data.(char(Exp_names(i))).Observables)

        A = [A Data.(char(Exp_names(i))).Custom_function];               

        else

        A = [A Data.(char(Exp_names(i))).Observables];

        end
    
    end

A = unique(A); % All observables over all experiments.

inputs.exps.obs_names{1} = char(strrep(strcat('Obs_', cellstr(num2str((1:length(A))'))), ' ', ''));

inputs.exps.obs{1} = A;

inputs.exps.n_obs{1} = length(A);

save(Save_loc, 'results', 'inputs')

end