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


function [Results, Processed_fitting_data, Processed_validation_data] = Run_GEARS(Results_folder, Model_information, Fitting_data, Validation_data, Alternative_options_script)
% The master function that runs the GEARS procedure. This runs the actual procedure.
% Results_folder   - A folder in which all results shall be saved. (string)
% Model_infomation - A structure containg the model infomation that will be used in the GEARS procedure.
    % .Model_name         - The name of the model. Will be used for identification. (string)
    % .Diff_eqns          - The model equations. Should be of the format {'dx = y^2'; 'dy = x'}. Equations not begining with d will be considered non-differential. (cell array) 
    % .Var_names          - The names of the states in the differential equations. (cell array)
    % .Fitting_params     - The names of the parameters whos value should be estimated. (cell array)
    % .Fixed params       - The names of fixed parameters. The value should be defined in data structures in the same order. (cell array) Optional
    % .Param_bounds_lower - The inital lower bounds for the fitting parameters should be in the same order as Fitting_params. (vector)
    % .Param_bounds_upper - The inital upper bounds for the fitting parameters should be in the same order as Fitting_params. (vector)
% Fitting_data - A structure containing the infomation regarding each
% Fitting_data."experiment name" holds the infomation for each experiment. Where the user names the experiment eg Data.exp1. (structure)
    % .Time_points        - The time point of each measurement (vector)
    % .Measurements       - Contains the measurements. Should be of dimensions # Time points by # Observables. Missing points can be denoted as NaN. (matrix)
    % .Standard_deviation - Contains the standard deviation of each point. Should be equal size to .Measurements. (matrix)
    % .Observables        - Contains the names of the observed states. If the observable is not solely a variable set as {'Custom_function'} and use .Custom_function (cell array).
    % .Custom_function    - For observable functions that are not solely states {'x1^2 - 1'; 'x2 - x1'} for example. (cell array) Optional
    % .Initial_conditions - The initial conditions of the experiment for each state. Please ensure the vector is in the same order as Var_names. (vector)
    % .Fixed_params       - The values of Fixed parameters. Should be in the same order as any Fixed_params (vector) Optional
% Validation_data            - Data structure that will be used for validation purposes. Exactly the same format as Fitting data should be used. If left out no model validation will be done (structure).  Optional
% Alternative_options_script - The path to an alternative options script. This script must include all the same options as "GEARS_options". (string) Optional


%% Matlab paths

	if isempty(which('amiwrap.m')) % Look for AMICI.

	error('AMICI needs to be in the Matlab path'); 

	end

	if isempty(which('MEIGO.m')) % Look for MEIGO.

	error('MEIGO needs to be in the Matlab path');

	end


%% Check inputs

Ins = who; % What is in the workspace.

Total_time_taken = cputime; % Timer

    if sum(ismember({'Results_folder', 'Model_information', 'Fitting_data'}, Ins)) ~= 3 % Check all required inputs are being used.
        
    error('The first 3 inputs are required')   
        
    end
    
    if ismember('Validation_data', Ins) % Find out if validation been provided.
     
        if ~isempty(Validation_data)

        Val = true;    

        else 

        Val = false;

        warning('You have not provided any validation data. We highly recommended that you provide GEARS with validation data.')

        Processed_validation_data = []; % For the output.

        end
        
    else 
        
    Val = false;
    
    warning('You have not provided any validation data. I highly recommended that you provide GEARS with validation data.')
        
    end
      
    if ~ischar(Results_folder) % Check the results folder is a string.
        
    error('Results_folder should be a string.')
        
    end
    
    if ismember('Alternative_options_script', Ins) % Find out if the default options file should be used or if annother has been provided.
        
        if ~isempty(Alternative_options_script)
            
        Use_alt_opts = true;    
                
        else 
        
        Use_alt_opts = false;    
            
        end        
        
    else
        
    Use_alt_opts = false;              
        
    end
    

%% Run the Options script

    if Use_alt_opts % If the default options should be used or a different options file.
        
    eval(Alternative_options_script)
        
    else 

    GEARS_Options

    end

    if ~isfinite(Max_number_subplots)

    Max_number_subplots = 10^6; % Basically never split up figures, avoids an error for inf or nan.

    end
    
Required_options = {'GEARS_diary','Calculate_NRMSE','Calculate_R2','Calculate_chi2',...
'Calculate_parameter_confidence','Plot_trajectories','Plot_samples','Plot_residuals',...
'Plot_pred_vs_meas','Plot_correlation_matrix','Calculate_parameter_confidence',...
'Plot_trajectories_with_uncertainty','Plot_convergence_curves','Plot_regularised_parameter_summary',...
'Max_number_open_figures','Max_number_subplots','Export_results_to_xls','Export_results_to_html',...
'Create_figure_report','Int_opts','MEIGO_options'}; % All required options

Workspace_with_options = who; % The current workspace.

    if sum(ismember(Required_options, Workspace_with_options)) ~= length(Required_options) % Check that all required options were found in the options file.
    
    Missing_options = Required_options(~ismember(Required_options, Workspace_with_options)); % Which options are missing.
        
    error('The following options were not found in your options script:\n%s', strjoin(Missing_options))    
        
    end
    
    if ismember('Seed_to_use', Workspace_with_options)
        
    Use_seed = true;    
        
    else
        
    Use_seed = false;            
        
    end
    
    if ismember('P_ref', Workspace_with_options)
        
    Use_P_ref = true;    
        
    else
        
    Use_P_ref = false;            
        
    end

    if ~ismember('Run_standard_optimisation', Workspace_with_options)

    Run_standard_optimisation = false;

    end

    if ~ismember('Run_multi_start', Workspace_with_options)

    Run_multi_start = false;

    end

    if Run_multi_start

    Run_standard_optimisation = true; % To stop parameter bounding and regularisation to being applied.

    Meigo_method = 'MULTISTART';

    else

    Meigo_method = 'ESS';

    end

    if Max_number_subplots == 0

    error('It is not possible to plot 0 subplots please change the Max_number_subplots option')

    end

    if ismember('Format_to_use', Workspace_with_options)

        if ~isempty(Format_to_use)

        eval(['format ' Format_to_use])

        end

    end


%% Diary
    
    if GEARS_diary 

    diary GEARS_diary % Start the diary
    
    end
    
    
%% Initial screen ouput

Line_spaces = @(x) fprintf(repmat('\n', [1, x])); % Creates x spaces on the screen output.

Line_spaces(1)

disp('%%%%% GEARS initiated %%%%%')

disp([char(169) ' 2018 Jake Alan Pitt and Julio R. Banga (IIM-CSIC)'])

Line_spaces(1)

    if Run_multi_start

    disp('Please note as the setting "Run_multi_start" is on global optimisation will not be used.')

    Line_spaces(1)

    end

    if Run_standard_optimisation

    disp('Please note as the setting "Run_standard_optimisation" is on no parameter bounding or regularisation will be performed.')

    Line_spaces(1)

    end


%% Create results folders

    if exist(Results_folder, 'dir') == 0 % If the results folder doesn't exist, create it.
        
    mkdir(Results_folder)   
            
    end
    
	if exist([Results_folder filesep 'Figures'], 'dir') == 0 % If some plotting is required and the figure folder doesn't exist, create it.
               
        if Plot_trajectories + Plot_samples + Plot_residuals + Plot_pred_vs_meas + Plot_correlation_matrix  + ...
        Plot_bounds_param_confidence + Plot_trajectories_with_uncertainty + Plot_convergence_curves + ...  
        Plot_regularised_parameter_summary ~= 0 

        mkdir([Results_folder filesep 'Figures']);    
        
        end

	end    

    if exist([Results_folder filesep 'Reports'], 'dir') == 0 % If the results markup folder doesn't exist and markup reports will be created, create it.

        if Export_results_to_xls + Export_results_to_html + Create_figure_report > 0

        mkdir([Results_folder filesep 'Reports'])

        end

    end
    

%% Check model infomation
% Reshape vectors and check input format.

    if ~ischar(Model_information.Model_name)
        
    error('Model_infomation.Model_name should be a string.')
        
    end
    
Info_fields = fieldnames(Model_information);

    if sum(ismember({'Model_name', 'Diff_eqns', 'Var_names', 'Fitting_params', 'Param_bounds_lower', 'Param_bounds_upper'}, Info_fields)) ~= 6
        
    error('Not all of the required model infomation has been provided')   
        
    end
    
Model_name = Model_information.Model_name;

Model_information.Param_bounds_lower = reshape(Model_information.Param_bounds_lower, [1 length(Model_information.Param_bounds_lower)]);

Model_information.Param_bounds_upper = reshape(Model_information.Param_bounds_upper, [1 length(Model_information.Param_bounds_upper)]);

Diff_eqns = reshape(Model_information.Diff_eqns, [length(Model_information.Diff_eqns), 1]);

Var_names = reshape(Model_information.Var_names, [length(Model_information.Var_names), 1]);

Fitting_params = reshape(Model_information.Fitting_params, [length(Model_information.Fitting_params), 1]);

    if ismember('Fixed_params', Info_fields)
        
    Fixed_params = reshape(Model_information.Fixed_params, [length(Model_information.Fixed_params), 1]);
        
    else
        
    Fixed_params = [];  
        
    end
    
    if sum(isfinite(Model_information.Param_bounds_lower)) ~= length(Model_information.Param_bounds_lower)
        
    error('All lower parameter bounds must be finite.')    
        
    end
    
    if sum(isfinite(Model_information.Param_bounds_upper)) ~= length(Model_information.Param_bounds_upper)
        
    error('All upper parameter bounds must be finite.')        
        
    end
    
    if sum(Model_information.Param_bounds_lower > 0) ~= length(Model_information.Param_bounds_lower)
        
    error('All lower parameter bounds must be strickly positive.')    
        
    end
    
    if sum(Model_information.Param_bounds_upper > 0) ~= length(Model_information.Param_bounds_upper)
        
    error('All upper parameter bounds must be strickly positive.')        
        
    end
    
    if sum(Model_information.Param_bounds_upper > Model_information.Param_bounds_lower) ~= length(Model_information.Param_bounds_upper)
        
    error('All upper bounds must be strickly greater than the lower bound for the same parameter.')        
        
    end

Param_bounds_lower = Model_information.Param_bounds_lower; % The parameter estimation will be performed in log scale.

Param_bounds_upper = Model_information.Param_bounds_upper; % The parameter estimation will be performed in log scale.


%% Create the AMICI files used to solve the internal IVP problem.

Parameters_in_eqns = Create_AMICI_files_GEARS(Model_name, Diff_eqns, Var_names, Fitting_params, [Results_folder filesep 'AMICI_files'], Fixed_params); % Creates all the files needed for solving the IVP problem with AMICI.
    
addpath(genpath([Results_folder filesep 'AMICI_files'])) % The AMICI files created are then added to the path.


%% Screen output

Line_spaces(1)

disp('Processing data ...')


%% Process data

[Processed_fitting_data, Param_initials_f] = Initialise_GEARS_data(Fitting_data, Var_names, Fitting_params, Parameters_in_eqns, Fixed_params); % Formats the Fitting data.

    if Val
        
    [Processed_validation_data, Param_initials_v] = Initialise_GEARS_data(Validation_data, Var_names, Fitting_params, Parameters_in_eqns, Fixed_params); % Formats any Validation data.
    
    Param_initials = Param_initials_f + Param_initials_v > 0; % If any of the initial conditions depend on parameters.

    else

    Param_initials = Param_initials_f; 

    end
    
    
%% Set up simulation function

Simulate = Generate_simulation_handle_GEARS(Model_name, Fixed_params, Param_initials); % Creates a handle for solving the IVP problem.


%% MEIGO problem definition

problem.f = 'Cost_function_GEARS'; % Use this cost function.

    if ismember(lower(MEIGO_options.local.solver), {'fmincon', 'fminsearch', 'nl2sol', 'lsqnonlin', 'solnp', 'n2fb'}) % The solvers that use high quality sensitivity infomation.

    problem.fjac = 'Jacobian_of_residuals_GEARS'; % Use this jacobian function.

    end

problem.x_L = Param_bounds_lower;

problem.x_U = Param_bounds_upper;

problem.x_0 = (problem.x_U - problem.x_L).*rand(size(problem.x_L)) + problem.x_L; % Random initial point.


%% Estimate completion time

Line_spaces(1)

Has_vtr = false;

    if ismember('vtr', fieldnames(problem))

        if ~isempty(problem.vtr)

        Has_vtr = true;

        end

    end

Has_maxeval = false;

    if ismember('maxeval', fieldnames(MEIGO_options))

        if ~isempty(MEIGO_options.maxeval)

        Has_maxeval = true;

        end

    end

Has_maxtime = false;

    if ismember('maxtime', fieldnames(MEIGO_options))

        if ~isempty(MEIGO_options.maxtime)

        Has_maxtime = true;

        end

    end

    if Has_vtr && ~Has_maxeval && ~Has_maxtime

    disp('Completion time can not be estimated when the convergance criteria is only a vtr.')

    else

        if Has_maxtime

        Expected_time_maxtime = 2*MEIGO_options.maxtime;

        end

        if Has_maxeval

        Expected_time_maxeval = cputime;

            for i = 1:30

            Testing_point = (problem.x_U - problem.x_L).*rand(size(problem.x_L)) + problem.x_L; % Random point.

            Cost_function_GEARS(Testing_point, Processed_fitting_data, Simulate, Int_opts, 0);

            end

        Expected_time_maxeval = cputime - Expected_time_maxeval;

        Expected_time_maxeval = MEIGO_options.maxeval*2*Expected_time_maxeval/30;

        end

        if Has_maxtime && Has_maxeval

        Expected_completion = ceil(min(Expected_time_maxtime, Expected_time_maxeval));

        else

            if Has_maxtime
            
            Expected_completion = ceil(Expected_time_maxtime);

            elseif Has_maxeval

            Expected_completion = ceil(Expected_time_maxeval);

            else

            error('No convergence criteria have been set')

            end

        end

    disp(['The expected time of completion of the GEARS computations is ' num2str(Expected_completion) ' seconds'])

    end


%% Screen ouput

Line_spaces(1)

disp('Performing optimisation ...')


%% Run MEIGO initial parameter estimation

	if Use_seed 

	rng(Seed_to_use) % Start at a requested seed.

	end

Seed = rng;
    
    if ~Run_standard_optimisation

    clear Cost_function_GEARS % To clear any persistent variables from the memory.

    Results_initial_run = MEIGO(problem, MEIGO_options, Meigo_method, Processed_fitting_data, Simulate, Int_opts, true); % Run exploratory estimation.

    else

    Results_initial_run = MEIGO(problem, MEIGO_options, Meigo_method, Processed_fitting_data, Simulate, Int_opts, false); % Run the estimation.

    end


%% Parameter sampling size correction

    if ~Run_standard_optimisation

    [~, ~, ~, Sampling_output] = Cost_function_GEARS(ones(length(problem.x_0), 1), Processed_fitting_data, Simulate, Int_opts, true); % Call the cost function to get the sampling results.

    clear Cost_function_GEARS % To clear the persistent variables from the memory.

    Num_sample_points = Sampling_output.Num_sample_points - 1; % The minus one accounts for the call above
    
    Param_sample_points = real(Sampling_output.Param_sample_points(1:Num_sample_points, :)); % Extract the sampling infomation. The real fixes a bug where tiny complex parts are added to parameter values.

    Param_sample_costs = real(Sampling_output.Param_sample_costs(1:Num_sample_points)); % Extract the sampling infomation

    Initial_cost = min(Param_sample_costs); % Cost solution to the first estimation.

    Dummy = find(Param_sample_costs == Initial_cost); 

    Initial_params = Param_sample_points(Dummy(1), :); % Parameters solution to the first estimation.


%% Automated tuning of regularisation parameters and parameter bounding

        if Use_P_ref

        [~, Cost_cut_offs, problem] = Analyse_samples(Param_sample_points, Param_sample_costs, problem, Initial_params); % Calculate P_ref, Cost_cut_offs and new parameter bounds in 2nd function in Run_GEARS.   

        P_ref = reshape(P_ref, length(P_ref), 1);

        else

        [P_ref, Cost_cut_offs, problem] = Analyse_samples(Param_sample_points, Param_sample_costs, problem, Initial_params); % Calculate P_ref, Cost_cut_offs and new parameter bounds in 2nd function in Run_GEARS.   

        end

    alpha = (median(Cost_cut_offs) - Initial_cost)/((diag(1./(P_ref))*(Initial_params' - P_ref))'*(diag(1./((P_ref)))*(Initial_params' - P_ref))); % Calculate the regularisation parameter alpha.


    %% Run MEIGO 2nd parameter estimation with regularisation and reduced bounds 

    problem.x_0 = (problem.x_U - problem.x_L).*rand(size(problem.x_L)) + problem.x_L; % Random initial point.

    problem.f = 'Cost_function_regularised_GEARS'; % Now use this cost function.

        if ismember(lower(MEIGO_options.local.solver), {'fmincon', 'fminsearch', 'nl2sol', 'lsqnonlin', 'solnp', 'n2fb'})

        problem.fjac = 'Jacobian_of_residuals_regularised_GEARS'; % Now use this jacobian function.

        end

    Results_2nd_run = MEIGO(problem, MEIGO_options, Meigo_method, Processed_fitting_data, Simulate, Int_opts, alpha, P_ref); % Run the final regularised estimation.

	else
    
    Initial_params = Results_initial_run.xbest;

    end
 
Total_time_taken = cputime - Total_time_taken; % Stop the timer after the computations are complete


%% Screen output 

Line_spaces(2)

disp('Optimisation complete')

    if Calculate_NRMSE + Calculate_R2 + Calculate_chi2 + Calculate_parameter_confidence > 0
       
    Line_spaces(1)
    
    disp('Performing post fit analyses ...')    
        
    end


%% Calculate NRMSEs

    if Calculate_NRMSE

    [Non_reg_NRMSE_fitting, Non_reg_NRMSE_fitting_exp_wise]= Calculate_NRMSE_GEARS(Initial_params, Processed_fitting_data, Simulate, Int_opts);

        if ~Run_standard_optimisation

        [Reg_NRMSE_fitting, Reg_NRMSE_fitting_exp_wise] = Calculate_NRMSE_GEARS(Results_2nd_run.xbest, Processed_fitting_data, Simulate, Int_opts);
        
        end

        if Val

        [Non_reg_NRMSE_validation, Non_reg_NRMSE_validation_exp_wise]= Calculate_NRMSE_GEARS(Initial_params, Processed_validation_data, Simulate, Int_opts);

            if ~Run_standard_optimisation

            [Reg_NRMSE_validation, Reg_NRMSE_validation_exp_wise] = Calculate_NRMSE_GEARS(Results_2nd_run.xbest, Processed_validation_data, Simulate, Int_opts);    

            end

        end
    
    end
    
    
%% Calculate R2s

    if Calculate_R2

    [Non_reg_R2_fitting, Non_reg_R2_fitting_exp_wise]= Calculate_R2_GEARS(Initial_params, Processed_fitting_data, Simulate, Int_opts);

        if ~Run_standard_optimisation

        [Reg_R2_fitting, Reg_R2_fitting_exp_wise] = Calculate_R2_GEARS(Results_2nd_run.xbest, Processed_fitting_data, Simulate, Int_opts);

        end

        if Val

        [Non_reg_R2_validation, Non_reg_R2_validation_exp_wise]= Calculate_R2_GEARS(Initial_params, Processed_validation_data, Simulate, Int_opts);

            if ~Run_standard_optimisation

            [Reg_R2_validation, Reg_R2_validation_exp_wise] = Calculate_R2_GEARS(Results_2nd_run.xbest, Processed_validation_data, Simulate, Int_opts);    

            end

        end       
                       
    end
    
    
%% Test_chi2

    if Calculate_chi2
        
        if license('test','statistics_toolbox') % The chi2 test relies on the statistics toolbox.

        [Non_reg_chi2_fitting_conclusion, Non_reg_chi2_fitting_p]= Calculate_chi2_GEARS(Initial_params, Processed_fitting_data, Simulate, Int_opts);

            if ~Run_standard_optimisation

            [Reg_chi2_fitting_conclusion, Reg_chi2_fitting_p] = Calculate_chi2_GEARS(Results_2nd_run.xbest, Processed_fitting_data, Simulate, Int_opts);

            end

            if Val

            [Non_reg_chi2_validation_conclusion, Non_reg_chi2_validation_p]= Calculate_chi2_GEARS(Initial_params, Processed_validation_data, Simulate, Int_opts);

                if ~Run_standard_optimisation

                [Reg_chi2_validation_conclusion, Reg_chi2_validation_p] = Calculate_chi2_GEARS(Results_2nd_run.xbest, Processed_validation_data, Simulate, Int_opts);    

                end

            end
        
        else
            
        Calculate_chi2 = false;    
        
        warning('"Calculate_chi2_GEARS" requires the Matlab statistics toolbox. This will be skipped.')   
            
        end
        
    end
  
    
%% Calculate parameter confidence

    if Calculate_parameter_confidence

        if ~Run_standard_optimisation
            	       
        [Parameter_confidence_intervals_regularised, Cond_FIM_regularised] = Calculate_parameter_confidence_GEARS(Results_2nd_run.xbest, Processed_fitting_data, Simulate, Int_opts, 1, alpha, P_ref);
    
        CV_regularised = (Parameter_confidence_intervals_regularised'*100/1.96)./Results_2nd_run.xbest;

        end
    
    [Parameter_confidence_intervals_non_regularised, Cond_FIM_non_regularised] = Calculate_parameter_confidence_GEARS(Initial_params, Processed_fitting_data, Simulate, Int_opts);
    
    CV_non_regularised = (Parameter_confidence_intervals_non_regularised'*100/1.96)./Initial_params;
        
    end

    
%% Screen output

    if Plot_trajectories + Plot_samples + Plot_residuals + Plot_pred_vs_meas + Plot_correlation_matrix + Plot_bounds_param_confidence + Plot_trajectories_with_uncertainty + Plot_convergence_curves + Plot_regularised_parameter_summary > 0
       
    Line_spaces(1)
    
    disp('Plotting figures ...')
        
    end
    

%% Plot trajectory figures

    if Plot_trajectories

        if ~Run_standard_optimisation

        Param_values_for_plot = [Initial_params' Results_2nd_run.xbest']; % Plot both solutions simultanously.

        Param_set_names = {'Non-regularised', 'Regularised'}; % To label each solution.

        Handle = Plot_trajectories_GEARS(Processed_fitting_data, Var_names, Param_values_for_plot, Param_set_names, Simulate, [Results_folder filesep 'Figures' filesep 'Fitting_figures'], Max_number_subplots, Int_opts); % Fitting figures

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            if Val

            Handle = Plot_trajectories_GEARS(Processed_validation_data, Var_names, Param_values_for_plot, Param_set_names, Simulate, [Results_folder filesep 'Figures' filesep 'Validation_figures'], Max_number_subplots, Int_opts); % Validation figures    

            Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            end

        else

        Param_values_for_plot = Initial_params'; % Plot both solutions simultanously.

        Param_set_names = {'Fit to data'}; % To label each solution.

        Handle = Plot_trajectories_GEARS(Processed_fitting_data, Var_names, Param_values_for_plot, Param_set_names, Simulate, [Results_folder filesep 'Figures' filesep 'Fitting_figures'], Max_number_subplots, Int_opts); % Fitting figures

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            if Val

            Param_set_names = {'Model prediction'}; % To label each solution.

            Handle = Plot_trajectories_GEARS(Processed_validation_data, Var_names, Param_values_for_plot, Param_set_names, Simulate, [Results_folder filesep 'Figures' filesep 'Validation_figures'], Max_number_subplots, Int_opts); % Validation figures    

            Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            end

        end
                
    end
       

%% Plot samples   

    if ~Run_standard_optimisation

        if Plot_samples

        Handle = Plot_samples_GEARS(Param_sample_points, Param_sample_costs, Cost_cut_offs, Fitting_params, problem, Param_bounds_upper, Param_bounds_lower, [Results_folder filesep 'Figures' filesep 'Sample_plot'], Max_number_subplots);

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        end

    end
    

%% Plot residuals

    if Plot_residuals

    Handle = Plot_residuals_GEARS(Initial_params, Processed_fitting_data, Simulate, [Results_folder filesep 'Figures' filesep 'Non_regularised_fitting_residuals'], Int_opts);

    Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        if ~Run_standard_optimisation

        Handle = Plot_residuals_GEARS(Results_2nd_run.xbest, Processed_fitting_data, Simulate, [Results_folder filesep 'Figures' filesep 'Regularised_fitting_residuals'], Int_opts);

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        end
    
        if Val

        Handle = Plot_residuals_GEARS(Initial_params, Processed_validation_data, Simulate, [Results_folder filesep 'Figures' filesep 'Non_regularised_validation_residuals'], Int_opts);
        
        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            if ~Run_standard_optimisation

            Handle = Plot_residuals_GEARS(Results_2nd_run.xbest, Processed_validation_data, Simulate, [Results_folder filesep 'Figures' filesep 'Regularised_validation_residuals'], Int_opts);    

            Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            end

        end
    
    end

    
%% Plot Pred VS Measure    

    if Plot_pred_vs_meas

    Handle = Plot_pred_vs_meas_GEARS(Initial_params, Processed_fitting_data, Simulate, [Results_folder filesep 'Figures' filesep 'Non_regularised_fitting_pred_vs_meas'], Int_opts);

    Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        if ~Run_standard_optimisation

        Handle = Plot_pred_vs_meas_GEARS(Results_2nd_run.xbest, Processed_fitting_data, Simulate, [Results_folder filesep 'Figures' filesep 'Regularised_fitting_pred_vs_meas'], Int_opts);

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        end

        if Val

        Handle = Plot_pred_vs_meas_GEARS(Initial_params, Processed_validation_data, Simulate, [Results_folder filesep 'Figures' filesep 'Non_regularised_validation_pred_vs_meas'], Int_opts);

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            if ~Run_standard_optimisation

            Handle = Plot_pred_vs_meas_GEARS(Results_2nd_run.xbest, Processed_validation_data, Simulate, [Results_folder filesep 'Figures' filesep 'Regularised_validation_pred_vs_meas'], Int_opts);

            Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            end

        end
        
    end

    
%% Plot correlation matrix 

    if Plot_correlation_matrix

    [Handle, Corr_matrix_initial] = Plot_correlation_matrix_GEARS(Initial_params, Processed_fitting_data, Simulate, Fitting_params, [Results_folder filesep 'Figures' filesep 'Non_regularised_fitting_corr_matrix'], Int_opts);

    Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        if ~Run_standard_optimisation
    
        [Handle, Corr_matrix_2nd] = Plot_correlation_matrix_GEARS(Results_2nd_run.xbest, Processed_fitting_data, Simulate, Fitting_params, [Results_folder filesep 'Figures' filesep 'Regularised_fitting_corr_matrix'], Int_opts, 1, alpha, P_ref);

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        end
    
    end
    
    
%% Plot bounds and parameter confindnce

    if ~Run_standard_optimisation

        if Plot_bounds_param_confidence

            if Calculate_parameter_confidence

            Handle = Plot_bounds_param_confidence_GEARS(Initial_params, Results_2nd_run.xbest, Fitting_params, problem, Param_bounds_upper, Param_bounds_lower, Parameter_confidence_intervals_non_regularised, Parameter_confidence_intervals_regularised, [Results_folder filesep 'Figures' filesep 'Bounds_param_confidence_plot']);

            Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            else

            warning('Plot_bounds_param_confidence_GEARS requires both Plot_bounds_param_confidence and Calculate_parameter_confidence to be true. This will be skipped') 

            end

        end
 
    end


%% Plot trajectories with uncerrtainty figures    
             
    if Plot_trajectories_with_uncertainty

    Handle = Plot_trajectories_with_uncertainty_GEARS(Processed_fitting_data, Initial_params, Var_names, Simulate, [Results_folder filesep 'Figures' filesep 'Non_regularised_fit_with_uncertainty'], Max_number_subplots, Int_opts); 
   
    Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        if ~Run_standard_optimisation

        Handle = Plot_trajectories_with_uncertainty_GEARS(Processed_fitting_data, Results_2nd_run.xbest, Var_names, Simulate, [Results_folder filesep 'Figures' filesep 'Regularised_fit_with_uncertainty'], Max_number_subplots,Int_opts); 

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        end
    
        if Val
            
        Handle = Plot_trajectories_with_uncertainty_GEARS(Processed_validation_data, Initial_params, Var_names, Simulate, [Results_folder filesep 'Figures' filesep 'Non_regularised_validation_with_uncertainty'], Max_number_subplots,Int_opts); 
   
        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

            if ~Run_standard_optimisation

            Handle = Plot_trajectories_with_uncertainty_GEARS(Processed_validation_data, Results_2nd_run.xbest, Var_names, Simulate, [Results_folder filesep 'Figures' filesep 'Regularised_validation_with_uncertainty'], Max_number_subplots,Int_opts); 

            Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.    

            end
    
        end    
    
    end
    
    
%% Plot convergence curves

    if ~Run_standard_optimisation

        if Plot_convergence_curves

        Handle = Plot_convergence_curves_GEARS(Results_2nd_run, [Results_folder filesep 'Figures' filesep 'Convergence_curves'], Max_number_subplots);

        Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.

        end      
     
    else

        if Plot_convergence_curves

            if ~Run_multi_start

            Handle = Plot_convergence_curves_GEARS(Results_initial_run, [Results_folder filesep 'Figures' filesep 'Convergence_curves'], Max_number_subplots);

            Num_figures_check(Handle, Max_number_open_figures) % If too many figures are open closes some.
            
            else

            warning('Convergence curves will not be plotted for multi-starts.')

            end
        end  

    end


%% Create results structure
% Save all the results into one results structure.

Results.Model_information = Model_information;

    if ~Run_standard_optimisation

    % Parameter_bounding
    Results.Parameter_bounding.Original_bounds.Lower = Param_bounds_lower;

    Results.Parameter_bounding.Original_bounds.Upper = Param_bounds_upper;

    Results.Parameter_bounding.Reduced_bounds.Lower = problem.x_L;

    Results.Parameter_bounding.Reduced_bounds.Upper = problem.x_U;

    % Regularisation
    Results.Regularisation.P_ref = P_ref;

    Results.Regularisation.alpha = alpha;

    % Sampling
    Results.Sampling.Parameter_sample_points = Param_sample_points;

    Results.Sampling.Parameter_sample_costs = Param_sample_costs;

    Results.Sampling.Sample_size = Num_sample_points;

    Results.Sampling.Cost_cut_offs = Cost_cut_offs;

    end

% Global Parameter estimation         
Results.Global_parameter_estimation.Non_regularised_estimation = Results_initial_run;    
        
[~, ~, Results.Global_parameter_estimation.Non_regularised_estimation.Jres] = Jacobian_of_residuals_GEARS(Initial_params, Processed_fitting_data, Simulate, Int_opts);
            
Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary = cell([length(Fitting_params) 4]);

Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 1) = Fitting_params;

Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 2) = cellstr(num2str(Results.Global_parameter_estimation.Non_regularised_estimation.xbest'));

Lower_bounds_active = Results.Global_parameter_estimation.Non_regularised_estimation.xbest' <= Param_bounds_lower' + 10^-3*Param_bounds_lower';
    
Upper_bounds_active = Results.Global_parameter_estimation.Non_regularised_estimation.xbest' >= Param_bounds_upper' - 10^-3*Param_bounds_upper';

    if Calculate_parameter_confidence          

    Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 3) = strcat('+-', cellstr(num2str(Parameter_confidence_intervals_non_regularised))); 
    
    Results.Statistics.Parameter_confidence.Non_regularised_estimation = Parameter_confidence_intervals_non_regularised;   
    
    Results.Statistics.Coefficents_of_variation.Non_regularised_estimation = reshape(CV_non_regularised, [length(CV_non_regularised), 1]);
    
    Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 4) = cellstr(num2str(reshape(CV_non_regularised, [length(CV_non_regularised), 1])));

    Results.Global_parameter_estimation.Non_regularised_estimation.FIM_cond_num = Cond_FIM_non_regularised;
    
    end

Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 5) = repmat({'Bounds not active'}, [length(Results.Global_parameter_estimation.Non_regularised_estimation.xbest), 1]);

Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(Lower_bounds_active, 5) = {'Lower bound active'};

Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(Upper_bounds_active, 5) = {'Upper bound active'};

    if ~Run_standard_optimisation

    Results.Global_parameter_estimation.Regularised_estimation = Results_2nd_run;

    [~, ~, Results.Global_parameter_estimation.Regularised_estimation.Jres] = Jacobian_of_residuals_regularised_GEARS(Results_2nd_run.xbest, Processed_fitting_data, Simulate, Int_opts, alpha, P_ref);

    Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary = cell([length(Fitting_params) 5]);

    Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 1) = Fitting_params;

    Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 2) =  cellstr(num2str(Results.Global_parameter_estimation.Regularised_estimation.xbest'));

	Lower_bounds_active = Results.Global_parameter_estimation.Regularised_estimation.xbest' <= problem.x_L' + 10^-3*problem.x_L';

	Upper_bounds_active = Results.Global_parameter_estimation.Regularised_estimation.xbest' >= problem.x_U' - 10^-3*problem.x_U';

        if Calculate_parameter_confidence

        Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 3) = strcat('+-', cellstr(num2str(Parameter_confidence_intervals_regularised))); 

        Results.Statistics.Parameter_confidence.Regularised_estimation = Parameter_confidence_intervals_regularised;

        Results_dummy.Statistics.Coefficents_of_variation.Regularised_estimation = reshape(CV_regularised, [length(CV_regularised), 1]);

        Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 4) =  cellstr(num2str((reshape(CV_regularised, [length(CV_regularised), 1]))));

        Results.Global_parameter_estimation.Regularised_estimation.FIM_cond_num = Cond_FIM_regularised;

        end

    Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 5) = repmat({'Bounds not active'}, [length(Results.Global_parameter_estimation.Regularised_estimation.xbest), 1]);

    Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(Lower_bounds_active, 5) = {'Lower bound active'};

    Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(Upper_bounds_active, 5) = {'Upper bound active'};

    end

% Statistics

    if Calculate_NRMSE

    Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.All_experiments = Non_reg_NRMSE_fitting;
    
    Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise = Non_reg_NRMSE_fitting_exp_wise;

        if ~Run_standard_optimisation
    
        Results.Statistics.NRMSE.Fitting.Regularised_estimation.All_experiments = Reg_NRMSE_fitting;
    
        Results.Statistics.NRMSE.Fitting.Regularised_estimation.Experiment_wise = Reg_NRMSE_fitting_exp_wise;

        end

        if Val
        
        Results.Statistics.NRMSE.Validation.Non_regularised_estimation.All_experiments = Non_reg_NRMSE_validation;

        Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise = Non_reg_NRMSE_validation_exp_wise;

            if ~Run_standard_optimisation

            Results.Statistics.NRMSE.Validation.Regularised_estimation.All_experiments = Reg_NRMSE_validation;

            Results.Statistics.NRMSE.Validation.Regularised_estimation.Experiment_wise = Reg_NRMSE_validation_exp_wise;

            end
        
        end       

    end
    
    if Calculate_R2

    Results.Statistics.R2.Fitting.Non_regularised_estimation.All_experiments = Non_reg_R2_fitting;
    
    Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise = Non_reg_R2_fitting_exp_wise;

        if ~Run_standard_optimisation
    
        Results.Statistics.R2.Fitting.Regularised_estimation.All_experiments = Reg_R2_fitting;
    
        Results.Statistics.R2.Fitting.Regularised_estimation.Experiment_wise = Reg_R2_fitting_exp_wise;

        end

        if Val
        
        Results.Statistics.R2.Validation.Non_regularised_estimation.All_experiments = Non_reg_R2_validation;

        Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise = Non_reg_R2_validation_exp_wise;

            if ~Run_standard_optimisation

            Results.Statistics.R2.Validation.Regularised_estimation.All_experiments = Reg_R2_validation;

            Results.Statistics.R2.Validation.Regularised_estimation.Experiment_wise = Reg_R2_validation_exp_wise;

            end
        
        end       

    end
    
	if Calculate_chi2

    Results.Statistics.chi2.Fitting.Non_regularised_estimation.Conclusion = Non_reg_chi2_fitting_conclusion;
    
    Results.Statistics.chi2.Fitting.Non_regularised_estimation.Probability = Non_reg_chi2_fitting_p;

        if ~Run_standard_optimisation
    
        Results.Statistics.chi2.Fitting.Regularised_estimation.Conclusion = Reg_chi2_fitting_conclusion;
    
        Results.Statistics.chi2.Fitting.Regularised_estimation.Probability = Reg_chi2_fitting_p;

        end
    
        if Val
        
        Results.Statistics.chi2.Validation.Non_regularised_estimation.Conclusion = Non_reg_chi2_validation_conclusion;

        Results.Statistics.chi2.Validation.Non_regularised_estimation.Probability = Non_reg_chi2_validation_p;

            if ~Run_standard_optimisation

            Results.Statistics.chi2.Validation.Regularised_estimation.Conclusion = Reg_chi2_validation_conclusion;

            Results.Statistics.chi2.Validation.Regularised_estimation.Probability = Reg_chi2_validation_p;

            end
        
        end       

	end
    
    if Plot_correlation_matrix
        
    Results.Statistics.Correlation_matrix.Non_regularised_estimation = Corr_matrix_initial;

        if ~Run_standard_optimisation
    
        Results.Statistics.Correlation_matrix.Regularised_estimation = Corr_matrix_2nd;

        end
    
    end   
    
% Settings used   

Results.Settings_used.GEARS_options.Calculate_NRMSE = Calculate_NRMSE;

Results.Settings_used.GEARS_options.Calculate_R2 = Calculate_R2;

Results.Settings_used.GEARS_options.Calculate_chi2 = Calculate_chi2;

Results.Settings_used.GEARS_options.Calculate_parameter_confidence = Calculate_parameter_confidence;

Results.Settings_used.GEARS_options.Plot_trajectories = Plot_trajectories;

Results.Settings_used.GEARS_options.Plot_samples = Plot_samples;

Results.Settings_used.GEARS_options.Plot_residuals = Plot_residuals;

Results.Settings_used.GEARS_options.Plot_pred_vs_meas = Plot_pred_vs_meas; 

Results.Settings_used.GEARS_options.Plot_correlation_matrix = Plot_correlation_matrix;

Results.Settings_used.GEARS_options.Plot_bounds_param_confidence = Plot_bounds_param_confidence;

Results.Settings_used.GEARS_options.Plot_trajectories_with_uncertainty = Plot_trajectories_with_uncertainty;

Results.Settings_used.GEARS_options.Plot_regularised_parameter_summary = Plot_regularised_parameter_summary;

Results.Settings_used.GEARS_options.Plot_convergence_curves = Plot_convergence_curves;

Results.Settings_used.GEARS_options.Max_number_open_figures = Max_number_open_figures;

Results.Settings_used.GEARS_options.Max_number_subplots = Max_number_subplots;

Results.Settings_used.GEARS_options.Export_results_to_xls = Export_results_to_xls;

Results.Settings_used.GEARS_options.Export_results_to_html = Export_results_to_html;

Results.Settings_used.GEARS_options.Create_figure_report = Create_figure_report;

    if Create_figure_report
        
    Results.Settings_used.GEARS_options.Ghostscript_path = Ghostscript_path;
        
    end

Results.Settings_used.AMICI_options = Int_opts;

Results.Settings_used.MEIGO_options = MEIGO_options;

    if Use_seed

    Results.Settings_used.Seed = Seed_to_use;

    end

    if Use_P_ref

    Results.Settings_used.P_ref = P_ref;
    
    end
    
% Misc
Full_results_path = what(Results_folder); % Full path to the results folder.

Full_results_path = Full_results_path.path; % Full path to the results folder.

Results.Results_folder = ['<a href="matlab: cd(''' Full_results_path '''); ">' Results_folder '</a>']; % Hyperlink to the results folder.

Results.Global_parameter_estimation.Seed = Seed;

Results.Simulation_handle = Simulate;


%% Plot_regularised_parameter_summary

    if ~Run_standard_optimisation

        if Plot_regularised_parameter_summary

        Handle = Plot_regularised_parameter_summary_GEARS(Results, [Results_folder filesep 'Figures' filesep 'Regularised_parameter_summary']);  

        Num_figures_check(Handle, Max_number_open_figures) % If too many figure are open, close some.

        end

    else

    Results_dummy = Results;

    Results_dummy.Global_parameter_estimation.Regularised_estimation.Parameter_summary = Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary;

        if Plot_regularised_parameter_summary

        Handle = Plot_regularised_parameter_summary_GEARS(Results_dummy, [Results_folder filesep 'Figures' filesep 'Non_regularised_parameter_summary']);  

        Num_figures_check(Handle, Max_number_open_figures) % If too many figure are open, close some.

        end

    end

    
%% Screen output

    if Export_results_to_xls + Create_figure_report + Export_results_to_html > 0
        
    Line_spaces(1)
    
    disp('Creating markup reports ...')  
    
    Results.Markup_tag = [Model_name ' ' char(datetime)]; % A tag used for identification.
        
    end
    
    
%% Export to xls

    if Export_results_to_xls

        if Val

        Write_GEARS_results_xls_summary(Results, Processed_fitting_data, [Results_folder filesep 'Reports' filesep 'GEARS_results_summary'], Processed_validation_data);    
              
        else

        Write_GEARS_results_xls_summary(Results, Processed_fitting_data, [Results_folder filesep 'Reports' filesep 'GEARS_results_summary']);    

        end
        
    end
  
    
%% Export to html 

    if Export_results_to_html

        if Val

        Write_GEARS_results_html_summary(Results, Processed_fitting_data, [Results_folder filesep 'Reports' filesep 'GEARS_results_summary_html'], Processed_validation_data, [Results_folder filesep 'Figures']);    

        else

        Write_GEARS_results_html_summary(Results, Processed_fitting_data, [Results_folder filesep 'Reports' filesep 'GEARS_results_summary_html'], [], [Results_folder filesep 'Figures']);    

        end       
        
    end
    
    
%% Create figure report

    if Create_figure_report
                    
    Figure_list = Find_files_in_folders_GEARS([Results_folder filesep 'Figures'], 'fig', 0); % Find all figures in the Results folder
        
    Prefered_order = {'Regularised_parameter_summary', 'Regularised_fit_with_uncertainty', 'Regularised_validation_with_uncertainty', 'Bounds_param_confidence_plot', 'Sample_plot', ...
    'Regularised_fitting_residuals', 'Regularised_validation_residuals', 'Regularised_fitting_pred_vs_meas', 'Regularised_validation_pred_vs_meas' , ...
    'Regularised_fitting_corr_matrix', 'Convergence_curves', 'Fitting_figures', 'Validation_figures', 'Non_regularised_fit_with_uncertainty', 'Non_regularised_validation_with_uncertainty', ...
    'Non_regularised_fitting_residuals', 'Non_Regularised_validation_residuals', 'Non_regularised_fitting_pred_vs_meas', ...
    'Non_Regularised_validation_pred_vs_meas', 'Non_regularised_fitting_corr_matrix'}'; % The best order for the figures to be in the report.

        if Run_standard_optimisation

        Prefered_order{1} = 'Non_regularised_parameter_summary';

        end

    Prefered_order = strcat(Results_folder, filesep, 'Figures', filesep, Prefered_order, '.fig'); % Full path for prefered order.

    In_list = ismember(Figure_list, Prefered_order); % What figures actually are found.

    Other_files = Figure_list(~In_list); % For files not in the prefered order.

    In_list = Figure_list(In_list);

    New_order = cell(size(In_list));

        for i = 1:length(In_list)                

        New_order(ismember(Prefered_order, In_list(i))) = In_list(i); % Reordering positions.

        end   
              
    New_order = [New_order; sort(Other_files)]; % Vector to reorder the figures.
        
    New_order = New_order(~cellfun(@isempty, New_order)); % Reorder the figures.

        if isunix && Plot_regularised_parameter_summary % Is linux and regularised parameter summary is to be plotted

        warning(['Printing of uitables is not possible in linux so the ' char(New_order(1)) ' figure will be left out of the figure report.'])

        New_order = New_order(2:end); % Printing of uitables is not possible in linux.
       
        end 

        if ~isempty(New_order)

            if isunix % Is linux

            Figure_report_GEARS(New_order, [Results_folder filesep 'Reports' filesep 'GEARS_Figure_report'], [], Results.Markup_tag)              

            else
                
            Figure_report_GEARS(New_order, [Results_folder filesep 'Reports' filesep 'GEARS_Figure_report'], Ghostscript_path, Results.Markup_tag)              
        
            end

        end
                
    end
   

%% Delete files created by eSS
% We remove the files created by eSS as they are always recreated in MEIGO anyway.

    if exist('ess_report.mat', 'file') ~= 0

    delete('ess_report.mat') 

    end

    if exist(['objf_' MEIGO_options.local.solver '.m'], 'file') ~= 0

    delete(['objf_' MEIGO_options.local.solver '.m'])

    end

    if exist(['fjac_' MEIGO_options.local.solver '.m'], 'file') ~= 0

    delete(['fjac_' MEIGO_options.local.solver '.m'])
    
    end


%% Save results 

Results.Total_time_taken = Total_time_taken;

    if Val

    save([Results_folder filesep 'GEARS_results'], 'Results', 'Processed_fitting_data', 'Processed_validation_data')

    else
        
    save([Results_folder filesep 'GEARS_results'], 'Results', 'Processed_fitting_data')    
        
    end

clear mex % So the AMICI files can be removed if the user wants to


%% Final output

Line_spaces(1)

    if length(Fitting_params) <= 10; % If less than 10 parameters display the results on screen             

        if ~Run_standard_optimisation
    
        disp(table(Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:,1), Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:,2), Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:,3), ... 
        Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:,4), Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 5), 'variablenames', {'Parameter', 'Value', 'Confidence_95', 'Coeff_of_var', 'Bounds_status'}))            

        else

        disp(table(Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:,1), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:,2), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:,3), ... 
        Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:,4), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 5), 'variablenames', {'Parameter', 'Value', 'Confidence_95', 'Coeff_of_var', 'Bounds_status'}))            

        end

    disp('Please note that the confidence here is calculated using the first order approximation via the FIM and may be inaccurate.')
       
    
    else 
        
    disp('Parameter values will not be show here due to the number of parameters fit')
        
    end

    if Calculate_parameter_confidence

        if ~Run_standard_optimisation

        Cond = Cond_FIM_regularised;
            
        else

        Cond = Cond_FIM_non_regularised;

        end

        if Cond > 10^5

        Cond_num_message = '. It is highly likely that a lack of identifiability exists here. Metrics calulated using the FIM are probably artefacts.';

        else

        Cond_num_message = '. The FIM is not singular.';

        end

    disp(['The FIM''s condition number is: ' num2str(Cond) Cond_num_message])

    end

    if Run_standard_optimisation

    Line_spaces(1)

    disp('Please note that as the setting "Run_standard_optimisation" is on no parameter bounding or regularisation was performed.')    

    end


Line_spaces(2)

disp(['GEARS analysis completed in ' num2str(Total_time_taken) ' seconds'])

Line_spaces(1)

disp('Please find all your results in the results folder requested:')

disp(['<a href="matlab: cd(''' Full_results_path '''); ">' Results_folder ' <-- Use this link to go there.</a>']) % Hyperlink to results

Line_spaces(1)

disp('. .   . .   . .')
disp(' |     |     | ')
disp('|_|   |_|   |_|')  

Line_spaces(2)

disp('%%%%% GEARS task complete %%%%%')

Line_spaces(3)


%% Diary

    if GEARS_diary

    diary off
    
        try

        movefile('GEARS_diary', [Results_folder filesep 'GEARS_diary']) % Move the diary to the results folder.
    
        catch            
        end

    end

rmpath(genpath([Results_folder filesep 'AMICI_files'])) % The AMICI files created are then added to the path.

end

function [P_ref, Cost_cut_offs, problem] = Analyse_samples(Param_sample_points, Param_sample_costs, problem, Initial_params)
% A function that calculates the reference parameter for the regularisation.
% Param_sample_points - The parameter sample points.(matrix)
% Param_sample_costs  - The parameter sample costs. (vector)
% problem             - The MEIGO problem definition. Updates to include reduced bounds. (structure)
% Initial_params      - The solution to the first optimisation step.
    

%% Use different number of boxes for histcounts to search for the cone

Number_boxes = [25, 50, 75, 100, 125, 150];

Cost_cut_offs_each_box = zeros(7, length(problem.x_U));

Cost_cut_offs = zeros(1, size(Cost_cut_offs_each_box, 2));

All_agree = @(Costs_here) 100*abs(median(Costs_here) - mean(Costs_here))./mean(Costs_here) <= 10 ... % Costs said to agree if the mean is similar to the median (within 10%)
& 100*abs(max(Costs_here) - min(Costs_here))/mean(Costs_here) <= 10; % and if the range is small relative to the mean (within 10%)

    for i = 1:size(Cost_cut_offs_each_box, 2) % Number of parameters.                   
        
        for z = 1:7
                                
        Previous_Cut_off_here = inf;

        Reduced_sample_params = (Param_sample_points(:, i)); % Initialise

        Reduced_sample_costs = (Param_sample_costs); % Initialise

        Num_equal_iterations = 0;

        Below_counter = 0;
        
            while Num_equal_iterations ~= 5 % Make sure it has converged
            
                if z ~= 7
                    
                [N, ~, Yedges] = histcounts2((Reduced_sample_params), (Reduced_sample_costs), Number_boxes(z));   

                else
                    
                [N, ~, Yedges] = histcounts2((Reduced_sample_params), (Reduced_sample_costs));   
                    
                end

            N = flipud(N');

            Sig_spread = N./repmat(sum(N, 1), size(N, 1), 1) >= 0.01;

            Size_spread = diff(sum(Sig_spread, 2));
            
            Largest_spread = find(sum(Sig_spread, 2) == max(sum(Sig_spread, 2)), 1, 'last');
                                  
                if max(Largest_spread) == size(Sig_spread, 1)                                                    
                                    
                Cut_off_here = Yedges(2);  
                    
                else   
                         
                Size_spread = Size_spread(Largest_spread:end, :);    
 
                pos = find(Size_spread == min(Size_spread), 1, 'last') + (Largest_spread);               
                    
                Yedges = fliplr(Yedges(2:end));    
                                
                Cut_off_here = Yedges(pos - 1);
                            
                end
                                            
                if Cut_off_here > Previous_Cut_off_here
                                     
                Below_counter = Below_counter + 1;
                                                                   
                end

            Reduced_sample_pos = Param_sample_costs <= Cut_off_here;

            Reduced_sample_params = Param_sample_points(Reduced_sample_pos, i);

            Reduced_sample_costs = Param_sample_costs(Reduced_sample_pos); 
            
                if Previous_Cut_off_here == Cut_off_here || Below_counter >= 3;

                Num_equal_iterations = Num_equal_iterations + 1;    

                end
                
                if Cut_off_here < Previous_Cut_off_here 
            
                Previous_Cut_off_here = Cut_off_here;                            

                end

            end        

        Cost_cut_offs_each_box(z, i) = Cut_off_here;
    
        end


%% Analyse which number of boxes makes most sens

    Failed_boxes = Cost_cut_offs_each_box(:, i) <= min(Param_sample_costs) + 0.1*min(Param_sample_costs);

        if sum(Failed_boxes) == 7 % If all agree that the cone is near the minima (possible lack of identifiabilty)

        Cost_cut_offs(i) = max(Cost_cut_offs_each_box(:, i));

        else

            if All_agree(Cost_cut_offs_each_box(~Failed_boxes, i)) % If all the histcounts agree no matter the box size

            Cost_cut_offs(i) = min(Cost_cut_offs_each_box(~Failed_boxes, i));

            else

            Number_bins = 1;

            Costs_dummy = Cost_cut_offs_each_box(~Failed_boxes, i);

            [Freq, ~, idx] = histcounts(Cost_cut_offs_each_box(~Failed_boxes, i), Number_bins); 

            Sorted_freq = sort(Freq);

                while ~All_agree(Costs_dummy(idx == find(Freq == Sorted_freq(end), 1, 'last')))

                Number_bins = Number_bins + 1;

                [Freq, ~, idx] = histcounts(Cost_cut_offs_each_box(~Failed_boxes, i), Number_bins);

                Sorted_freq = sort(Freq);

                end
                
            Idx_indices = unique(idx);

            Groups_agree = false([length(Idx_indices), 1]);

            Groups_range = zeros([length(Idx_indices), 1]);

                for j = 1:length(Idx_indices)

                Costs_dummy_2 = Costs_dummy(idx == Idx_indices(j));

                Groups_agree(j) = All_agree(Costs_dummy_2);

                Points_meeting_this_condition = Param_sample_points(Param_sample_costs <= min(Costs_dummy_2), i);
                    
                Groups_range(j) = max(Points_meeting_this_condition) - min(Points_meeting_this_condition);

                end

            Groups_range(~Groups_agree) = inf; % Don't use a group that doesn't agree with itself.

            Groups_range_min = find(Groups_range == min(Groups_range));            

            Cost_cut_offs(i) = min(Costs_dummy(idx == Idx_indices(Groups_range_min(1))));

            end
        
        end
   
    end

    
%% Create the p_ref and new parameter bounds   
  
P_ref = zeros(length(Cost_cut_offs), 1);

New_upper_bounds = problem.x_U;

New_lower_bounds = problem.x_L; 

    for i = 1:length(Cost_cut_offs)
        
    Accepted_pos = (Param_sample_costs < Cost_cut_offs(i)); % Parameters with valid cost
    
    Accepted_params = Param_sample_points(Accepted_pos, i);
       
    [N, Edges] = histcounts(Accepted_params, 'normalization', 'probability'); 
    
    Start_pos = find(N == max(N));
    
    End_pos = Start_pos(1) + 1;
    
    Prob = max(N);

        while Prob < 0.9          

            if Start_pos(1) ~= 1 && End_pos ~= length(N) + 1
                
               if sum(N(1:Start_pos - 1)) > sum(N(End_pos + 1:end))
                   
               Start_pos = Start_pos - 1;  
                   
               else
                   
               End_pos = End_pos + 1;      
                   
               end
                    
            else
                
                if Start_pos == 1
                
                End_pos = End_pos + 1; 
                
                elseif End_pos == length(N) + 1
                
                Start_pos = Start_pos - 1;
                
                end    
   
            end
            
        Prob = sum(N(Start_pos:End_pos - 1));    
                        
        end               
                      
    New_upper_bounds(i) = ceil(Edges(End_pos));
    
        if New_upper_bounds(i) > problem.x_U(i) % Don't round outside of original bounds.
            
        New_upper_bounds(i) = problem.x_U(i);
        
        end
    
    New_lower_bounds(i) = floor(Edges(Start_pos));
        
        if New_lower_bounds(i) < problem.x_L(i) % Don't round outside of original bounds.
            
        New_lower_bounds(i) = problem.x_L(i);
        
        end        

    Accepted_params = Accepted_params(Accepted_params >= New_lower_bounds(i) & Accepted_params <= New_upper_bounds(i));
        
    P_ref(i) = median(Accepted_params);

        if All_agree([P_ref(i) Initial_params(i)]) % To avoid p_ref and the initial solution being the same.

        Shift = sqrt(sum((Accepted_params - mean(Accepted_params)).^2)/(length(Accepted_params) - 1))./(max(Accepted_params) - min(Accepted_params)); % Standard deviation relative to the parameter bounds;

            if sum(Accepted_params <= P_ref(i) + Shift) < sum(Accepted_params >= P_ref(i) + Shift)

            Direction = 1;
    
            else

            Direction = -1;

            end

        P_ref(i) = P_ref(i) + Direction*Shift;

        end
    
    end
    
problem.x_L = New_lower_bounds;

problem.x_U = New_upper_bounds;

end

function [] = Num_figures_check(Figure_handle, Max_number_open_figures)
% This function ensures that large numbers of figures are not opened at the same time.
% Figure_handle           - The figure's handle that will be closed if we have a large number of figures already open.
% Max_number_open_figures - The maximum number of figures that should be open at the same time.

Figure_array =  findobj('type','figure');

Num_open_figures = length(Figure_array);

    if Num_open_figures > Max_number_open_figures

    close(Figure_handle)
    
%     warning('In order to stop a large number of figures being opened, some figures have been closed. All figures have still been saved in the Results_folder.')
        
    end
    
clear Figure_handle % Used to aviod saving the figures in the final mat file    

end

