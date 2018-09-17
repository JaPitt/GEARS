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


%% AMICI options 

Int_opts.atol = 1E-8; % Absolute integration tolerance. (scaler)

Int_opts.rtol = 1E-8; % Relative integration tolerance. (scaler)

Int_opts.maxsteps = 1E5; % Maximum number of integration steps. (scaler)


%% MEIGO options

%%% Convergance options (At least one is required, when more than one is provided MEIGO will stop once one of said criteria is met)

MEIGO_options.maxeval = 25000; % Maximum number of function evaluations. (scaler)

% problem.vtr = []; % The cost value of the optima. If the cost is known this can be used. (scaler) [Optional]

% MEIGO_options.maxtime = []; % Maximum CPU time in seconds. If not finite custom convergence criteria will be used. (scaler) [Optional]
%%%

MEIGO_options.log_var = 1:length(Model_information.Fitting_params); % The indices of variables that should be searched for in log space.

MEIGO_options.log_level = 2; % What should be considered in log space: 1 - Only ref set, 2 - Everything expect the local search, 3 - Everything.

MEIGO_options.local.solver = 'nl2sol'; % Local solver. (string)

MEIGO_options.local.n1 = 1; % Number of iterations before the first local search. (scaler)

MEIGO_options.local.n2 = 10; % Number of iterations between the two local searches. (scaler)

MEIGO_options.local.finish = MEIGO_options.local.solver; % Local solver for the final search. (string)


%% GEARS options

% Run_standard_optimisation = false; % If GEARS should only perform a singular (non-regularised) global optimisation without parameter bounding. This can be used to get insights for the convergence criterion. (logical) [Optional]

% Run_multi_start = false; % If GEARS should only perform a singular (non-regularised) multi-start without parameter bounding. This can be comparisons with global optimisation. (logical) [Optional]

%%% Plotting options

Plot_trajectories = true; % If the trajectories should be plotted. (logical)

Plot_samples = true; % If the parameter samples should be plotted. (logical)

Plot_residuals = true; % If the normalised residuals should be plotted. (logical)

Plot_pred_vs_meas = true; % If the predictions should be plotted against the measurements. (logical)

Plot_correlation_matrix = true; % If the correlation matrix should be calculated. (logical)

Plot_bounds_param_confidence = true; % If the parameter confidence intervals and parameter bounds should be plotted, turn on/off here. Requires Calculate_parameter_confidence to be true. (logical) [Optional]
    
Plot_trajectories_with_uncertainty = true; % % If the traject trajectories should be plotted with uncertainty intervals. (logical) 

Plot_convergence_curves = true; % If convergence curves should be plotted. (logical)

Plot_regularised_parameter_summary = true; % If the summary of the estimation should be created in a uitable, recommended if using Create_figure_report. (logical)

Max_number_open_figures = 0; % The maximum number of figures that will be open at one time (to prevent huge numbers of figures being opened). All figures are still saved. (scaler)

Max_number_subplots = 8; % The maximum number of subplots that should be contained in a single figure. (scaler)

%%% Post analysis options

Calculate_NRMSE = true; % If normalised root mean square error should be calculated. (logical)

Calculate_R2 = true; % If R2 goodness of fit metrics should be calculated. (logical)

Calculate_chi2 = true; % If the chi2 hypothesis should be tested. Requires the statistics toolbox (logical)

Calculate_parameter_confidence = true; % If the confidence for the parameters should be calculated using the FIM (logical)

%%% Markup options

GEARS_diary = true; % If GEARS should save a diary of the procedure. (logical)

Export_results_to_xls = true; % If GEARS results should be exported to xls. (logical)

Export_results_to_html = true; % If GEARS results should be exported to html. (logical)

Create_figure_report = true; % If GEARS should create a pdf report containing all figures. Requires Ghostscript. (logical)

    if Create_figure_report % This option is only needed if Create_figure_report is true
               
    Ghostscript_path = 'C:\Program Files (x86)\gs\gs9.19\bin\gswin32c.exe'; % The path to the Ghostscript exe. Only required for the figure report. Not required for Linux. (string)   
        
    end
    
    
%% Other options

% P_ref = []; % If you already have a specific reference value you want to use. (vector) [Optional]
 
% Seed_to_use = []; % If you want to start the optimisation at a specific seed. (seed) [Optional]

% Format_to_use = 'shorte'; % The display format matlab should use for the GEARS procedure. The default is 'shorte'. (string) [Optional] 

