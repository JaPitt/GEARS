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


function [] = Write_GEARS_results_html_summary(Results, Fitting_data, Save_loc, Validation_data, Figures_folder)
% A function that writes the results of the GEARS procedure into a html file.
% Results         - The GEARS results should be in the format all GEARS results hold. (structure)
% Fitting_data    - The data that has been fit to. Should be in the processed GEARS format. (structure)
% Save_loc        - The name of where the hmtl file will be saved. Please dont use an extension. (string)
% Validation_data - The data used to validate to. Should be in the processed GEARS format. (structure) [Optional]
% Figures_folder  - The folder containing figures to be included in the html report. [Optional]


%% Check inputs

Ins = who;

    if sum(ismember({'Results', 'Fitting_data', 'Save_loc'}, Ins)) ~= 3
        
    error('The first 3 inputs are required to use this function')   
        
    end 
    
    if ismember('Validation_data', Ins) 
        
        if ~isempty(Validation_data)
        
        Val = 1; 

        else

        Val = 0;    

        Validation_Exp_names = [];
    
        end
        
    else

    Val = 0;    

    Validation_Exp_names = [];

    end

    if ismember('Figures_folder', Ins) 
        
    Include_figs = 1;

    else

    Include_figs = 1;

    end

Fitting_params = reshape(Results.Model_information.Fitting_params, [length(Results.Model_information.Fitting_params) 1]);

   
%% Global parameter estimation

    if ismember('Regularised_estimation', fieldnames(Results.Global_parameter_estimation))

    Reg = 1;

    else

    Reg = 0;

    end

% Non_reg_param_summary_ = table({Results.Markup_tag});

Headers = {'Parameter_value', 'Confidence_95', 'Coeff_of_var_percentage', 'Bounds_status'};

Non_reg_param_summary = table(Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 2), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(: ,3), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 4), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 5), 'variablenames', Headers, 'rownames', Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 1));  

Headers =[Headers 'P_ref'];

    if Reg

    Reg_param_summary = table(Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 2), Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(: ,3), Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 4), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 5), Results.Regularisation.P_ref, 'variablenames', Headers, 'rownames', Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 1));  

    Reg_param_summary.Properties.DimensionNames{1} = 'Parameter_names';

    end

Non_reg_param_summary.Properties.DimensionNames{1} = 'Parameter_names';


%% Statistics

Results_fields = fieldnames(Results);

Input_names_extension = [];

    if ismember('Statistics', Results_fields)
                                
        if Results.Settings_used.GEARS_options.Calculate_NRMSE

            if Reg
            
            Fitting = [Results.Statistics.NRMSE.Fitting.Regularised_estimation.All_experiments Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.All_experiments];

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.NRMSE.Fitting.Regularised_estimation.Experiment_wise  Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_NRMSE_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            else

            Fitting = Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.All_experiments;

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_NRMSE_table = table(Fitting, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            end
        
        assignin('base','Fitting_NRMSE_table', Fitting_NRMSE_table);
        
        Input_names_extension = [Input_names_extension ' Fitting_NRMSE_table'];

            if Val

                if Reg
            
                Validation = [Results.Statistics.NRMSE.Validation.Regularised_estimation.All_experiments Results.Statistics.NRMSE.Validation.Non_regularised_estimation.All_experiments];

                Validation_Exp_names = fieldnames(Validation_data);

                    if length(Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise) > 1

                    Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                    Validation = [Validation; Results.Statistics.NRMSE.Validation.Regularised_estimation.Experiment_wise  Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise];

                    end

                Validation_NRMSE_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 
            
                else

                Validation = Results.Statistics.NRMSE.Validation.Non_regularised_estimation.All_experiments;

                Validation_Exp_names = fieldnames(Validation_data);

                    if length(Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise) > 1

                    Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                    Validation = [Validation; Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise];

                    end

                Validation_NRMSE_table = table(Validation, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                end

            Validation_NRMSE_table.Properties.DimensionNames{1} = 'Experiment';
            
            assignin('base','Validation_NRMSE_table', Validation_NRMSE_table);
        
            Input_names_extension = [Input_names_extension ' Validation_NRMSE_table'];
                                                                    
            else

            assignin('base','Validation_NRMSE_table', []);
        
            Input_names_extension = [Input_names_extension ' Validation_NRMSE_table'];

            end
                
        Fitting_NRMSE_table.Properties.DimensionNames{1} = 'Experiment';
                                                                  
        end
        
        if Results.Settings_used.GEARS_options.Calculate_R2

            if Reg
            
            Fitting = [Results.Statistics.R2.Fitting.Regularised_estimation.All_experiments Results.Statistics.R2.Fitting.Non_regularised_estimation.All_experiments];

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.R2.Fitting.Regularised_estimation.Experiment_wise  Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_R2_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            assignin('base','Fitting_R2_table', Fitting_R2_table);

            Input_names_extension = [Input_names_extension ' Fitting_R2_table'];

            else

            Fitting = Results.Statistics.R2.Fitting.Non_regularised_estimation.All_experiments;

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_R2_table = table(Fitting, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            assignin('base','Fitting_R2_table', Fitting_R2_table);

            Input_names_extension = [Input_names_extension ' Fitting_R2_table'];

            end
        
            if Val

                    if Reg

                    Validation = [Results.Statistics.R2.Validation.Regularised_estimation.All_experiments Results.Statistics.R2.Validation.Non_regularised_estimation.All_experiments];

                    Validation_Exp_names = fieldnames(Validation_data);

                        if length(Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise) > 1

                        Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                        Validation = [Validation; Results.Statistics.R2.Validation.Regularised_estimation.Experiment_wise  Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise];

                        end

                    Validation_R2_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                    Validation_R2_table.Properties.DimensionNames{1} = 'Experiment';

                    assignin('base','Validation_R2_table', Validation_R2_table);

                    Input_names_extension = [Input_names_extension ' Validation_R2_table'];

                    else

                    Validation = Results.Statistics.R2.Validation.Non_regularised_estimation.All_experiments;

                    Validation_Exp_names = fieldnames(Validation_data);

                        if length(Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise) > 1

                        Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                        Validation = [Validation; Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise];

                        end

                    Validation_R2_table = table(Validation, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                    Validation_R2_table.Properties.DimensionNames{1} = 'Experiment';

                    assignin('base','Validation_R2_table', Validation_R2_table);

                    Input_names_extension = [Input_names_extension ' Validation_R2_table'];

                    end
                                            
                else
                                    
                assignin('base','Validation_R2_table', []);

                Input_names_extension = [Input_names_extension ' Validation_R2_table'];

                end

        Fitting_R2_table.Properties.DimensionNames{1} = 'Experiment';
          
        end
        
        if Results.Settings_used.GEARS_options.Calculate_chi2

            if Reg
            
            Fitting = [{num2str(Results.Statistics.chi2.Fitting.Regularised_estimation.Probability)} {num2str(Results.Statistics.chi2.Fitting.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Fitting.Regularised_estimation.Conclusion}, {Results.Statistics.chi2.Fitting.Non_regularised_estimation.Conclusion}];

            Fitting_Exp_names = {'p'; 'Conclusion'};

            Fitting_chi2_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            assignin('base','Fitting_chi2_table', Fitting_chi2_table);

            Input_names_extension = [Input_names_extension ' Fitting_chi2_table'];

            else
            
            Fitting = [{num2str(Results.Statistics.chi2.Fitting.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Fitting.Non_regularised_estimation.Conclusion}];

            Fitting_Exp_names = {'p'; 'Conclusion'};

            Fitting_chi2_table = table(Fitting, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            assignin('base','Fitting_chi2_table', Fitting_chi2_table);

            Input_names_extension = [Input_names_extension ' Fitting_chi2_table'];

            end
        
            if Val

                if Reg
                
                Validation = [{num2str(Results.Statistics.chi2.Validation.Regularised_estimation.Probability)} {num2str(Results.Statistics.chi2.Validation.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Validation.Regularised_estimation.Conclusion} {Results.Statistics.chi2.Validation.Non_regularised_estimation.Conclusion}];

                Validation_Exp_names = {'p'; 'Conclusion'};

                Validation_chi2_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                assignin('base','Validation_chi2_table', Validation_chi2_table);

                Input_names_extension = [Input_names_extension ' Validation_chi2_table'];

                else

                Validation = [{num2str(Results.Statistics.chi2.Validation.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Validation.Non_regularised_estimation.Conclusion}];

                Validation_Exp_names = {'p'; 'Conclusion'};

                Validation_chi2_table = table(Validation, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                assignin('base','Validation_chi2_table', Validation_chi2_table);

                Input_names_extension = [Input_names_extension ' Validation_chi2_table'];

                end
                                            
            else

            assignin('base','Validation_chi2_table', []);

            Input_names_extension = [Input_names_extension ' Validation_chi2_table'];

            end
                                                      
        end
          
    end


%% Parameter bounding

    if Reg

    Parameter_bounds_original = table(Results.Parameter_bounding.Original_bounds.Lower', Results.Parameter_bounding.Original_bounds.Upper', 'variablenames', {'Lower', 'Upper'}, 'rownames', Fitting_params);

    Parameter_bounds_reduced = table(Results.Parameter_bounding.Reduced_bounds.Lower', Results.Parameter_bounding.Reduced_bounds.Upper', 'variablenames', {'Lower', 'Upper'}, 'rownames', Fitting_params);

    Parameter_bounds_original.Properties.DimensionNames{1} = 'Parameter_names';

    Parameter_bounds_reduced.Properties.DimensionNames{1} = 'Parameter_names';


    %% Sampling

    Sample_table = table(Results.Sampling.Sample_size/prod(Results.Parameter_bounding.Original_bounds.Upper - Results.Parameter_bounding.Original_bounds.Lower), ...
    median(Results.Sampling.Cost_cut_offs), Results.Sampling.Sample_size, ...
    'variablenames', {'Sample_dens_approx', 'Avg_Cost_cut_off', 'Sample_size'});

    end
  

%% Publish results 

pub_opt.showCode = false;

pub_opt.outputDir = Save_loc;

r = which('html_report_format_GEARS.m');

Markup_tag = Results.Markup_tag;

assignin('base', 'Markup_tag', Markup_tag);

    if Reg

    alpha1 = Results.Regularisation.alpha;

    assignin('base', 'alpha1', alpha1);

    assignin('base', 'Reg_param_summary', Reg_param_summary);

    assignin('base', 'Parameter_bounds_reduced', Parameter_bounds_reduced);

    assignin('base', 'Parameter_bounds_original', Parameter_bounds_original);

    assignin('base', 'Sample_table', Sample_table);

    end

assignin('base', 'Non_reg_param_summary', Non_reg_param_summary);

    if Include_figs

    Figures_to_include = Find_files_in_folders_GEARS(Figures_folder, '.fig', 0, 1);

    else 

    Figures_to_include = [];

    end

assignin('base', 'Figures_to_include', Figures_to_include);

    if Reg

        if ~isempty(Input_names_extension)

        pub_opt.codeToEvaluate = sprintf('html_report_format_GEARS(%s)',['alpha1, Markup_tag, Reg_param_summary, Non_reg_param_summary, Parameter_bounds_reduced, Parameter_bounds_original, Sample_table, Figures_to_include' strrep(Input_names_extension, ' ', ', ')]);

        else

        pub_opt.codeToEvaluate = sprintf('html_report_format_GEARS(%s)','alpha1, Markup_tag, Reg_param_summary, Non_reg_param_summary, Parameter_bounds_reduced, Parameter_bounds_original, Sample_table, Figures_to_include');

        end

    else

        if ~isempty(Input_names_extension)

        pub_opt.codeToEvaluate = sprintf('html_report_format_GEARS(%s)',['[], Markup_tag, [], Non_reg_param_summary, [], [], [], Figures_to_include' strrep(Input_names_extension, ' ', ', ')]);

        else

        pub_opt.codeToEvaluate = sprintf('html_report_format_GEARS(%s)','[], Markup_tag, [], Non_reg_param_summary, [], [], [], Figures_to_include');

        end

    end

pub_opt.useNewFigure = false;

pub_opt.maxWidth = 700;

pub_opt.createThumbnail = false;

Figures_open_before = findobj('type','figure'); % The handles of all figures open before publishing.

publish(r, pub_opt);

Figures_open_after =  findobj('type','figure'); % The handles of all figures open after publishing

close(Figures_open_after(~ismember(Figures_open_after, Figures_open_before))) % Close figures opened in publishing

% movefile([Save_loc filesep 'html_report_format_GEARS.html'], [Save_loc '.html'])

% rmdir(Save_loc)

Assigned_vars = {'Markup_tag' 'Reg_param_summary' 'alpha1' 'Parameter_bounds_reduced' 'Parameter_bounds_original' 'Sample_table' 'Non_reg_param_summary' ...
'Validation_chi2_table' 'Validation_chi2_table' 'Fitting_chi2_table' 'Fitting_chi2_table' 'Validation_R2_table' ...
'Validation_R2_table' 'Fitting_R2_table' 'Fitting_R2_table' 'Validation_NRMSE_table' 'Fitting_NRMSE_table' 'Figures_to_include'};

evalin('base',['clear ' strjoin(Assigned_vars)] ); % To clear the assignin variables clear the variables from this function


end
    
