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


function [] = Figure_report_GEARS(Figure_list, Save_loc, Ghostscript_path, Markup_tag)
% A function that publishes a set of figures into a single pdf (requires ghostscript).
% Figure_list       - The list of figures that should be published. (cell array)
% Save_loc          - The pdf save location. Please do not include the .pdf extension. (string)
% Ghost_script_path - The path of the Ghostscript exe file. Not required for linux (input []). (string) 
% Markup_tag        - A tag used for identification. (string) [Optional]


%% Check inputs

Ins = who;

    if sum(ismember({'Figure_list', 'Save_loc', 'Ghostscript_path'}, Ins)) ~= 3
        
    error('The first 3 inputs are required for this function')   
        
    end

    if isunix

    Is_linux = 1;

    else

    Is_linux = 0;

    end

    if isempty(Ghostscript_path)

    Ghostscript_path = 'dummy'; % To avoid an error in fileparts

    end

	if ismember('Markup_tag', Ins)

        if ~isempty(Markup_tag)
        
        Tag = 1;

        else

        Tag = 0;

        end

    else

    Tag = 0;

	end


%% Create post script with all figures 

    for i = 1:length(Figure_list)
        
    Handle = openfig(char(Figure_list(i)), 'invisible');

        for j = 1:length(Handle)
            
        Subhandle = Handle(length(Handle) - j + 1);    
        
        [~, File_name, ~] = fileparts(Subhandle.FileName);

            if Tag

            ann1 = annotation('textbox', [0 0.97 1 .05], 'String', {strrep(File_name, '_', '\_')}, 'Linestyle', 'none', 'HorizontalAlignment', 'center');

            ann1.FontSize = 16;

            ann1.FontWeight = 'bold';            
                
            ann2 = annotation('textbox', [0 0.97 1 .05], 'String', {Markup_tag}, 'FitBoxToText', 'on', 'Linestyle', 'none');
                        
            ann2.FontSize = 12;

            end
                
        orient landscape            
            
            try
        
                if i == 1 && j == 1                                

                print('Dummy', '-dpsc2', '-bestfit')

                else

                print('Dummy', '-dpsc2', '-append', '-bestfit')    

                end
            
            catch
                
                if i == 1 && j == 1                                

                print('Dummy', '-dpsc2', '-loose')

                else

                print('Dummy', '-dpsc2', '-append', '-loose')    

                end    
                
            end
    
        close(Subhandle)
    
        end
                     
    end   

    
%% Check the Ghostscript path leads to a exe file

    if ~Is_linux

    Ghostscript_fun_found = 0; % Initialise

    [Path_used, ~, ext] = fileparts(Ghostscript_path);

    Extension_used = 1;

        if isempty(ext)

        ext = '.exe'; 

        Extension_used = 0;

        end

        if ~Extension_used

        Ghostscript_path = [Ghostscript_path '.exe']; % If no extension is given check if file given is a exe anyway.

        end

        if strcmp(ext, '.exe') % If it is a exe file

            if isempty(Path_used) % Check this folder if no path is given

            Files_that_exist = Find_files_in_folders_GEARS(cd, '.exe', 0, 1);        

            else

            Files_that_exist = Find_files_in_folders_GEARS(Path_used, '.exe', 0, 1);    

            end

            if ismember(Ghostscript_path, Files_that_exist)

            Ghostscript_fun_found = 1;       

            end

        end

        if ~Ghostscript_fun_found

            if ispc

            warning('The ghostscript exe was not found, searching default location')    

            Ghostscript_path = Search_default_win_paths;    

                if ~isempty(Ghostscript_path)

                Ghostscript_fun_found = 1;

                display('Ghostscript found in the default location')                   

                end

            end

        end

        if ~Ghostscript_fun_found

        warning('The ghostscript exe was not found, the figure report will be left in postscript format.')  

        end

    else

    Ghostscript_path = 'gs';

    Ghostscript_fun_found = 1;

    end

    
%% Convert Postscript to pdf

    if Ghostscript_fun_found
        
        try    

        Ghost_script_options = '-dBATCH -dNOPAUSE -q';

        Command = sprintf('"%s" %s -sDEVICE=pdfwrite -dPDFFitPage -o %s %s', Ghostscript_path, Ghost_script_options, [Save_loc '.pdf'],  'Dummy.ps');   

        system(Command);

        delete('Dummy.ps')   

        catch
            
        movefile('Dummy.ps', [Save_loc '.ps'])    
            
        warning('There was a problem using Ghostscript, the figure report will be left in postscript format.')    
    
        end
    
    else 
        
    movefile('Dummy.ps', [Save_loc '.ps'])
                        
    end

end

function [Path_to_use] = Search_default_win_paths
% A function that searches default program paths for ghostscript if the provided path does not find a function.

Exe_files = [Find_files_in_folders_GEARS('C:\Program Files\gs', 'exe', 1, 0); Find_files_in_folders_GEARS('C:\Program Files (x86)\gs', 'exe', 1, 0)];

Files_wanted = {'gswin32c', 'gswin64c'};

File_names = cell(size(Exe_files));

    for i = 1:length(Exe_files)
    
    [~, File_names{i}] = fileparts(Exe_files{i});
        
    end

Options = Exe_files(ismember(File_names, Files_wanted));

    if ~isempty(Options)
        
        if length(Options) > 1
            
        Versions = cell(size(Options));
        
            for i = 1:length(Options)
                
            dummy = strsplit(Options{i}, filesep);
                
            Versions{i} = dummy{end - 2};    
                
            end
            
        Versions = sort(Versions); 
        
        Version_to_use = Versions(end); % Use the most recent version by default.
        
        Latest_versions = sort(Options(strcmp(Versions, Version_to_use)));
        
        Path_to_use = Latest_versions{1}; % Use 32 version by default
                                             
        else
            
        Path_to_use = Options{1};
                
        end        
        
    else
        
    Path_to_use = {};    
                
    end
     
end

