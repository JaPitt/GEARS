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


function[Files] = Find_files_in_folders_GEARS(Folder_name, Extension, Include_subfolders, Case_sensitive_extension)
% A function that creates a list of all the files within a given folder and its subfolders. Can additionally be used to single out a specific type of file.
% Folder_name              - The folder you want to search in (string).
% Extension                - The extension to be searched for. Leave empty not to search. (string or []). Optional: default [].
% Include_subfolders       - If the content of subfolders of Folder_name should be included in the search. (logical) Optional: default 0.
% Case_sensitive_extension - If the extension search should be case sensitive or not (logical value) Optional: default 0.


%% Check inputs

Ins = who;

    if ~ismember('Folder_name', Ins)
        
    error('A Folder name is required')
        
    end
    
    if ~ismember('Extension', Ins)
        
    Extension_search = 0;            
        
    else
        
        if isempty(Extension)
            
        Extension_search = 0;   
            
        else
            
        Extension_search = 1;    
        
        end
   
    end
    
    if ~ismember('Include_subfolders', Ins)
        
    Include_subfolders = 0;            
        
    else
        
        if isempty(Include_subfolders)
            
        Include_subfolders = 0;               
        
        end
   
    end
    
    if ~ismember('Case_sensitive_extension', Ins)
        
    Case_sensitive_extension = 0;            
        
    else
        
        if isempty(Case_sensitive_extension)
            
        Case_sensitive_extension = 0;               
        
        end
   
    end
    
    if isempty(dir(Folder_name))
        
    warning('The requested search folder was not found')
    
    Files = {};
    
    return
        
    end
   

%% Find files

Start_time = datetime;
	
	if Include_subfolders

	Files = getAllFiles(Folder_name); % Find all files in both folders and subfolders
	
	else
	
	Dir_info = dir(Folder_name);      % Get the data for the current directory
  
	Is_dir = [Dir_info.isdir];  % Find the index for directories
  
	Files = {Dir_info(~Is_dir).name}';  % All files in Folder_name (no subfolders)
    
        if ~isempty(Files)
      
        Files = cellfun(@(x) fullfile(Folder_name, x), Files, 'UniformOutput', false);

        end
		
    end
  
    
%% Find files with specific extensions

    if Extension_search % Decide if we are looking for specific extension
        
    Num_pot_files = length(Files); % To see the number of total files.
        
    Extension = strrep(Extension, '.', ''); % Ensure inputs can use . or not
    
    Files = Files(~cellfun(@isempty, strfind(upper(Files), ['.' upper(Extension)]))); % Reduce number of files that need to checked more seriously. To avoid confusion with things such as .m and .mat more is needed.

    Extension_finder = @(x) strsplit(x, '.');
    
    End_fun = @(x) x(end);
        
    Extensions_found = cellfun(@(x) End_fun(Extension_finder(x)), Files);
    
        if Case_sensitive_extension % Find files with wanted extension
              
        Files = Files(strcmp(Extensions_found, Extension)); % Case sensitive  
            
        else 
            
        Files = Files(strcmpi(Extensions_found, Extension)); % Case insensitive    
            
        end
        

%% Outputs

        if isempty(Files) 

        disp(['No files with the extension ' Extension ' were found.']);        

        else
            
         Time_taken = seconds(datetime - Start_time);

        disp([num2str(length(Files)) ' out of a potential ' num2str(Num_pot_files) ' files were found with the extension ' Extension ' in ' num2str(Time_taken) ' seconds.']);             
            
        end

    else 
        
        if isempty(Files) 

        disp('This folder is empty.');        
          
        else 
            
        Time_taken = seconds(datetime - Start_time);    

        disp([num2str(length(Files)) ' files were found in ' num2str(Time_taken) ' seconds.']);
    
        end
     
    end    
       
end

function File_list = getAllFiles(Dir_name)
% A function that finds all files in the folder and its subfolders.
% Dir_name - The name of the dir that will be searched for files

Dir_info = dir(Dir_name);      % Get the data for the current directory
  
Is_dir = [Dir_info.isdir];  % Find the index for directories
  
File_list = {Dir_info(~Is_dir).name}';  % Get a list of the files
  
	if ~isempty(File_list)
      
	File_list = cellfun(@(x) fullfile(Dir_name, x), File_list,'UniformOutput', false);

	end
  
Sub_dirs = {Dir_info(Is_dir).name};  % Get a list of the subdirectories
  
Valid_index = ~ismember(Sub_dirs,{'.','..'});  % Find index of subdirectories
                                               %  that are not '.' or '..'
	for i = find(Valid_index)                  % Loop over valid subdirectories
      
	Next_dir = fullfile(Dir_name,Sub_dirs{i});    % Get the subdirectory path
    
	File_list = [File_list; getAllFiles(Next_dir)];  % Recursively call getAllFiles
    
    end

end

