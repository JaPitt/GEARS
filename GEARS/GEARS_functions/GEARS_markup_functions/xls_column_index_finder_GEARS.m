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


function [index] = xls_column_index_finder_GEARS(index)
% A function that transforms numbers into excel column titles and vice versa.
% index - If scaler will be transformed into the index ith excel column title. If char it will transform the title into the column number. (scaler or char)


%% Check inputs

    if ~ischar(index) && ~isnumeric(index)

    error('index should be either a string or a number')
        
    end
    
    
%% letter -> numbers    

    if ischar(index)

    index = upper(index);

    Alphabet = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U' 'V', 'W', 'X', 'Y', 'Z'};

    Number_letters = length(index); % Number of letters in the input

    Letter_pos = zeros(1, length(Number_letters));

        for i = 1:Number_letters

        Letter_pos(i) = find(strcmp(Alphabet, index(i))); % Find the pos of each letter     

        end

    index = sum((26*ones(size(Letter_pos))).^(Number_letters-1:-1:0).*Letter_pos); % Calculate the number index

    else
    
    
%% numbers -> letters

        if isnumeric(index)

        Letter_index = [];

        Alphabet = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U' 'V', 'W', 'X', 'Y', 'Z'};

        Keep_going = 1; % Initialise

            while Keep_going

                if index < 27 % If we have reached the last letter to be added

                Keep_going = 0; 

                end    

            Last_letter = mod(index, 26); % The position of the current letter being considered

                if Last_letter == 0 % Using mod give 26 as zero.

                Last_letter = 26; % Correct for z

                end

            Letter_index = [Alphabet(Last_letter) Letter_index]; % Extend the index found

            index = (index - Last_letter)/26; % Reduce the index to account for the letter index just found

            end

        index = strrep(strjoin(Letter_index), ' ', ''); % Join the letters together and elimate any spaces   

        end
    
    end

end

