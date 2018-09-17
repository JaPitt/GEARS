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


function [] = Add_GEARS_paths()
% Adds the GEARS package to the Matlab path. This will add all paths required to run the GEARS package.

addpath(cd) % This folder.

addpath(genpath([cd filesep 'GEARS_functions']))

addpath(genpath([cd filesep 'GEARS_iterating_functions']))

addpath(genpath([cd filesep 'GEARS_format_for_external_packages']))

addpath(genpath([cd filesep 'Third_party_tools']))

disp('The GEARS package has been added to the Matlab path. Enjoy!')

end

