function setup()
% Adds paths of HBF to Matlab.

%  Copyright (c) 2017 Haizhao Yang
%  National University of Singapore 
%  This file is distributed under the terms of the MIT License.

file_path = mfilename('fullpath');
tmp = strfind(file_path,'setup');
file_path = file_path(1:(tmp(end)-1));

% Foulder for all soource files recursively
addpath(genpath([file_path 'inv']));
addpath(genpath([file_path 'src']));
addpath(genpath([file_path 'kernels']));

end
