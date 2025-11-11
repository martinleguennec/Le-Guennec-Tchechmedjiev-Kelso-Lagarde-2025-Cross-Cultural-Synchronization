function [plat_beg, plat_end] = determine_plat_bounds(i_plat, stim, stim_plat_indices, num_plat, plat_period)
% DETERMINE_PLATEAU_BOUNDS Determines the start and end of each plateau.
%
% INPUTS:
%   i_plat - Integer, current plateau index.
%   stim - Array, size S-by-1, time of the onset of the S stimuli.
%   stim_plat_indices - Array, size P-by-2, each row contains the index of 
%                       the stimulus at the beginning and at the end of the
%                       P-th plateau.
%   num_plat - Integer, number of plateaus in the current trial.
%   plat_period - Double, period of the plateau in milliseconds.
%
% OUTPUTS:
%   plat_beg - Double, start time of the plateau.
%   plat_end - Double, end time of the plateau.
%
% The function determines the start and end of each plateau based on the
% stimulli and the indices of the plateaus.
%
% Version 1.0, aug. 2024 by XXX, tested with Matlab R2021b 
% (9.11.0.1769968) running on a maci64 computer
% Contact: XXX

if i_plat == 1
    plat_beg = stim(stim_plat_indices(i_plat, 1)) - 0.5;
    plat_end = mean([stim(stim_plat_indices(i_plat, 2)), stim(stim_plat_indices(i_plat + 1, 1))]);
elseif i_plat == num_plat
    plat_beg = mean([stim(stim_plat_indices(i_plat - 1, 2)), stim(stim_plat_indices(i_plat, 1))]);
    plat_end = stim(stim_plat_indices(i_plat, 2)) + 0.5 * plat_period / 1000;
else
    plat_beg = mean([stim(stim_plat_indices(i_plat - 1, 2)), stim(stim_plat_indices(i_plat, 1))]);
    plat_end = mean([stim(stim_plat_indices(i_plat, 2)), stim(stim_plat_indices(i_plat + 1, 1))]);
end
