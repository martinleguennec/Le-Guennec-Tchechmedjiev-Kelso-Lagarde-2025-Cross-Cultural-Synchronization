function [stim_times, stim_plat_indices, num_plat] = create_stim_indexes(stim_signal, samp_freq, subj, task)
% CREATE_STIM_INDEXES Returns an array with the time at which occur the
% stimuli, the stimuli indexes in each plateau and the number of plateaus
% in the trial
%
% INPUTS:
%   stim_signal - Array, size S-by-1, voltage of the signal sent to the
%                 speakers.
%   samp_freq - Double, the sampling frequency of the signal in Hz.
%   subject - String, the subject identifier.
%   task - String, the task identifier ('sync' or 'syncslow').
%
% OUTPUTS:
%   stim_times - Array, size S-by-1, time of the onset of the S stimuli.
%   stim_plat_indices - Array, size P-by-2, each row contains the index of 
%                       the stimulus at the beginning and at the end of the
%                       P-th plateau.
%   num_plat - Double, number of plateaus of frequency performed by the 
%              subject (depends on the task).
%
% This function identifies the onset of the stimuli, then creates an array
% containing the indices of the stimulus at the beginnind and ending of
% each plateau with the number of plateaus varying depending on the task.
%
% Version 1.0, aug. 2024 by XXX, tested with Matlab R2021b 
% (9.11.0.1769968) running on a maci64 computer
% Contact: XXX


% Identify stimuli indexes
stim_logical = find(abs(movmean(stim_signal, 3)) > 500);
diff_logical = find(diff(stim_logical) > 100);
diff_logical = [0; diff_logical] + 1;
stim_onset = stim_logical(diff_logical);

% The first two triggers and the last triggers aren't sound signals
stim_onset = stim_onset(3:end);

% Transform indexes into time
stim_times = stim_onset ./ samp_freq;

% Determine the number of plateaus and the total number of stimuli
if task == "sync"
    num_plat = 18;
    num_stim = 276;
    if numel(stim_times) < 276
        error("Not the right number of plateaus")
    end
elseif task == "syncslow"
    num_plat = 11;
    num_stim = 171;
    if numel(stim_times) < 171
        error("Not the right number of plateaus")
    end
end

% If not the right number of stimuli, delete the first detected when the
% period between them is not right (is_right_first_stim) or if it is 
% variable (is_not_variable), then delete the last in excess
if numel(stim_times) ~= num_stim

    diff_stim = diff(stim_times);  
    is_right_first_stim = false;
    i_beg = 1;

    while is_right_first_stim == false

        is_right_period = diff_stim(1) < 0.9;
        is_not_variable = abs(diff_stim(1) - diff_stim(2)) > 0.002;

        % Handle particular cases
        if subj == "IN04" || subj == "IN06"
            is_right_period = false;
            is_not_variable = false;
        end

        % Itteratively increase i_beg if the stimuli are not wanted or exit 
        % the loop
        if (is_right_period || is_not_variable)
            diff_stim = diff_stim(2:end);
            i_beg = i_beg + 1;  
        else
            is_right_first_stim = true;
        end
    end

    stim_times = stim_times(i_beg : i_beg + num_stim - 1);
end

% Prepare matrix with first and last stim indexes for each plateau
stim_plat_indices = zeros(num_plat, 2);
stim_plat_indices(1:2, :) = [1, 21 ; 22, 36];
for plateau_idx = 3:num_plat
    stim_plat_indices(plateau_idx, :) = stim_plat_indices(plateau_idx - 1, :) + 15;
end

