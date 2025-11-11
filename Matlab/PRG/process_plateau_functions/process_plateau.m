function process_plateau(plat_idx, num_plat, stim_struct, mov_phase_struct, infos, display_fig)
% PROCESS_PLATEAU Processes each plateau and stores the movement data between peaks.
%
% INPUTS:
%   plat_idx - Integer, index of the plateau that this function is 
%              processing.
%   num_plat - Integer, number of plateaus in the current trial.
%   stim_struct - Structure containing:
%       stim - Array, size S-by-1, time of the onset of the S stimuli.
%       stim_plat_indices - Array, size P-by-2, each row contains the index
%                           of the stimulus at the beginning and at the end 
%                           of the P-th plateau.
%   mov_phase_struct - Structure containing:
%       mov - Array, size N-by-2, first column for the N movement values,
%             second column for corresponding time.
%       peaks - Array, size R-by-2, first column for the peak value, second
%               column for the time of occurrence.
%       relphases - Array, size R-by-2, first column for the wrapped 
%                   relative phase values (bounded between -pi and +pi), 
%                   second column for the time of occurrence.
%       relphases_unwrap - Array, size R-by-2, first column for the 
%                          unbounded relative phase values, second column 
%                          for the time of occurrence.
%   infos - Structure containing subject information:
%       subject - String, the subject identifier.
%       task - String, the task identifier ('sync' or 'syncslow').
%       task_number - Integer, the task number.
%       group - String, group of the subject ('French' or 'Indian').
%   display_fig - Logical, true to display the figure; false otherwise.
%
% This function identifies the stimuli, movement, peaks, relative phase
% values (wrapped and unwrapped) for the current plateau. It then
% calculates the dwell time (an index of the synchronization coupling); and
% identifies the time spent in flexion, in extension, or dwelling for each
% movement.
%
% Version 1.1, Sep. 2024 by XXX, tested with Matlab R2021b 
% (9.11.0.1769968) running on a maci64 computer
% Contact: XXX

global tbl_relphase;

% Calculate plateau frequency and period
plat_freq = 1 + (plat_idx-1) * 0.3;
plat_period = 1 / plat_freq * 1000;

% Determine start and end of the plateau
[plat_beg, plat_end] = determine_plat_bounds(plat_idx, ...
    stim_struct.stim, stim_struct.stim_plat_indices, num_plat, plat_period);

% Identify data for the current plateau
plat_stim = stim_struct.stim(...
    stim_struct.stim_plat_indices(plat_idx, 1) : ... 
    stim_struct.stim_plat_indices(plat_idx, 2));
plat_movement = mov_phase_struct.mov(...
    (mov_phase_struct.mov(:,2) >= plat_beg) & ...
    (mov_phase_struct.mov(:,2) <= plat_end), :);
plat_peaks = mov_phase_struct.peaks(...
    (mov_phase_struct.peaks(:,2) >= plat_beg) & ...
    (mov_phase_struct.peaks(:,2) <= plat_end), :);
plat_relphase = mov_phase_struct.relphases(...
    (mov_phase_struct.relphases(:,2) >= plat_beg) & ...
    (mov_phase_struct.relphases(:,2) <= plat_end), :);
plat_relphase_unwrap = mov_phase_struct.relphases_unwrap(...
    (mov_phase_struct.relphases_unwrap(:,2) >= plat_beg) & ...
    (mov_phase_struct.relphases_unwrap(:,2) <= plat_end), :);

% Calculate dwell time from unwrapped relative phase values
plat_DT = dwell_time(plat_relphase_unwrap(4:end, 1), 0.2);


% Save results in tbl_relphase
% One line of the table corresponds to one tap
% Obviously, the first tap of each plateau won't have a period assigned to
% it since there is no previous tap to calculate it
size_tbl = size(plat_relphase, 1);
tbl_relphase_var_names = {'group', 'subject', 'task', 'task_number', ...
                          'frequency', 'tap', 'DT', ...
                          'period', 'relphase'};
tbl_relphase = [tbl_relphase; ...
    table(repmat(infos.group, size_tbl, 1), ...         % group
          repmat(infos.subject, size_tbl, 1), ...       % participant ID
          repmat(infos.task, size_tbl, 1), ...          % task
          repmat(infos.task_number, size_tbl, 1), ...   %trial
          repmat(plat_freq, size_tbl, 1), ...           % frequency of metronome
          (1:size_tbl)', ...                            % tap number within the plateau
          repmat(plat_DT, size_tbl, 1), ...             % dwell time
          [nan; diff(plat_relphase_unwrap(:, 2))], ...  % period between consecutive taps
          plat_relphase(:, 1), ...                      % relative phase
          'VariableNames', tbl_relphase_var_names)];

if display_fig
    plot_plateau(plat_idx, plat_stim, plat_movement, plat_peaks, ...
        plat_relphase, plat_relphase_unwrap, infos.subject, infos.task, ...
        infos.task_number, plat_beg, plat_end, plat_freq, plat_DT);
    pause, close
end
