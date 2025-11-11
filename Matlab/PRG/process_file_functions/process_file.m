function process_file(files, file_idx, subj, group)
% PROCESS_FILE Processes files within a subject's directory.
%
% INPUTS:
%   files - Structure array, information about the files in the directory.
%   file_idx - Integer, index of the file to be processed.
%   subj - String, subject identifier.
%   group - String, group of the subject ('French' or 'Indian').
%
% This function loops through each file in the subject's directory,
% processes the data, and stores it in the global structure.
%
% Version 1.1, Sep. 2024 by XXX, tested with Matlab R2021b
% (9.11.0.1769968) running on a maci64 computer
% Contact: XXX

global SAMP_FREQ_STIM SAMP_FREQ_MOV;

% Save the current file path, test name, and test number
path_file = fullfile(files(file_idx).folder, files(file_idx).name);
file_name_parts = split(files(file_idx).name, ["_", "."]);
task = convertCharsToStrings(file_name_parts{2});
task_number = str2double(file_name_parts{3});

[stim_signal, mov] = read_file(path_file, SAMP_FREQ_MOV);

[stim, stim_plat_indices, num_plat] = ...
    create_stim_indexes(stim_signal, SAMP_FREQ_STIM, subj, task);

[peaks, relphases, relphases_unwrap] = ...
    compute_relphase(mov, subj, task, task_number, ...
                     stim_plat_indices, stim, SAMP_FREQ_MOV);

% Create structure for stimulus information
stim_struct.stim = stim;
stim_struct.stim_plat_indices = stim_plat_indices;

% Create structure for movement and phase information
mov_phase_struct.mov = mov;
mov_phase_struct.peaks = peaks;
mov_phase_struct.relphases = relphases;
mov_phase_struct.relphases_unwrap = relphases_unwrap;

% Create structure for the process_plateau function
infos.subject = subj;
infos.task = task;
infos.task_number = task_number;
infos.group = group;

for plat_idx = 1:num_plat
    process_plateau(plat_idx, num_plat, stim_struct, ...
                    mov_phase_struct, infos, false);
end

