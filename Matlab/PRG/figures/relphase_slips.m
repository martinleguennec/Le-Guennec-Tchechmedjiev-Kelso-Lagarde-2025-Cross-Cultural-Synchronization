function relphases_slips = relphase_slips(relphases)
% relphase_slips : prepare relative phases data so that it can be plotted
%                  in an interval [-pi, pi]
%
% Calling Sequence
%   relphases_slips = relphase_slips(relphases)
%
% Parameters
%   relphases        : matrix,  of dimension Nx2, first column has relative
%                               phase data, second column is the
%                               corresponding time
%
% Output
%   relphases_slips  : matrix,  of dimension Mx2, same as relphases but
%                               with added NaN and extrem values
%
% Description
%   Detects when there is a phase slip and adds NaN with -pi or pi values 
%   before and after, so that when plotted the values before and after a 
%   phase slips are not linked

% Authors
%   XXX
%


%% Code
phase_wrap = unwrap(relphases(:,1));
t_phase = relphases(:,2);

% Phase slips correspond to change by +/- 2 pi, hence we round the phases
% to 2 pi then find when the phase increases or decreases by that amount
phase_round = round(phase_wrap ./ (2*pi));
disp(phase_round)
diff_phase = diff(phase_round);
phase_slip = find(abs(diff_phase)>0);
disp(phase_slip)

% Because we want the phases slips to be visible, we add a NaN when there
% is one so that when we plot it, the two dots around the NaN wont be
% linked together
% We also add pi or -pi before and after the NaN so that when we plot and
% we set ylim([-pi pi]), the line will go to the extremity of the plot
for i = 1:numel(phase_slip)
    iBeg = phase_slip(i)+(i-1)*3;    % Index before phase slip
    iEnd = phase_slip(i)+1+(i-1)*3;  % Index after phase slip

    % We must add or substract a little amount to be sure that when we
    % use mod afterwards it is on the right side
    if diff_phase(phase_slip(i)) > 0
        last_before_slip = pi - .001;
        first_after_slip = (-1)*pi + .001;
    elseif diff_phase(phase_slip(i)) < 0
        last_before_slip = (-1)*pi + .001;
        first_after_slip = pi - .001;
    end
    
    % Add the prepared values to the phases vector
    phase_wrap = [phase_wrap(1:iBeg); ...
        last_before_slip; NaN; first_after_slip;  ...
        phase_wrap(iEnd:end)];

    % Add values to the time vector too
    t_diff = t_phase(iEnd) - t_phase(iBeg);
    t_phase = [t_phase(1:iBeg); ...
        t_phase(iBeg) + t_diff/4; t_phase(iBeg) + t_diff/2; ...
        t_phase(iBeg) + 3*t_diff/4; ...
        t_phase(iEnd:end)];
end

% Put the phases values between -pi and pi
phase_wrap = mod(phase_wrap + pi, 2*pi) -pi;

% Output matrix
relphases_slips = [phase_wrap, t_phase];
