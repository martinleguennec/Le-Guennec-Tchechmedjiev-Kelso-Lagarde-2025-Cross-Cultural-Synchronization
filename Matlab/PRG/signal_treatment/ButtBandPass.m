function filteredSignal = ButtBandPass(order, sampFreq, cutoffFreq, signal)
    % buttBandPass : band pass filter of a signal with Butterworh
    %
    % Calling Sequence
    %   filteredSignal = buttBandPass(order, samplingFreq, cutoffFreq, signal)
    %
    % Parameters
    %   order           : number,  the order of the filter
    %   sampFreq        : number,  sampling frequency in Hz
    %   cutoffFreq      : vector,  dimension 1 x 2, low and high cutoff 
    %                              frequency in Hz
    %   signal          : vector,  the signal to filter
    %
    % Output
    %   filteredSignal  : vector,  the signal filtered, the same size as
    %                              signal (input)
    %
    % Description
    %   buttBandPass : band pass filter of a signal with Butterworh. The
    %   Sf is a vector (same size as S)

    
    % Authors
    %   XXX



    
    % The signal must be a horizontal vector
    if size(signal, 1) < size(signal, 2)
        signal = signal';
    end
    
    nyquistWs = sampFreq/2; % We apply the Nyquist frequency
    [b, a] = butter(order, [(cutoffFreq(1) ./ nyquistWs) (cutoffFreq(2) ./ nyquistWs)], 'bandpass');
    
    % Duplicate the signal before and after it to delete the starting effect
    duplique = [rot90(rot90(signal)) ; signal ; rot90(rot90(signal))] ;
    
    % Apply filter on the duplicated signal
    filtrage = filtfilt(b, a, duplique) ; 
    
    % Cut the added values before and after signal
    coupe = filtrage(numel(signal)+1:2*numel(signal)) ;
    
    filteredSignal = coupe;
    
    clear duplique filtrage coupe;
end