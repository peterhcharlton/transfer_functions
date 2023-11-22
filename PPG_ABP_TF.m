function transformed_sig = PPG_ABP_TF(S, sig_type, options)
% PPG_ABP_TF    Transforms digital PPG signals into either digital or radial
%               ABP signals (and vice-versa) using an empirical transfer
%               function.
%   
%  Inputs:
%
%    S         -  a pulsatile signal, consisting of either a signal pulse, or a
%                   signal containing several pulses. S should be a structure, containing
%                   a vector of amplitudes, S.v, and the sampling frequency (in Hz), S.fs.
%    sig_type  -  a string containing the input signal type (either 'digPPG',
%                   'digABP', or 'radABP'
%    options   -  (optional) a structure of options, which may contain any of:
%                    options.do_plot                    - a logical (true or false)
%                    options.retain_onset_time          - a logical
%
%  Outputs:
%
%    transformed_sig  -  the transformed pulsatile signal.
%
%  Exemplary usage:
%
%    transformed_sig = PPG_ABP_TF(S, 'digPPG')             transforms the digital PPG signal, S, into digital and radial BP signals
%    transformed_sig = PPG_ABP_TF(S, 'digABP')              transforms the digital ABP signal, S, into a digital PPG signal
%    transformed_sig = PulseAnalyse(___, options)          uses options to specify the analysis.
%
%  For further information please see the accompanying manual.
%
%  This script contains items either copied or modified from the RRest
%  toolbox which is covered by the GNU public licence (<a href="http://github.com/peterhcharlton/RRest/">link</a>).
%  It also contains items either copied or modified from PulseAnalyse.m .
%
% Peter H. Charlton, King's College London, August 2017

% orientate signal vector as a column
S.v = S.v(:);

% setup univeral parameters (for use throughout the script)
if nargin < 3
    options.do_plot = false;
    options.retain_onset_time = false;
end
up = setup_up(S, sig_type, options);

% pre-process input signal
[S_proc, up] = pre_process_input_signal(S, up);

% find FFT of input signal
S_proc_fft = fft_input_signal(S_proc, up);

% define transfer function
bode_plot_data = define_transfer_function(S_proc_fft, up);

% transform signal
transformed_sig = transform_signal(S_proc, S_proc_fft, bode_plot_data, up);

% post-process output signal
if strcmp(up.no_of_pulses, 'single')
    transformed_sig = post_process_output_signal(transformed_sig, up, options, S_proc.time_till_onset);
else
    transformed_sig = post_process_output_signal(transformed_sig, up, options);
end

end

function up = setup_up(S, sig_type, options)

%% Analysis settings

% Threshold signal duration to distinguish between a single pulse or multiple pulses: 
up.analysis.max_duration_of_single_pulse = 2.5;   % in secs

% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.157;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.02;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.Fpass = 38.5;  % in HZ
up.paramSet.elim_vhf.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
up.paramSet.elim_vhf.Dpass = 0.05;
up.paramSet.elim_vhf.Dstop = 0.01;

% Number of points for FFT
min_no_secs = 4;
actual_no_secs = (length(S.v)-1)/S.fs;
no_secs = max(min_no_secs, actual_no_secs);
up.n_fft = 2^nextpow2(no_secs*S.fs);

% Number of times to repeat a single pulse for analysis
up.no_beats = 3;

%% Determine analysis to be conducted

% Determine whether this is a single pulse, or several pulses
up.no_of_pulses = determine_no_of_pulses(S, up);

% Determine the transformations to be performed
switch sig_type
    case 'digPPG'
        up.transformations = {'digPPG_to_digABP', 'digPPG_to_radABP'};
    case 'digABP'
        up.transformations = {'digABP_to_digPPG'};   % , 'digPPG_to_radABP'
    case 'radABP'
        up.transformations = {'radABP_to_digPPG'};   % , 'digPPG_to_digABP'
end
up.input_sig = sig_type;

end

function no_of_pulses = determine_no_of_pulses(S, up)
% DETERMINE_NO_OF_PULSES  Determines whether this is a single pulse wave,
% of a pulsatile signal containing multiple pulse waves.

signal_duration = (length(S.v)-1)/S.fs;

% If duration of signal is greater than a threshold then assume this is
% multiple pulse waves:
if signal_duration > up.analysis.max_duration_of_single_pulse
    no_of_pulses = 'multiple';
else
    no_of_pulses = 'single';
end

end

function [S_proc, up] = pre_process_input_signal(S, up)

% Pre-processing is peformed according to whether this is a single pulse or multiple pulses.

switch up.no_of_pulses
    
    case 'multiple'
        
        % Eliminate very low frequency content
        %S_elim_lf = elim_vlfs(S, up);
        S_elim_lf = S;
        
        % Eliminate very high frequency content
        S_proc = elim_vhfs(S_elim_lf, up);
        
    case 'single'
        
        % Eliminate low frequency content
        S_elim_lf = eliminate_low_freq_from_single_beat(S, up);
        
        % Ensure that signal commences at start of systole
        [S_proc, S_proc.time_till_onset] = align_pulse(S_elim_lf, up);
        
        % Repeat pulse
        if S_proc.v(end) == S_proc.v(1)
            S_proc.v = S_proc.v(1:end-1);
        end
        up.no_samples_in_one_beat = length(S_proc.v);
        S_proc.v = repmat(S_proc.v, [up.no_beats,1]);

end

% % Normalise PPG
% ppg_norm = ppg.v - linspace(ppg.v(1), ppg.v(end), length(ppg.v));
% ppg_norm = ppg_norm - mean(ppg_norm);
% 

end

function s_filt = elim_vlfs(s, up)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vlf.Fstop up.paramSet.elim_vlf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vlf.Dstop up.paramSet.elim_vlf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0266;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

s_filt.v = filtfilt(AMfilter.numerator, 1, s.v);
s_filt.v = s.v-s_filt.v;
s_filt.fs = s.fs;
end

function s_filt = elim_vhfs(s, up)
%% Filter signal to remove VHFs
% Adapted from RRest

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Check to see if sampling freq is at least twice the freq of interest
if (up.paramSet.elim_vhf.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([up.paramSet.elim_vhf.Fstop up.paramSet.elim_vhf.Fpass]/(s.fs/2), [1 0], [up.paramSet.elim_vhf.Dstop up.paramSet.elim_vhf.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.139;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

%% Remove VHFs
s_dt=detrend(s.v);
s_filt.v = filtfilt(AMfilter.numerator, 1, s_dt);
end

function S_elim_lf = eliminate_low_freq_from_single_beat(sig, up)

% Correct for low frequency baseline drift in a single beat

diff_1 = sig.v(2) - sig.v(1);
desired_val_end = sig.v(1) - diff_1;  % old setting
desired_val_end = sig.v(1);
correction_line = linspace(0, sig.v(end)-desired_val_end, length(sig.v));
S_elim_lf.v = sig.v - correction_line';
S_elim_lf.v = S_elim_lf.v - mean(S_elim_lf.v);
S_elim_lf.fs = sig.fs;

end

function [S_aligned, time_till_onset] = align_pulse(sig, up)

% Ensure that signal commences at start of systole

[~, align_el] = min(sig.v);
S_aligned.v = sig.v([align_el:end, 1:(align_el-1)]);
S_aligned.fs = sig.fs;
time_till_onset = (align_el-1)/sig.fs;

end

function S_proc_fft = fft_input_signal(S_proc, up)

% find FFT of input signal
S_proc_fft.raw = fft(S_proc.v, up.n_fft);
S_proc_fft.v = S_proc_fft.raw(1:up.n_fft/2+1);
S_proc_fft.f = S_proc.fs*(0:(up.n_fft/2))/up.n_fft; S_proc_fft.f = S_proc_fft.f(:);

end

function bode_plot_data = define_transfer_function(S_proc_fft, up)

% insert data from published Bode Plot
bode_data = insert_bode_data;

% find values of Bode data at these frequencies
bode_plot_data = calc_rel_bode_data(bode_data, S_proc_fft.f);

end

function bode_data = insert_bode_data

% insert data from Fig. 2 (b) of "Noninvasive assessment of the digital
% volume pulse. Comparison with the peripheral pressure pulse." (used the
% solid line)

%% "Tonometer / Photoplethysmograph" Amplitude
% x-coordinates in cm
x = [0, 0.9, 1.8, 3.5, 5.3, 7.1, 7.9];
% y-coordinates in cm
y = [1.5, 1.6, 2.25, 2.9, 3.3, 2.3, 1.9];
% x-coordinates in Hz
bode_data.digPPG_to_radABP.amp.f = x*10/8.8;    % 8.8cm is 10 Hz
% y-coordinates in scaling factor
bode_data.digPPG_to_radABP.amp.v = y*3/6.0;      % 6.0cm is 3
clear x y

%% "Tonometer / Photoplethysmograph" Phase
% x-coordinates in cm
x = [0,   0.95, 1.9, 3.6, 4.7,  5.5, 8.4];
% y-coordinates in cm
y = [3.8, 4.3,  4.2, 3.8, 3.65, 3.5, 3.25];
% x-coordinates in Hz
bode_data.digPPG_to_radABP.phase.f = x*10/9.2;    % 9.2cm is 10 Hz
% y-coordinates in deg
bode_data.digPPG_to_radABP.phase.v = (y-3.8)*135/3.4;      % The x-axis is 3.8cm above the measurement line, and 135 deg is 3.4 cm
clear x y

%% rad ABP to dig PPG
bode_data.radABP_to_digPPG = bode_data.digPPG_to_radABP;
bode_data.radABP_to_digPPG.amp.v = 1./bode_data.digPPG_to_radABP.amp.v;
bode_data.radABP_to_digPPG.phase.v = -1*bode_data.digPPG_to_radABP.phase.v;

%% "Finapress / Photoplethysmograph" Amplitude
% x-coordinates in cm
x = [0.0,  0.85, 1.75, 2.6, 3.6,  5.25, 6.15, 7.0, 7.85];
% y-coordinates in cm
y = [1.45, 1.5,  2.1,  2.7, 3.25, 4.3,  4.0,  2.9, 2.15];
% x-coordinates in Hz
bode_data.digPPG_to_digABP.amp.f = x*10/8.75;    % 8.75cm is 10 Hz
% y-coordinates in scaling factor
bode_data.digPPG_to_digABP.amp.v = y*3/5.9;      % 5.9cm is 3
clear x y

%% "Finapress / Photoplethysmograph" Phase
% x-coordinates in cm
x = [0,   0.8, 2.15, 2.75, 3.65, 4.6, 6.4, 7.15, 8.25];
% y-coordinates in cm
y = [3.8, 4.3, 4.35, 4.2,  3.9,  3.8, 3.1, 2.95, 2.7];
% x-coordinates in Hz
bode_data.digPPG_to_digABP.phase.f = x*10/9.25;    % 9.25cm is 10 Hz
% y-coordinates in deg
bode_data.digPPG_to_digABP.phase.v = (y-3.8)*135/3.4;      % The x-axis is 3.8cm above the measurement line, and 135 deg is 3.4 cm
clear x y

%% dig ABP to dig PPG
bode_data.digABP_to_digPPG = bode_data.digPPG_to_digABP;
bode_data.digABP_to_digPPG.amp.v = 1./bode_data.digPPG_to_digABP.amp.v;
bode_data.digABP_to_digPPG.phase.v = -1*bode_data.digPPG_to_digABP.phase.v;

end

function transformed_sig = transform_signal(S_proc, S_proc_fft, bode_plot_data, up)

% store processed input signal
eval(['transformed_sig.' up.input_sig ' = S_proc;']);

% make input signal into a variable
eval([up.input_sig ' = S_proc;']);
eval([up.input_sig '.fft = S_proc_fft;']);

% carry out each transformation in turn
for transformation_no = 1 : length(up.transformations)
    
    output_sig_name = up.transformations{transformation_no}(length(up.input_sig)+5:end);
    eval(['rel_bode_data = bode_plot_data.' up.transformations{transformation_no} ';'])
    eval(['rel_sigs.input = ' up.transformations{transformation_no}(1:6) ';']);
    
    % Recover output signal from FFT of input signal
    N = up.n_fft;
    x = nan(N,1);
    x_bar = calc_x_bar(S_proc_fft);
    rel_sigs.output.v = nan(N-1,1);
    rel_sigs.output.fs = rel_sigs.input.fs;
    temp.real_x_bar = real(x_bar(1:(N/2)+1));
    temp.imag_x_bar = imag(x_bar(1:(N/2)+1));
    temp.bode_phase = (rel_bode_data.phase.v*2*pi/360);
    for i = 1 : N-1
        
        % Find sum of cosines
%         sum_cosines = 0;
%         for k = 0 : N/2
%             amplitude_factor = rel_bode_data.amp.v(k+1);
%             phase_offset = rel_bode_data.phase.v(k+1)*2*pi/360;  % since original is in degrees
%             temp_to_add = real(x_bar(k+1)) * amplitude_factor * cos((2*pi*k*i/N) + phase_offset);
%             sum_cosines = sum_cosines + temp_to_add;
%             clear temp_to_add
%         end
%         clear k
        
        % sum_cosines = sum( real(x_bar(1:(N/2)+1)) .* rel_bode_data.amp.v .* cos((2*pi*([0 : N/2]')*i/N) + (rel_bode_data.phase.v*2*pi/360)) );
        
        sum_cosines = sum( temp.real_x_bar .* rel_bode_data.amp.v .* cos((2*pi*([0 : N/2]')*i/N) + temp.bode_phase) );
        
        % Find sum of sines
%         sum_sines = 0;
%         for k = 0 : N/2
%             amplitude_factor = rel_bode_data.amp.v(k+1);
%             phase_offset = rel_bode_data.phase.v(k+1)*2*pi/360;  % since original is in degrees
%             temp_to_add = imag(x_bar(k+1)) * amplitude_factor * sin((2*pi*k*i/N) + phase_offset);
%             sum_sines = sum_sines + temp_to_add;
%             clear temp_to_add
%         end
%         clear k
        
        % sum_sines = sum( imag(x_bar(1:(N/2)+1)) .* rel_bode_data.amp.v .* sin((2*pi*([0 : N/2]')*i/N) + (rel_bode_data.phase.v*2*pi/360)) );
        
        sum_sines = sum( temp.imag_x_bar .* rel_bode_data.amp.v .* sin((2*pi*([0 : N/2]')*i/N) + temp.bode_phase) );
        
        % add sines and cosines
        rel_sigs.output.v(i) = sum_cosines + sum_sines;
        clear sum_sines sum_cosines
        
    end
    clear i
    
    % store sampling frequency
    rel_sigs.output.fs = rel_sigs.input.fs;

    % Store this transformed signal
    eval([output_sig_name ' = rel_sigs.output;']);
    eval(['transformed_sig.' output_sig_name ' = ' output_sig_name ';']);
    clear rel_sigs
    
end


end

function x_bar = calc_x_bar(fft)

% using methodology at: http://www.dspguide.com/ch8/5.htm

x = fft.raw;
N = length(x);
real_x_bar = nan(size(x));
imag_x_bar = nan(size(x));
for s = 1 : N
    
    if s == 1 || s == N/2
        real_x_bar(s) =    real(x(s)) / N;
        imag_x_bar(s) = -1*imag(x(s)) / N;
    else
        real_x_bar(s) =    real(x(s)) / (N/2);
        imag_x_bar(s) = -1*imag(x(s)) / (N/2);
    end
    
end

x_bar = real_x_bar + 1i*imag_x_bar;

end

function rel_bode_data = calc_rel_bode_data(bode_data, req_freqs)

% evaluate the bode data at the required frequencies

% cycle through each of the bode plots
bode_plot_names = fieldnames(bode_data);

for bode_plot_no = 1 : length(bode_plot_names)
    
    % extract data for this bode plot
    eval(['curr_bode_data = bode_data.' bode_plot_names{bode_plot_no} ';'])
    
    % amplitudes
    curr_rel_bode_data.amp.f = req_freqs;
    req_freqs_within_rel_range = req_freqs(req_freqs<=curr_bode_data.amp.f(end));
    %curr_rel_bode_data.amp.v = zeros(size(req_freqs));
    curr_rel_bode_data.amp.v = curr_bode_data.amp.v(end)*ones(size(req_freqs));
    curr_rel_bode_data.amp.v(1:length(req_freqs_within_rel_range)) = interp1(curr_bode_data.amp.f, curr_bode_data.amp.v, req_freqs_within_rel_range);
    
    % phase
    curr_rel_bode_data.phase.f = req_freqs;
    req_freqs_within_rel_range = req_freqs(req_freqs<=curr_bode_data.phase.f(end));
    % curr_rel_bode_data.phase.v = zeros(size(req_freqs));
    curr_rel_bode_data.phase.v = curr_bode_data.phase.v(end)*ones(size(req_freqs));
    curr_rel_bode_data.phase.v(1:length(req_freqs_within_rel_range)) = interp1(curr_bode_data.phase.f, curr_bode_data.phase.v, req_freqs_within_rel_range);
    
    % re-insert data for this bode plot
    eval(['rel_bode_data.' bode_plot_names{bode_plot_no} ' = curr_rel_bode_data;'])    
    
end

end

function transformed_sig = post_process_output_signal(transformed_sig, up, options, time_till_onset)

% Post-processing is peformed according to whether this is a single pulse or multiple pulses.

switch up.no_of_pulses
    
    case 'multiple'
        
        % Do nothing
        
    case 'single'
        
        %% Extract one pulse from the multiple pulses used in the analysis
        
        % Determine which samples to extract
        samples_in_one_pulse = up.no_samples_in_one_beat;
        rel_pulse_no = ceil(up.no_beats/2);
        rel_pulse_start_el = ((rel_pulse_no-1)*samples_in_one_pulse) + 1;
        rel_pulse_end_el = (rel_pulse_no*samples_in_one_pulse);
        
        % extract these samples
        sig_names = fieldnames(transformed_sig);
        for sig_no = 1 : length(sig_names)
            eval(['curr_sig = transformed_sig.' sig_names{sig_no} ';'])            
            curr_sig.v = curr_sig.v(rel_pulse_start_el : rel_pulse_end_el);
            % align beat
            curr_sig = align_pulse(curr_sig, up);
            eval(['transformed_sig.' sig_names{sig_no} ' = curr_sig;'])
        end
        
        %% Scale to occupy range of 0 to 1
        for sig_no = 1 : length(sig_names)
            eval(['curr_sig = transformed_sig.' sig_names{sig_no} ';'])            
            curr_sig.v = curr_sig.v-min(curr_sig.v);
            curr_sig.v = curr_sig.v / max(curr_sig.v);
            eval(['transformed_sig.' sig_names{sig_no} ' = curr_sig;'])
        end
        
        %% Shift to retain original onset time (if specified in options)
        temp = fieldnames(options);
        if sum(strcmp(temp, 'retain_onset_time')) && options.retain_onset_time && strcmp(up.no_of_pulses, 'single')
            for sig_no = 1 : length(sig_names)
                eval(['curr_sig = transformed_sig.' sig_names{sig_no} ';'])
                no_samples_to_shift = round(curr_sig.fs * time_till_onset);
                curr_sig.v = [curr_sig.v(end-no_samples_to_shift+1 : end); curr_sig.v(1: end-no_samples_to_shift)];
                eval(['transformed_sig.' sig_names{sig_no} ' = curr_sig;'])
            end
        end

end

if options.do_plot
    ftsize = 14;
    a = transformed_sig;
    sigs = fieldnames(a);
    for sig_no = 1 : length(sigs)
        eval(['curr_sig = a.' sigs{sig_no} ';'])
        plot([0:length(curr_sig.v)-1]/curr_sig.fs, curr_sig.v), hold on,
    end
    legend(sigs(:)')
    set(gca, 'YTick', [], 'FontSize', ftsize)
    xlabel('Time [s]', 'FontSize', ftsize)
    ylab = ylabel({'Signal', '[au]'}, 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'position', get(ylab,'position')-[0.06,0,0]);
    ylim([0 1.1])
    box off
end

end

function extras

%     
%     % Recover PPG from FFT
%     N = n_fft;
%     x = nan(N,1);
%     x_bar = calc_x_bar(ppg_fft);
%     for i = 1 : N-1
%         
%         % Find sum of cosines
%         sum_cosines = 0;
%         for k = 0 : N/2
%             if ppg_fft.f(k+1) < 10
%                 scale_factor = 1;
%             else
%                 scale_factor = 0;
%             end
%             
%             temp_to_add = real(x_bar(k+1)) * cos(2*pi*k*i/N) * scale_factor;
%             sum_cosines = sum_cosines + temp_to_add;
%             clear temp_to_add
%         end
%         clear k
%         
%         % Find sum of sines
%         sum_sines = 0;
%         for k = 0 : N/2
%             if ppg_fft.f(k+1) < 10
%                 scale_factor = 1;
%             else
%                 scale_factor = 0;
%             end
%             
%             temp_to_add = imag(x_bar(k+1)) * sin(2*pi*k*i/N) * scale_factor;
%             sum_sines = sum_sines + temp_to_add;
%             clear temp_to_add
%         end
%         clear k
%         
%         % add sines and cosines
%         x(i) = sum_cosines + sum_sines;
%         clear sum_sines sum_cosines
%         
%     end
%     clear i
%     
%     % find FFT of ABP
%     abp_fft = ppg_fft;
%     
%     % find ABP
%     abp = ifft(abp_fft);


end