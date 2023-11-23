function transformed_sig = per_cen_TF(S, sig_type, options)
% PER_CEN_TF    Transforms peripheral ABP signals into central signals (and
% vice-versa) using an empirical transfer function.
%   
%  # Inputs:
%
%    S         -  a pulsatile signal, consisting of either a signal pulse, or a
%                   signal containing several pulses. S should be a structure, containing
%                   a vector of amplitudes, S.v, and the sampling frequency (in Hz), S.fs.
%    sig_type  -  a string containing the input signal type (either
%                   'brachBP')
%    options   -  (optional) a structure of options, which may contain any of:
%                    options.do_plot                    - a logical (true or false)
%                    options.retain_onset_time          - a logical
%
%  # Outputs:
%
%    transformed_sig  -  the transformed pulsatile signal.
%
%  # Exemplary usage:
%
%    transformed_sig = per_cen_TF(S, 'brachABP')             transforms the brachial arterial blood pressure (ABP) signal, S, into a central BP signal
%    transformed_sig = per_cen_TF(S, 'radiaABP')             transforms the radial arterial blood pressure (ABP) signal, S, into a central BP signal
%    transformed_sig = per_cen_TF(S, 'carotABP')             transforms the carotid arterial blood pressure (ABP) signal, S, into a central BP signal
%    transformed_sig = per_cen_TF                            runs example
%
%   # Author
%   Peter H. Charlton
%   
%   # Documentation
%   <https://github.com/peterhcharlton/transfer_functions/>
%   
%   # License - MIT
%      Copyright (c) 2023 Peter H. Charlton
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% add signal vector and example options if required
if nargin < 2
    sig_type = 'brachABP';
end
if nargin < 3
    options.do_plot = false;
    options.retain_onset_time = false;
end
if nargin == 0
    S.v = [71.3192829820745;71.8659405509528;74.2534441409816;79.1625107700863;86.3067832793381;95.2293502459513;104.652179719127;112.451218597538;117.723444305620;120.914367981151;122.654755393275;123.356366836000;123.238957803213;122.384280338716;120.817835080023;118.611623827181;115.923150517787;112.957920289471;109.910187998602;106.926207301524;104.093730392040;101.459258432680;99.0415418922058;96.8435808456187;94.8571250866637;93.0604240715774;91.4017265065773;89.7880300666053;88.0805812828248;86.1586271364101;84.0949194461569;82.2989684498210;81.2192853601828;81.0276244710996;81.7089854075620;83.0158619819154;84.4599919876046;85.6468655618829;86.4759801921877;87.0373381285752;87.4516923898707;87.7895447386185;88.0651455310773;88.2552441859827;88.3238398033121;88.2491818393021;88.0425205752095;87.7391068923064;87.3749416906154;86.9762756264029;86.5558590184267;86.1121918291860;85.6347737961742;85.1108546006333;84.5314340175578;83.8942619906961;83.2083387450540;82.4886646556408;81.7547402099687;81.0275659330508;80.3213921811460;79.6452191792600;79.0027970211450;78.3933756880507;77.8139551049752;77.2645352719184;76.7443661701299;76.2549478371106;75.7962802728605;75.3668634398787;74.9644472819138;74.5852817052134;74.2253916104000;73.8822269337221;73.5539126283035;73.2394736697686;72.9383100431169;72.6492217183478;72.3701836448350;72.0981957475766;71.8293579290702;71.5621701518148;71.3470336758419]; % data from 25-yr old baseline subject of PWDB
    S.v = [67.6921173029326;69.5339784614018;77.4237658180260;90.3545041984259;102.521723490848;110.823346143355;116.671657463078;121.334939157061;124.674437756466;126.290393267294;126.458812589718;125.755710124096;124.403841439320;122.523959554213;120.440822587729;118.459688171309;116.706559455031;115.264688520198;114.240578029376;113.636478038817;113.293384573420;113.034293208076;112.796952386496;112.585862221182;112.330520949591;111.855424184112;110.935566299605;109.306940695905;106.617038560292;102.588352955093;98.0541547120768;96.3617775145871;98.0752354629761;99.9139465426934;99.5821033585541;98.6764958303062;98.1908988023208;97.9055567807105;97.2609557776258;96.2488455868113;95.2314852647435;94.2876267802217;93.3287679206906;92.3676590049079;91.3975498641197;90.3636891295417;89.3103279074514;88.3664694229296;87.5786148385053;86.8755123728839;86.1664097572589;85.4325565228683;84.7039534197310;83.9926007478545;83.2782480009761;82.5488948790884;81.8232918509530;81.1231894603335;80.4545878572337;79.8099868541490;79.1796362073231;78.5567857480019;77.9391854199340;77.3260852043689;76.7212351950590;76.1223853357528;75.5295356264503;74.9412610315258;74.3592115922301;73.7842873310639;73.2161882405269;72.6540892999937;72.0980655113393;71.5493919064395;71.0093435171700;70.4782203510311;69.9554973948976;69.4404996318938;68.9325520451444;68.4310546196491;67.9384074154093]; % data for 75-year old from PWDB
    S.fs = 100; % Hz
    sig_type = 'radiaABP';
    options.do_plot = true;
end

% orientate signal vector as a column
S.v = S.v(:);

% add ABP to the end of the sig type if necessary
if isempty(strfind(sig_type, 'BP'))
    sig_type = [sig_type, 'ABP'];
elseif isempty(strfind(sig_type, 'ABP'))
    sig_type = strrep(sig_type, 'BP', 'ABP');
end
% adjust artery names
sig_type = strrep(sig_type, 'braABP', 'brachABP');
sig_type = strrep(sig_type, 'carABP', 'carotABP');
sig_type = strrep(sig_type, 'radABP', 'radiaABP');

% setup univeral parameters (for use throughout the script)
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

% Filter characteristics: Eliminate VLFs (below cardiac freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.6;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.2;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
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
up.n_fft = 2^nextpow2(no_secs*S.fs+1);

% Number of times to repeat a single pulse for analysis
up.no_beats = 3;

%% Determine analysis to be conducted

% Determine whether this is a single pulse, or several pulses
up.no_of_pulses = determine_no_of_pulses(S, up);

% Determine the transformations to be performed
switch sig_type
    case {'brachBP', 'brachABP'}
        up.transformations = {'brachABP_to_centrABP'};
    case {'carotBP', 'carotABP', 'carotidBP', 'carotidABP'}
        up.transformations = {'carotABP_to_centrABP'};
    case {'centrABP', 'centrBP'}
        up.transformations = {'centrABP_to_brachABP', 'centrABP_to_radiaABP'};
    case {'radiaBP', 'radiaABP', 'radialBP', 'radialABP', 'radBP', 'radABP'}
        up.transformations = {'radiaABP_to_centrABP'};
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

init.mean = mean(S.v);
init.sd = std(S.v);

% Pre-processing is peformed according to whether this is a single pulse or multiple pulses.

switch up.no_of_pulses
    
    case 'multiple'
        
        % Eliminate very low frequency content
%         S_elim_lf = elim_vlfs(S, up);
        
        % Eliminate very high frequency content
%         S_proc = elim_vhfs(S_elim_lf, up);
          S_proc = elim_vhfs(S, up);
        
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

S_proc.v = S_proc.v-mean(S_proc.v);
S_proc.v = init.sd*S_proc.v/std(S_proc.v);
S_proc.v = S_proc.v + init.mean;

end

function s_filt = elim_vlfs(s, up)
%% Filter pre-processed signal to remove frequencies below cardiac
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

% insert data from Fig. 2 (b) of:
% Karamanoglu, M.; O?Rourke, M. F.; Avolio, A. P.; Kelly, R. P. An analysis of the relationship between central aortic and peripheral upper limb pressure waves in man. Eur. Heart J. 1993, 14, 160?7.
% (used control data)

%% "Brachial / Ascending Aortic" Amplitude
% x-coordinates in Hz
bode_data.centrABP_to_brachABP.amp.f = 0:9;
% y-coordinates (rough guess)
y = [1, 1.05, 1.2, 1.8, 2.5, 2.25, 2.2, 2.55, 2.3, 2.1];
% y-coordinates in cm
y = [0, 0.04, 0.1, 0.42, 0.8, 0.68, 0.63, 0.88, 0.72, 0.53];
% y-coordinates in scaling factor
bode_data.centrABP_to_brachABP.amp.v = 1+(y*3/1.6);  % 1.6 cm between 1 and 4
clear y

%% "Brachial / Ascending Aortic" Phase
% x-coordinates in Hz
bode_data.centrABP_to_brachABP.phase.f = 0:9;
% y-coordinates (rough guess)
y = [0, -0.1, -0.25, -0.5, -1.4, -1.9, -2.4, -2.5, -2.9,-3.5];
% y-coordinates in cm
y = [0, -0.05, -0.15, -0.32, -0.75, -0.95, -1.23, -1.39, -1.63, -1.93];
% y-coordinates in rad
y = y*6/3.25;  % 3.25 cm between 0 and -6.
% y-coordinates in deg
bode_data.centrABP_to_brachABP.phase.v = y*360/(2*pi);
clear y

%% Brachial ABP to Central ABP
bode_data.brachABP_to_centrABP = bode_data.centrABP_to_brachABP;
bode_data.brachABP_to_centrABP.amp.v = 1./bode_data.centrABP_to_brachABP.amp.v;
bode_data.brachABP_to_centrABP.phase.v = -1*bode_data.centrABP_to_brachABP.phase.v;

%% "Radial / Ascending Aortic" Amplitude
% x-coordinates in Hz
bode_data.centrABP_to_radiaABP.amp.f = 0:9;
% y-coordinates in cm
y = [0, 0.05, 0.22, 0.48, 0.98, 0.64, 0.38, 0.45, 0.32, 0.22];
% y-coordinates in scaling factor
bode_data.centrABP_to_radiaABP.amp.v = 1+(y*3/1.4);  % 1.4 cm between 1 and 4
clear y

%% "Radial / Ascending Aortic" Phase
% x-coordinates in Hz
bode_data.centrABP_to_radiaABP.phase.f = 0:9;
% y-coordinates in cm
y = [0, -0.16, -0.3, -0.55, -1.1, -1.64, -1.98, -2.28, -2.73, -2.98];
% y-coordinates in rad
y = y*8/3.62;  % 3.62 cm between 0 and -8.
% y-coordinates in deg
bode_data.centrABP_to_radiaABP.phase.v = y*360/(2*pi);
clear y

%% Radial ABP to Central ABP
bode_data.radiaABP_to_centrABP = bode_data.centrABP_to_radiaABP;
bode_data.radiaABP_to_centrABP.amp.v = 1./bode_data.centrABP_to_radiaABP.amp.v;
bode_data.radiaABP_to_centrABP.phase.v = -1*bode_data.centrABP_to_radiaABP.phase.v;


% insert data from Fig. 4 of:
% ?1. Karamanoglu, M.; Feneley, M. P. Derivation of the ascending aortic-carotid pressure transfer function with an arterial model. Am. J. Physiol. Circ. Physiol. 1996, 271, H2399?H2404, doi:10.1152/ajpheart.1996.271.6.H2399.
% (used averaged bold black lines)

%% "Carotid / Ascending Aortic" Amplitude
% x-coordinates in Hz
bode_data.centrABP_to_carotABP.amp.f = 0:12;
% y-coordinates in cm
y = [0, 0.1, 0.2, 0.5, 0.8, 1.02, 1.0, 0.75, 0.45, 0.15, -0.05, -0.22, -0.32];
% y-coordinates in scaling factor
bode_data.centrABP_to_carotABP.amp.v = 1+(y*1/1.9);  % 1.9 cm between 1 and 2
clear y

%% "Carotid / Ascending Aortic" Phase
% x-coordinates in Hz
bode_data.centrABP_to_carotABP.phase.f = 0:12;
% y-coordinates in cm
y = [0, -0.5, -1.2, -2.5, -5.8, -7.5, -11.3, -14.4, -17, -19, -20.5, -21.5, -23]./10;
% y-coordinates in rad
y = y*3/3.27;  % 3.27 cm between 0 and -3.
% y-coordinates in deg
bode_data.centrABP_to_carotABP.phase.v = y*360/(2*pi);
clear y

%% Carotid ABP to Central PPG
bode_data.carotABP_to_centrABP = bode_data.centrABP_to_carotABP;
bode_data.carotABP_to_centrABP.amp.v = 1./bode_data.centrABP_to_carotABP.amp.v;
bode_data.carotABP_to_centrABP.phase.v = -1*bode_data.centrABP_to_carotABP.phase.v;

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
    eval(['rel_sigs.input = ' up.transformations{transformation_no}(1:8) ';']);
    
    % Recover output signal from FFT of input signal
    N = up.n_fft;
    x = nan(N,1);
    x_bar = calc_x_bar(S_proc_fft);
    rel_sigs.output.v = nan(N-1,1);
    rel_sigs.output.fs = rel_sigs.input.fs;
    temp.real_x_bar = real(x_bar(1:(N/2)+1));
    temp.imag_x_bar = imag(x_bar(1:(N/2)+1));
    temp.bode_phase = (rel_bode_data.phase.v*2*pi/360);
        
    % ~~~~ new version 1 ~~~~
    const = (2*pi/N)*([0 : N/2]');
    const2 = temp.real_x_bar .* rel_bode_data.amp.v;
    const3 = temp.imag_x_bar .* rel_bode_data.amp.v;
    for i = 1 : N-1
        
        sum_cosines = sum( const2  .* cos((const*i) + temp.bode_phase) );
        
        sum_sines = sum( const3 .* sin((const*i) + temp.bode_phase) );
        
        % add sines and cosines
        rel_sigs.output.v(i) = sum_cosines + sum_sines;
        
    end
    % ~~~~ end new version 1 ~~~~
    
%     % ~~~~ old version ~~~~
%     tic
%     for i = 1 : N-1
%         
%         sum_cosines = sum( temp.real_x_bar .* rel_bode_data.amp.v .* cos((2*pi*([0 : N/2]')*i/N) + temp.bode_phase) );
%         
%         sum_sines = sum( temp.imag_x_bar .* rel_bode_data.amp.v .* sin((2*pi*([0 : N/2]')*i/N) + temp.bode_phase) );
%         
%         % add sines and cosines
%         rel_sigs.output.v(i) = sum_cosines + sum_sines;
%         
%     end
%     toc
%     clear i
%     clear sum_sines sum_cosines
%     % ~~~~ end old version ~~~~
    
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
%     curr_rel_bode_data.amp.v = zeros(size(req_freqs));
    curr_rel_bode_data.amp.v = ones(size(req_freqs));
    curr_rel_bode_data.amp.v(1:length(req_freqs_within_rel_range)) = interp1(curr_bode_data.amp.f, curr_bode_data.amp.v, req_freqs_within_rel_range);
    
    % phase
    curr_rel_bode_data.phase.f = req_freqs;
    req_freqs_within_rel_range = req_freqs(req_freqs<=curr_bode_data.phase.f(end));
%     curr_rel_bode_data.phase.v = zeros(size(req_freqs));
    curr_rel_bode_data.phase.v = curr_bode_data.phase.v(end)*ones(size(req_freqs));  % added 27th Sept 2019
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
        eval(['req_len = length(transformed_sig.' up.input_sig '.v);']);
        actual_len = length(transformed_sig.centrABP.v);
        if actual_len > req_len
            transformed_sig.centrABP.v = transformed_sig.centrABP.v(1:req_len);
        end
        
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
        
%         %% Scale to occupy range of 0 to 1
%         for sig_no = 1 : length(sig_names)
%             eval(['curr_sig = transformed_sig.' sig_names{sig_no} ';'])            
%             curr_sig.v = curr_sig.v-min(curr_sig.v);
%             curr_sig.v = curr_sig.v / max(curr_sig.v);
%             eval(['transformed_sig.' sig_names{sig_no} ' = curr_sig;'])
%         end
        
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
    %ylim([0 1.1])
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