%Spatial Audio Rendering of Underwater Audio Signals
% uses function and code from the following libraries:
% Acoustical Spherical Array Processing by Archontis Politis
% Array-Response-Simulator by Archontis Politis
% COMPASS by Archontis Politis
% Stefan Schucker, 2023
clear all; close all;

%Simulation settings for the tetrahedron hydrophone array

%Define the array properties of the hydrophone array
hydrophone_dir_deg = [120 -19.46;240 -19.46;360 -19.46; 0 90]; %placement of hydrophones on the sphere
hydrophone_dir_rad = hydrophone_dir_deg*pi/180; %convert to radians

%Define the numbers of hydrophone sensors
nHydrophone = 4;

%Define the radius of tetrahedron
edge_length = 0.5547; %distance between two hydrophones
R=sqrt(6)/4*edge_length; %Using the formula for the radius of a sphere that encloses a tetrahedron

%get a visual interpretation of the hydrophone array
plotHydrophoneArray(hydrophone_dir_deg, R); view(25,10);
h = gcf; h.Position(3) = 1.5*h.Position(3); h.Position(4) = 1.5*h.Position(4);

c = 1500; %speed of sound in salt water
f_max = 20000; %maximum frequency that can be heard by human auditory system
kR_max = 2*pi*f_max*R/c; %wavevector

% frequency vector
fs = 48000; %sampling frequency of the audio file
Lfilt = 2048; %length of STFT 
f = (0:Lfilt/2)'*fs/Lfilt; %discrete sampling points in frequency space
kR = 2*pi*f*R/c; %array of wave vector
nBins = Lfilt/2+1; %numbers of frequency bins
sht_order = floor(sqrt(nHydrophone-1));
arrayType = 'open'; %Open sphere Array

%%
%Read in the HDF5 input files that can be downloaded from zenodo
data = h5read('Source_Files/Raw_HDF5_Zenodo/6th January/2023-01-06--09-16-18--00-27-41.hdf5','/Signals')';
data = double(data);
% Resample to 
A_Format = resample(data,48000,192000);
%check for recording length in minutes
t_length_min = length(A_Format)/fs/60;

%%
% If you wish to use the raw input from zenodo as wave files in A-Format
% for example to process them in reaper or audacity you can save the
% downsampled version of it here

out_file_name = 'A-Format/Skjervoy_06_01_23_R3.wav';  
audiowrite(out_file_name, A_Format, fs, 'BitsPerSample', 32,'Title', 'Recording Skjervoy','Artist','Stefan Schucker','Comment', 'Mai 2023');


%%
% Before you process the whole recording you should be aware that this
% implementation is not optimized for run time. So it is highly recomended
% to limit the A-Format to about 5min length otherwise the parametric
% decoding will take very long time since it's run speed is growing not linear
% with the input size but exponentially

%set up a time window you want to render:

t_start_min = fs*60*5
t_end_min = fs*60*6
A_Format = A_Format(t_start_min:t_end_min,:);
length(A_Format)/fs/60


%%
% Create Encoding Filter Matrix
f_alias= 703;
%create  % real SH matrix for hydrophones
aziElev2aziPolar = @(dirs) [dirs(:,1) pi/2-dirs(:,2)]; % function to convert from azimuth-inclination to azimuth-elevation
Y_hydro = sqrt(4*pi) * getSH(sht_order, aziElev2aziPolar(hydrophone_dir_rad), 'real');
M_hydro2sh_sht = (1/nHydrophone)*Y_hydro';
[H_radInverse,h_radInverse]=arraySHTfiltersTheory_radInverse(R,4,1,2048,48000,0);
for kk=1:nBins
    M_hydro2sh_radinv(:,:,kk) = diag(replicatePerOrder(H_radInverse(kk,:),2))*M_hydro2sh_sht;
end

% Create Encoding Filter Matrix with Diffuse SoundField Equalization
M_diffcoh = getDiffCohMtxTheory(hydrophone_dir_rad, arrayType, R, 1, f, []);
% Apply diffuse-field equalization to the encoding filter matrix
M_hydro2sh_diffeq = arraySHTfilters_diffEQ(M_hydro2sh_radinv, M_diffcoh, f_alias, fs);

%%
addpath('COMPASS-ref-main/compass-lib')
addpath('COMPASS-ref-main/ext-lib/afSTFT')
pars.fs = fs;
% hop size, determines resolution - the number of bands will be hopsize+1
% (for uniform mode)
pars.hopSize = 128;
% length of buffer for real-time processing, or how many windows to 
% process at once and make available at each iteration for offline processing
pars.frameSize = 16*128;
% options {'uniform','hybrid','low_delay'} for 
% a) normal STFT bin resolution, 
% b) for hybrid with LF additional bands, 
% c) for hybrid with lower processing delay for real-time apps (but with 
% reduced alias suppression than 'hybrid', only for real-time apps)
pars.STFTmode = 'uniform';
switch pars.STFTmode
    case 'uniform'
        pars.centerfreq = (0:pars.hopSize)'*pars.fs/(2*pars.hopSize);
    case {'hybrid','low_delay'}
        switch pars.fs
            case 44100
                load('./data/afSTFTCenterFreq133.mat', 'afCenterFreq44100');
                pars.centerfreq = afCenterFreq44100;
            case 48000
                load('./data/afSTFTCenterFreq133.mat', 'afCenterFreq48000');
                pars.centerfreq = afCenterFreq48000;
        end
end

%  linear processing doing nothing
nCH_in = size(A_Format,2);
nBands = length(pars.centerfreq);
M_dec = M_hydro2sh_diffeq;
pars.M_dec = M_dec; % filter matrix mapping inputs to outputs

% process input
B_Format = afSTFT_template(A_Format, pars);

%%

%Parameter Definition for Decoder
addpath('./COMPASS-ref-main')
% Decoding for Headphones or speaker arrangement
tic
SHorder=1;
   % Configure analysis
    input_struct.fs         = 48000;
    input_struct.SHorder    = SHorder;
    input_struct.AMBformat  = 0; % 0:N3d, 1:SN3D, 2: FuMa  
    input_struct.rotation_ypr = [0 0 0];    
    analysis_struct = compass_analysis_init(input_struct);
    
    % Load input recording
    nSH = (input_struct.SHorder+1)^2;
    %ambisig = audioread('COMPASS-ref-main/resources/03_choir_st350_foa_fuma.wav');
    %ambisig = ambisig(:, 1:nSH);
    ambisig = B_Format;
    
    % Apply analysis
    [compass_signals, compass_parameters] = compass_analysis(ambisig, analysis_struct);
    
    % Configure synthesis
    load('COMPASS-ref-main/resources/Lspkr_example.mat', 'ls_dirs_deg_aziElev')
    output_struct.ls_dirs = ls_dirs_deg_aziElev;
    output_struct.mode = 1; % 0: loudspeaker rendering, 1: headphone monitoring of loudspeaker rendering
    output_struct.eq = 1;
    output_struct.streamBalance = 1;
    output_struct.decodeBalance = 1;
    output_struct.diffusionLevel = 0.5;
    output_struct.vbapNorm = 0;
    output_struct.vbapSpread = 30;
    if output_struct.mode
        sofa_struct = loadSofaFile('COMPASS-ref-main/resources/ownsurround2016_short_48k.sofa');
        output_struct.hrirs = sofa_struct.IR;
        output_struct.hrtf_fs = sofa_struct.IR_fs;
        output_struct.hrtf_dirs = sofa_struct.SourcePosition(1:2,:).';
    end
    synthesis_struct = compass_synthesis_init(analysis_struct, output_struct);
    
    % Apply synthesis
    [output_signals, synthesis_struct] = compass_synthesis(compass_signals, synthesis_struct, compass_parameters);
    if output_struct.mode==0
        output_signals = 0.5*output_signals/max(max(abs(output_signals)));
        audiowrite(['choir_o' num2str(SHorder) '_ls.wav'], output_signals, input_struct.fs);
    elseif output_struct.mode==1
        output_signals = 0.5*output_signals/max(max(abs(output_signals)));
        audiowrite(['Binaural_Rendering/Skjervoy_1' num2str(SHorder) '_hp.wav'], output_signals, input_struct.fs);
    toc
    end
   







