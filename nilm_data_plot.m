% NILM data analyze
% the main driver program for NILM analysis and plotting

% Bruce Madigan, Gray Beards Engineering, May, 2023
% This software was written for the specific purpose of 
% simple processing and plotting of NILM power data for reporting purposes.
% No warranty is provided whatsoever as to the correctness nor usefulness of this code
% Original intention was to have the code run under both Octave and Matlab without
% any modifications, but many changes that occurred likely broke the matlab compatibility
% therefore code should always run under octave

close all % close open plots
clear all % clear all variables
%see if we are running in Matlab or Octave 
try 
    OCTAVE_VERSION;
    inOctave = 1;
catch
    inOctave = 0;
end

if(inOctave==1)
  pkg load image
  pkg load signal
  pkg load io
  %graphics_toolkit ('gnuplot');
  %graphics_toolkit ('fltk');
  graphics_toolkit ('qt');
end

graphics_toolkit

b_save_plots=1; %flag to write plot images to disk
%user='bmadigan'; % until a better way is found, change this to the username being used
%[FNAME, FPATH, FLTIDX] = uigetfile (['/home/' user '/Documents/' '*.csv']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single file mode  (comment out for directory mode below)
[FNAME, FPATH, FLTIDX] = uigetfile ('*.csv');
% display the instantanious voltage and current data plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% my decision to put the column labels as the first row causes issues to automatically parse data
% so you have to know the structure of the file to make the code below work
% issues araise in that changes to the data collection must be known and changed here
% I have a version using cells and structures but the code is horrible so I punted and 
% hardcoded the labels - beware if the data collection file structure changes!


% for development it is easier just to hard code a data file to work with
   %FNAME='2023_01_08_17_41_17.csv';
   %raw_file = csv2cell('/home/bmadigan/Documents/2023_01_08_17_41_17.csv', ',');
   
% the data files are stored by default in the current users 'Documents' directory
%user='bmadigan'; % until a better way is found, change this to the username being used
%[FNAME, FPATH, FLTIDX] = uigetfile (['/home/' user '/Documents/' '*.csv']);

fid = fopen ([FPATH FNAME]);
txt = fgetl (fid) % display the labels in the output terminal
fclose (fid);

raw_data=dlmread([FPATH FNAME],','); % dlmread works nice but the first row must be removed
% get rid of the first row 
raw_data(1,:) = [];

% get the numeric data and split into vectors, here is where you need to know the file structure
t=raw_data(:,1);
PH1V = raw_data(:,2);
PH1I = raw_data(:,3);
PH2V = raw_data(:,4);
PH2I = raw_data(:,5);
PH3V = raw_data(:,6);
PH3I = raw_data(:,7);
ArcV = raw_data(:,8);
ArcI = raw_data(:,9);

% set power source name for title
idx = strfind (FPATH, '455');
if(~isempty(idx))
  PSNAME='Lincoln PowerWave 455 Inverter GMAW ';
end
idx = strfind (FPATH, '300');
if(~isempty(idx))
  PSNAME='Miller Dynasty 300 Inverter GTAW ';
end
idx = strfind (FPATH, '355');
if(~isempty(idx))
  PSNAME='Lincoln SquareWave 355 GTAW SCR ';
  % the lincoln is a single phase machine so weld only need the phase legs being used
  PH1V = PH2V-PH3V;
  PH1I = PH2I;
  PH2V = PH2V*0.0;
  PH3V = PH3V*0.0;
  PH3I = PH3I*0.0;
  ArcI = ArcI*-1.0; % had wire in sensor the other way
end

% during the arc start we may get negative ArcI/V values so clamp to 0
idx=find(ArcV<0);
ArcV(idx)=0;
%idx=find(ArcI<0);
%ArcI(idx)=0;

% limit the plot ranges so all plots have same range
PHV_min=-480;
PHV_max=480;
PHI_min=-35;
PHI_max=35;
ArcV_min=0;
ArcV_max=45;
ArcI_min=0;
ArcI_max=350;
% plot the input voltage and current data
% switch b_print_plot to 1 if you want to produce images for a report
b_print_plot=b_save_plots;
%FONTSIZE=12;

waveform_color=[0.85 0.85 0.85]; % color of instantaneous data (light gray)

% could use the data labels from the file, but since we know what they are, simpler to just use text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 1
subplot (2, 4, 1)
plot(t,PH1V,'color',waveform_color);
% calc moving RMS
Fs    = 1/diff(t(1:2)); % sampling frequency
% get window length in terms of sample frequency Fs
% and known data frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 
PH1V_rms = movfun (@rms, PH1V, num_pts_in_cycle);
hold on
plot(t,PH1V_rms,'-b','linewidth',4);
area (t,PH1V_rms,'facecolor',[0.85 0.85 1],'edgecolor','b','linewidth',2,'facealpha', 0.25);
xlabel('time (s)');
%ylabel('PH1V (V),PH1VRMS (V)');
ylabel('PH1V (V,VRMS)');
title_text=sprintf('%s','Phase 1 Voltage');
title(title_text,'interpreter', 'none');
grid on
ylim([PHV_min PHV_max]); % clamp the y range
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot (2, 4, 5)
plot(t,PH1I,'color',waveform_color);
% calc moving RMS
Fs    = 1/diff(t(1:2)); % sampling frequency
% get window length in terms of sample frequency Fs
% and known data frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 
PH1I_rms = movfun (@rms, PH1I, num_pts_in_cycle);
hold on
plot(t,PH1I_rms,'-b','linewidth',4);
%area (t,PH1I_rms,'facecolor',[0.85 0.85 1],'edgecolor','b','linewidth',2,'facealpha', 0.25);
xlabel('time (s)');
%ylabel('PH1I (A),PH1IRMS (A)');
ylabel('PH1I (A,ARMS)');
title_text=sprintf('%s','Phase 1 Current');
title(title_text,'interpreter', 'none');
grid on
ylim([PHI_min PHI_max]); % clamp the y range
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 2
subplot (2, 4, 2)
plot(t,PH2V,'color',waveform_color);
% calc moving RMS
Fs    = 1/diff(t(1:2)); % sampling frequency
% get window length in terms of sample frequency Fs
% and known data frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 
PH2V_rms = movfun (@rms, PH2V, num_pts_in_cycle);
hold on
plot(t,PH2V_rms,'-b','linewidth',4);
area (t,PH2V_rms,'facecolor',[0.85 0.85 1],'edgecolor','b','linewidth',2,'facealpha', 0.25);
xlabel('time (s)');
%ylabel('PH2V (V),PH2VRMS (V)');
ylabel('PH2V (V,VRMS)');
title_text=sprintf('%s','Phase 2 Voltage');
title(title_text,'interpreter', 'none');
grid on
ylim([PHV_min PHV_max]); % clamp the y range
text (0, 750, [PSNAME FNAME],'fontsize',40,'interpreter', 'none','fontweight','bold');

hold off


subplot (2, 4, 6)
plot(t,PH2I,'color',waveform_color);
% calc moving RMS
Fs    = 1/diff(t(1:2)); % sampling frequency
% get window length in terms of sample frequency Fs
% and known data frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 
PH2I_rms = movfun (@rms, PH2I, num_pts_in_cycle);
hold on
plot(t,PH2I_rms,'-b','linewidth',4);
%area (t,PH2I_rms,'facecolor',[0.85 0.85 1],'edgecolor','b','linewidth',2,'facealpha', 0.25);
xlabel('time (s)');
%ylabel('PH2I (A),PH2IRMS (A)');
ylabel('PH2I (A,ARMS)');
title_text=sprintf('%s','Phase 2 Current');
title(title_text,'interpreter', 'none');
grid on
ylim([PHI_min PHI_max]); % clamp the y range
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 3
subplot (2, 4, 3)
plot(t,PH3V,'color',waveform_color);
% calc moving RMS
Fs    = 1/diff(t(1:2)); % sampling frequency
% get window length in terms of sample frequency Fs
% and known data frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 
PH3V_rms = movfun (@rms, PH3V, num_pts_in_cycle);
hold on
plot(t,PH3V_rms,'-b','linewidth',4);
area (t,PH3V_rms,'facecolor',[0.85 0.85 1],'edgecolor','b','linewidth',2,'facealpha', 0.25);
xlabel('time (s)');
%ylabel('PH3V (V),PH3VRMS (V)');
ylabel('PH3V (V,VRMS)');
title_text=sprintf('%s','Phase 3 Voltage');
title(title_text,'interpreter', 'none');
grid on
ylim([PHV_min PHV_max]); % clamp the y range
hold off

subplot (2, 4, 7)
plot(t,PH3I,'color',waveform_color);
% calc moving RMS
Fs    = 1/diff(t(1:2)); % sampling frequency
% get window length in terms of sample frequency Fs
% and known data frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 
PH3I_rms = movfun (@rms, PH3I, num_pts_in_cycle);
hold on
plot(t,PH3I_rms,'-b','linewidth',4);
%area (t,PH3I_rms,'facecolor',[0.85 0.85 1],'edgecolor','b','linewidth',2,'facealpha', 0.25);
xlabel('time (s)');
%ylabel('PH3I (A),PH3IRMS (A)');
ylabel('PH3I (A,ARMS)');
title_text=sprintf('%s','Phase 3 Current');
title(title_text,'interpreter', 'none');
grid on
ylim([PHI_min PHI_max]); % clamp the y range
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arc
subplot (2, 4, 4)
plot(t,ArcV,'color',waveform_color);
% calc moving RMS
Fs    = 1/diff(t(1:2)); % sampling frequency
% get window length in terms of sample frequency Fs
% and known data frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 
ArcV_rms = movfun (@rms, ArcV, num_pts_in_cycle);
ArcV_mean = movfun (@rms, ArcV, num_pts_in_cycle);
hold on
plot(t,ArcV_rms,'-b','linewidth',4);

%plot(t,ArcV_mean,'-b','linewidth',4);
%area (t,ArcV_rms,'facecolor',[0.85 0.85 1],'edgecolor','b','linewidth',2);
xlabel('time (s)');
%ylabel('ArcV (V),ArcVRMS (V)');
ylabel('ArcV (V,VRMS)');
title_text=sprintf('%s','Arc Voltage');
title(title_text,'interpreter', 'none');
grid on
ylim([ArcV_min ArcV_max]); % clamp the y range
hold off


subplot (2, 4, 8)
plot(t,ArcI,'color',waveform_color);
% calc moving RMS
Fs    = 1/diff(t(1:2)); % sampling frequency
% get window length in terms of sample frequency Fs
% and known data frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 
ArcI_rms = movfun (@rms, ArcI, num_pts_in_cycle);
ArcI_mean = movfun (@mean, ArcI, num_pts_in_cycle);
hold on
plot(t,ArcI_rms,'-b','linewidth',4);
%plot(t,ArcI_mean,'-b','linewidth',4);
%area (t,ArcI_rms,'facecolor',[0.85 0.85 1],'edgecolor','b','linewidth',2);
xlabel('time (s)');
ylabel('ArcI (A,ARMS)');
%ylabel('ArcI (A),ArcIRMS (A)');
title_text=sprintf('%s','Arc Current');
title(title_text,'interpreter', 'none');
grid on
ylim([ArcI_min ArcI_max]); % clamp the y range
hold off


% save the plot to an image
if(b_print_plot==1)
plot_file_name=sprintf('%s.png',[FPATH FNAME '_all_data']);
print (plot_file_name, '-dpng','-F:8','-S1200,600');
  pause(0.1);
  fflush(stdout);

end


%now plot the input and output power
figure

% Instantaneous Input power for each phase
PH1Power=PH1V.*PH1I;
PH2Power=PH2V.*PH2I;
PH3Power=PH3V.*PH3I;

% get window length of a single 60 Hz cycle in terms of sample frequency Fs
Fs    = 1/diff(t(1:2)); % sampling frequency
% power frequency 60Hz
num_pts_in_cycle=round(Fs/60.0); 

% real power ie mean of Instantaneous
PH1Preal = movfun (@mean, PH1Power, num_pts_in_cycle); % Real power for each cycle
PH2Preal = movfun (@mean, PH2Power, num_pts_in_cycle); % Real power for each cycle
PH3Preal = movfun (@mean, PH3Power, num_pts_in_cycle); % Real power for each cycle

%total real power
if(inOctave==1)
    PHTPreal = PH1Preal .+ PH2Preal .+ PH3Preal;
else
    PHTPreal = PH1Preal + PH2Preal + PH3Preal;
end 
% call the input real power the input power
Input_Power=PHTPreal; %real
%Input_Power=PHTPavg; %average or instantenious


% for the arc, changes in I and V happen, so most info about welding is in output power
% but we are trying to use only input power

Output_Power=ArcV.*ArcI; %instantaneous
% get the output power with the same averaging as the input power
Output_Power=movfun (@mean, Output_Power, num_pts_in_cycle); % mean power for each cycle

% limit the plot ranges so all plots have same range
PHV_min=-480;
PHV_max=480;
PHI_min=-30;
PHI_max=30;
ArcV_min=0;
ArcV_max=40;
ArcI_min=0;
ArcI_max=300;
Input_Power_max=8000;
Output_Power_max=8000;
% plot the data
% there are many combinations to plot, just pick a few for a look, then develop more
% switch b_print_plot to 1 if you want to produce images for a report
b_print_plot=b_save_plots;
%FONTSIZE=12;

waveform_color=[0.0 0.0 0.85];

% set the horizontal axis to exclude data where arc is not on, like when we forgot to click stop
% find the indices where the Output power > 100 W
idx=find(Output_Power>100);
% get the time back a little from where the arc started
t_start_limit=t(idx(1))*0.9;
t_start_limit=floor(t(idx(1)));
% get the time aahead a little from where the arc stopped
t_stop_limit=ceil(t(idx(end)));

% get the peak output power from the center half portion
center=round((idx(end)-idx(1))/4);
Output_Power_rms = movfun (@rms, Output_Power, num_pts_in_cycle);
%Output_Power_max=max(Output_Power(idx(center:3*center)));
Output_Power_max=max(Output_Power_rms(idx(center:3*center)));

Input_Power_max=max(Input_Power(idx(center:3*center)));
%max_power=max(Output_Power_max,Input_Power_max);
Output_Power_max=Output_Power_max*1.1; % give a little head room on plot
Input_Power_max=Input_Power_max*1.1;



plot(t,Input_Power,t,Output_Power);
legend("Input Power","Output Power");
ylabel('Power(W)');
xlabel('time (s)','fontweight','bold');
grid on
%title_text=sprintf('%s','Input/Output Power (W)',[PSNAME FNAME]);
title_text=sprintf('%s',[PSNAME FNAME]);
title(title_text,'interpreter', 'none','fontsize',16);


% save the plot to an image
if(b_print_plot==1)
  plot_file_name=sprintf('%s.png',[FPATH FNAME '_weld_power_data']);
  print (plot_file_name, '-dpng');
  %print (plot_file_name, '-dpng','-F:8','-S1200,600');
  %print (plot_file_name, '-dpng','-F:8','-S850,1100');
  %print (plot_file_name, '-dpng','-F:8','-S850,1100');
  pause(0.1);
  fflush(stdout);

%  set(gcf, 'PaperPosition', [0.75 0.11 7 10]);
%  plot_file_name=sprintf('%s.pdf',[FPATH FNAME '_appendix_weld_data']);
%  print (plot_file_name, '-dpdf','-F:12','-append','-opengl');
end

% design a low pass filter
% the sample frequency
Fs    = 1/diff(t(1:2)); % sampling frequency from sampling period
% Nyquist frequency, in Hz.
% the Nyquist frequency is half the sampling frequency.
Fnyq = Fs/2

% the cut-off frequency of the Low pass filter in Hz.
% Fc frequency must be greater than 0 and less than Fnyq.
Fc=3

% create a first-order Butterworth low pass.
% The returned vectors are of legth n.
% Thus a first order filter is created with n = 2.
[b,a]=butter(2, Fc/Fnyq);


% Apply the filter to the input signal 
lp_Input_Power=filter(b,a,Input_Power);
lp_Output_Power=filter(b,a,Output_Power);

figure
plot(t,Input_Power,t,lp_Input_Power);
legend("Input Power","Filtered Input Power");
ylabel('Power(W)');
xlabel('time (s)','fontweight','bold');
grid on
%title_text=sprintf('%s','Input/Output Power (W)',[PSNAME FNAME]);
title_text=sprintf('%s',[PSNAME FNAME]);
title(title_text,'interpreter', 'none','fontsize',16);


% save the plot to an image
if(b_print_plot==1)
  plot_file_name=sprintf('%s.png',[FPATH FNAME '_input_power_data']);
  print (plot_file_name, '-dpng');
  %print (plot_file_name, '-dpng','-F:8','-S1200,600');
  %print (plot_file_name, '-dpng','-F:8','-S850,1100');
  %print (plot_file_name, '-dpng','-F:8','-S850,1100');
  pause(0.1);
  fflush(stdout);

%  set(gcf, 'PaperPosition', [0.75 0.11 7 10]);
%  plot_file_name=sprintf('%s.pdf',[FPATH FNAME '_appendix_weld_data']);
%  print (plot_file_name, '-dpdf','-F:12','-append','-opengl');
end

figure
plot(t,Input_Power,t,lp_Input_Power);
legend("Output Power","Filtered Output Power");
ylabel('Power(W)');
xlabel('time (s)','fontweight','bold');
grid on
%title_text=sprintf('%s','Input/Output Power (W)',[PSNAME FNAME]);
title_text=sprintf('%s',[PSNAME FNAME]);
title(title_text,'interpreter', 'none','fontsize',16);


% save the plot to an image
if(b_print_plot==1)
  plot_file_name=sprintf('%s.png',[FPATH FNAME '_output_power_data']);
  print (plot_file_name, '-dpng');
  %print (plot_file_name, '-dpng','-F:8','-S1200,600');
  %print (plot_file_name, '-dpng','-F:8','-S850,1100');
  %print (plot_file_name, '-dpng','-F:8','-S850,1100');
  pause(0.1);
  fflush(stdout);

%  set(gcf, 'PaperPosition', [0.75 0.11 7 10]);
%  plot_file_name=sprintf('%s.pdf',[FPATH FNAME '_appendix_weld_data']);
%  print (plot_file_name, '-dpdf','-F:12','-append','-opengl');
end
