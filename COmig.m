%{

Comig.m - Common offset Kirchhoff migration in time and depth.
Copyright (C) 2013 Jesper S Dramsch, Matthias Schneider, Dela Spickermann, Jan Walda

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
%}

clear all        % Clear workspace
close all        % Close figures
clc              % Clear command line window

format long      % Double precision

tStart = tic;    % runtime measurement

dcmp = 20;       % CMP-Distance [m]
dt = 2e-3;       % Samplinginterval [s]
nt = 1001;       % Number of samples
ns = 101;        % Number of traces
nh = 5;          % Number of offsets
Fs = 1/dt;       % Frequency sampling [Hz]
hmax = 1000;     % Maximum half-offset [m]
dh = 250;        % half-offset increment [m]
vmin = 2754;     % Minimum test velocity [m/s]
vmax = 2754;     % Maximum test velocity [m/s]
vfinal = 2754;   % Final migration velocity [m/s]
dv = 100;        % Velocity increment [m/s]
aper = 400;      % Aperturewidth [m]
dz = vfinal/1e3; % Depthsampling increment [m]
flag_interp = 0; % 1 = use interpolation, 0 = use rounding
kirch_time=0;    % Time Migration
kirch_depth=1;   % Depth Migration

%% Open file
% Original data
filename = 'data-6/SEIS-orig';

fid = fopen(filename,'r');
data = reshape(fread(fid, [nt nh*ns],'single'),nt,ns,nh);
fclose(fid);

% With sqrt{-i omega} filtered data
filenamefilt = 'data-6/SEIS-filt';

fidfilt = fopen(filenamefilt,'r');
filtdata = reshape(fread(fidfilt, [nt nh*ns],'single'),nt,ns,nh);
fclose(fidfilt);


%% Input Plots
% Compare max amplitudes of each CMP in original data
mig_graphs('OffsetLine',data,((1:ns)-1)*dcmp,'ampsorig')
% Compare max amplitudes of each CMP in original data
mig_graphs('OffsetLine',filtdata,((1:ns)-1)*dcmp,'ampsfilt')
% Input signal normalized
mig_graphs('CompLine','Filtered data',filtdata(:,51,1)/max(filtdata(:,51,1)),'Original data',data(:,51,1)/max(data(:,51,1)),((1:nt)-1)*dcmp,'Zeit [s]','Normalisierte Amplitude','waveletnorm')
% Input signal not normalized
mig_graphs('CompLine','Filtered data',filtdata(:,51,1),'Original data',data(:,51,1),((1:nt)-1)*dcmp,'Zeit [s]','Amplitude','waveletorig')
% Input signal original
mig_graphs('COG',data(:,:),((1:ns*nh)-1)*dcmp,(0:nt-1)','Time [s]','InputCOG')
% Input signal filtered
mig_graphs('COG',filtdata(:,:),((1:ns*nh)-1)*dcmp,(0:nt-1)','Time [s]','InputCOGfilt')

%% Frequency analysis
NFFT = 2^nextpow2(nt);                    % calculate next 2^n to prepare adta for FFT
fdata = fft(mean(data(:,:,1),2),NFFT)/nt; % FFT of extended dataset
faxis = Fs/2*linspace(0,1,NFFT/2+1);      % Skaling of the x-axis

%% Frequency Plot
mig_graphs('SingleLine',abs(fdata(1:length(faxis)))/max(abs(fdata(1:length(faxis)))),faxis,'Frequencz [Hz]','Normalized Amplitude','freq')

%% Kirchhoff migration

%% discretizing and initializing
h=0:dh:hmax;                       % half offset
t_orig=0:dt:((nt-1)*dt);
i_v=0;
aper_half = round(.5*aper/dcmp);      % half aperture
% x-sampling = cmp-sampling

%% Loop over velocities
for v = vmin:dv:vmax;
    i_v = i_v+1;
    
    t_depth=t_orig*v*0.5;                   % TWT-time to depth conversion
    zmax = max(t_depth);                    % Max depth [m]
    z=0:dz:zmax;                            % Depthsampling
    z_len = length(z);
    Kirchhofftime(1:nt,1:ns,1:nh)=0;        % (time, CMP, halfoffset)
    %Kirchhoffdepth(1:z_len,1:ns,1:nh)=0;% (depth, CMP, halfoffset)
    mig(1:z_len,1:ns)=0;                % stacking result
    filt_interp(1:z_len,1:ns,1:nh)=0;
    
    %% loop over half offsets
    for i_h = 1:nh
        disp('half offset'); disp((i_h-1)*dh);
        % to get the feeling it still runs when using
        % interpolation at zdiff
        if kirch_time == 1
            [Kirchhofftime(:,:,i_h), Skala(:,i_h)] = CO_kirch_time(filtdata(:,:,i_h), v, h(i_h), dt, dcmp, aper_half, flag_interp);
            Kirchhoffdepth(:,:,i_h) = interp1(Skala(:,1),Kirchhofftime(:,:,i_h),Skala(:,i_h),'spline');
            z = Skala(:,1);
            z_len = length(z);
            z_max=max(z);
            COG = Kirchhoffdepth;
        end
        if kirch_depth == 1
            [Kirchhoffdepth(:,:,i_h)] = CO_kirch_depth(filtdata(:,:,i_h), v, h(i_h), dt, dz, dcmp, aper_half, flag_interp);
            COG = Kirchhoffdepth/(2*pi); %Was schlaueres ueberlegen
        end
    end
    
    %% CO-Gather for each velocity
    mig_graphs('COG',COG(:,:),((1:ns*nh)-1)*dcmp,z,'Depth [km]',sprintf('COGv%g',v))
    % Compare max amplitude in each CMP in COG
    mig_graphs('OffsetLine',COG,((1:ns)-1)*dcmp,sprintf('ampsv%g',v))
    
    if v == vfinal % If loop reaches the correct velocity (estimated with constant velocity scan)
        % (estimated with constant velocity scan)
        mig(1:z_len,1:ns) = mean(COG,3); % summing CO-Gather
        
        %% Plot of the migration result
        mig_graphs('PolarPlot',mig(:,:),((1:ns)-1)*dcmp/1000,z,'CMP [km]','Depth [km]',sprintf('sum_vm%g',v))
        % max amplitude of each CMP in migrated section
        mig_graphs('SingleLine',max(abs(mig(:,:))),((1:ns)-1)*dcmp,'CMP','Maximum Amplitude','migamp')
        
        %% Calculation of Signal-to-Noise-Ratio
        SNRin = log(max(max(filtdata(:,:,1)))/mean(mean(abs(filtdata(100:200,:)))));
        SNRorig = log(max(max(abs(data(:,:,1))))/mean(mean(abs(data(100:200,:)))));
        SNRout = log(max(max(mig(:,:)))/mean(mean(abs(mig(100:200,:)))));
        if kirch_time == 1
            % Plot trace 51 normalized
            mig_graphs('CompLine','SNR Input',filtdata(:,51,1)/max(filtdata(:,51,1)),'SNR Migriert',mig(:,51)/max(mig(:,51)),((1:nt)-1)*dt,'Zeit [s]','Normalized Amplitude','SNRnorm')
            % Plot trace 51 not normalized
            mig_graphs('CompLine','SNR Input',filtdata(:,51,1),'SNR Migriert',mig(:,51),((1:nt)-1)*dt,'Zeit [s]','Amplitude','SNRreal')
            
            % Comparison results normalized
            mig_graphs('CompLine','Original Data',data(:,51,1)/max(data(:,51,1)),'Migrated Data',mig(:,51,1)/max(mig(:,51,1)),((1:nt)-1)*dt,'Zeit [s]','Normalized Amplitude','waveletNorm')
            % Comparison results not normalized
            mig_graphs('CompLine','Original Data',data(:,51,1),'Migrated Data',mig(:,51,1),((1:nt)-1)*dt,'Zeit [s]','Amplitude','wavelet')
            
        end
        % Output to screen
        fprintf('Signal-to-Noise ratio von %f2 (filtered) bzw. %f2 (original) auf %f2\n',SNRin,SNRorig,SNRout)
        
        %% Ergebnis Plots
        % Input signal normalized
        mig_graphs('CompLine','Original Data',data(:,51,1)/max(data(:,51,1)),'Filtered Data',filtdata(:,51,1)/max(filtdata(:,51,1)),((1:nt)-1)*dt,'Zeit [s]','Normalized Amplitude','inoutNorm')
        % Input signal not normalized
        mig_graphs('CompLine','Original Data',data(:,51,1),'Filtered Data',filtdata(:,51,1),((1:nt)-1)*dt,'Zeit [s]','Amplitude','inout')
        
        
        %% Frequency analysis of migrated and summed data
        %Use FFT and faxis vectors from above
        migdata = fft(mean(mig,2),NFFT)/nt; % FFT of migrated dataset
        
        %% Frequency Plot of migration
        mig_graphs('SingleLine',abs(migdata(1:length(faxis)))/max(abs(migdata(1:length(faxis)))),faxis,'Frequencz [Hz]','Normalized Amplitude','freq_mig')
        % Frequency Comparison
        mig_graphs('CompLine','Original data',abs(fdata(1:length(faxis)))/max(abs(fdata(1:length(faxis)))),'Migrated data',abs(migdata(1:length(faxis)))/max(abs(migdata(1:length(faxis)))),faxis,'Frequency [Hz]','Normalized Amplitude','freq_comp')
        % Frequency Comparison
        mig_graphs('CompLine','Original data',abs(fdata(1:length(faxis))),'Migrated data',abs(migdata(1:length(faxis))),faxis,'Frequency [Hz]','Amplitude','freq_comp')
        % Fileoutput of datamatrices
        dlmwrite('output/mig.dat',mig)
        dlmwrite('output/COGatherh0.dat',COG(:,:,1));
        dlmwrite('output/COGatherh250.dat',COG(:,:,2));
        dlmwrite('output/COGatherh500.dat',COG(:,:,3));
        dlmwrite('output/COGatherh750.dat',COG(:,:,4));
        dlmwrite('output/COGatherh1000.dat',COG(:,:,5));
        
        %%% Für alte Amplitudenvergleiche müsste zurück interpoliert werden, da
        %%% z_len != length(t) und z_len hängt von v ab, siehe
        %%% t_depth=t_orig*v*0.5;
        %%% zmax = max(t_depth);     % zmax nimmt mit steigendem v zu, aber
        %%% z=0:dz:zmax;             % dz bleibt gleich !
        
    end
    
end
tElapsed = toc(tStart);               % calculation time