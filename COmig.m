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
vmin = 2625;     % Minimum test velocity [m/s]
vmax = 2625;     % Maximum test velocity [m/s]
vfinal = 2625;   % Final migration velocity [m/s]
dv = 100;        % Velocity increment [m/s]
aper = 200;      % Aperturewidth [m]
dz = 4;          % Depthsampling increment [m]
interpolate=0;   % 1 = use interpolation, 0 = use rounding
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

%% Frequency analysis
NFFT = 2^nextpow2(nt);                    % calculate next 2^n to prepare adta for FFT
fdata = fft(mean(data(:,:,1),2),NFFT)/nt; % FFT of extended dataset
faxis = Fs/2*linspace(0,1,NFFT/2+1);      % Skaling of the x-axis

%% Frequency Plot


%% Kirchhoff migration

%% discretizing and initializing
h=0:dh:hmax;                       % half offset
t_orig=0:dt:((nt-1)*dt);
i_v=0;
aper_half = round(.5*aper/dcmp);      % half aperture
z_len = length(z);
% x-sampling = cmp-sampling

%% Loop over velocities
for v = vmin:dv:vmax;
    i_v = i_v+1;
    
    t_depth=t_orig*v*0.5;                   % TWT-time to depth conversion
    zmax = max(t_depth);                    % Max depth [m]
    z=0:dz:zmax;                            % Depthsampling
    Kirchhofftime(1:nt,1:ns,1:nh)=0;        % (time, CMP, halfoffset)
    Kirchhoffdepth(1:z_len,1:ns,1:nh)=0;% (depth, CMP, halfoffset)
    mig(1:z_len,1:ns)=0;                % stacking result
    filt_interp(1:z_len,1:ns,1:nh)=0;
    
    %% loop over half offsets
    for i_h = 1:nh
        disp('half offset'); disp((i_h-1)*dh);
        % to get the feeling it still runs when using
        % interpolation at zdiff
        if kirch_time == 1
        Kirchhofftime(:,:,i_h) = CO_kirch_time(filtdata(:,:,i_h), v, h(i_h), dt, dcmp, half_aper);
        z = sqrt((h/(v)).^2+((0:nt-1)'*dt).^2)*v*0.5; % Depth skaling
        Kirchhoffdepth(:,:,i_h) = interp1(z(:,1),Kirchhoff(:,:,i_h),z(:,i_h),'spline');
        end
        if kirch_depth == 1
        [Kirchhoffdepth(:,:,i_h)] = CO_kirch_depth(filtdata(:,:,i_h), v, h(i_h), dt, dcmp, half_aper);
        end
    end
    
    %% CO-Gather for each velocity
    
    if v == vfinal % If loop reaches the correct velocity (estimated with constant velocity scan)
        % (estimated with constant velocity scan)
        mig(1:z_len,1:ns) = sum(COG,3); % summing CO-Gather
        
        %% Plot of the migration result
        
        %% SNR plot
        
        % Calculation of Signal-to-Noise-Ratio
        SNRin = log(max(max(filtdata(:,:,1)))/mean(mean(abs(filtdata(100:200,:)))));
        SNRorig = log(max(max(abs(data(:,:,1))))/mean(mean(abs(data(100:200,:)))));
        SNRout = log(max(max(mig(:,:)))/mean(mean(abs(mig(100:200,:)))));
        
        fprintf('Signal-to-Noise ratio von %f2 (filtered) bzw. %f2 (original) auf %f2\n',SNRin,SNRorig,SNRout)
        
        %% Ergebnis Plots
        
        
        %% Frequency analysis of migrated and summed data
        %Use FFT and faxis vectors from above
        migdata = fft(mean(mig,2),NFFT)/nt; % FFT of migrated dataset
        
        %% Frequency Plot of migraTION
        
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