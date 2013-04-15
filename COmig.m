 %{

Comig.m - Constant offset Kirchhoff migration in time and depth.
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

dcmp = 20;       % CMP-Distance [m]
dt = 2e-3;       % Samplinginterval [s]
nt = 1001;       % Number of samples
ns = 101;        % Number of traces
nh = 5;          % Number of offsets
Fs = 1/dt;       % Frequency sampling [Hz]
hmax = 1000;     % Maximum half-offset [m]
dh = 250;        % Offset increment [m]
vmin = 1850;     % Minimum test velocity [m/s]
vmax = 1850;     % Maximum test velocity [m/s]
vfinal = 1850;   % Final migration velocity [m/s]
dv = 100;         % Velocity increment [m/s]
aper = 1500;      % Aperturewidth [m]

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

%% Frequency analysis
NFFT = 2^nextpow2(nt);                    % calculate next 2^n to prepare adta for FFT
fdata = fft(mean(data(:,:,1),2),NFFT)/nt; % FFT of extended dataset
faxis = Fs/2*linspace(0,1,NFFT/2+1);      % Skaling of the x-axis

% Plot of the normalized frequency spektrum
fx = figure(1);
plot(faxis,abs(fdata(1:length(faxis)))/max(abs(fdata(1:length(faxis)))));
xlabel('Frequency [Hz]','Fontsize',24);
ylabel('Normalized Amplitude','Fontsize',24);
%title('Frequency analysis','Fontsize',24)
set(gca,'Fontsize',24)
set(fx, 'Position', [0 0 1280 1024] );    % Size of the new frame
axis ([0 75 0 1])
print('-dpng','freq.png');                % Outputfile of figure


%% Kirchhoff migration

% Half aperture
half_aper = round(.5*aper/dcmp);

% Initialising
h = 0:dh:hmax;
Kirchhoffdepth(1:nt,1:ns,1:nh)=0;
i_v = 0;
Kirchhoff(1:nt,1:ns,1:nh)=0;              %(time, CMP, offset)
Skala(1:nt,1:nh) = 0;
t=(0:nt-1)'*dt;

%% Loop over velocities
for v = vmin:dv:vmax;
    i_v = i_v+1;
    Kirchhoff(1:nt,1:ns,1:nh)=0;          %(time, CMP, offset)
    Skala(1:nt,1:nh) = 0;
    
    %% loop over half offsets
    for i_h = 1:nh
        [Kirchhoff(:,:,i_h), Skala(:,i_h)] = CO_kirch(filtdata(:,:,i_h), v, h(i_h), dt, dcmp, half_aper);
        Kirchhoff(1,:,i_h) = 0;
        Kirchhoffdepth(:,:,i_h) = interp1(Skala(:,1),Kirchhoff(:,:,i_h),Skala(:,i_h),'spline');
    end
    
    % CO-Gather for each velocity
    fx=figure(v);
    set(fx, 'Position', [0 0 1280 1024] );
    imagesc(((1:ns*nh)-1)*dcmp,Skala(:,1),Kirchhoffdepth(:,:),[-max(max(max(abs(Kirchhoffdepth(1:nt-5,:,:))))) max(max(max(abs(Kirchhoffdepth(1:nt-5,:,:)))))])
    %title('Tiefenmigration','Fontsize',24)
    xlabel('CMP','Fontsize',24)
    ylabel('Depth [m]','Fontsize',24)
    set(gca,'Fontsize',24)
    colormap([ones(101,1),(0:.01:1)',(0:.01:1)';(1:-.01:0)',(1:-.01:0)',ones(101,1)])    % polarized plot
    colorbar
    set(gca,'Fontsize',24)
    set(gca,'XTickLabel',{' 0 ','2 / 0','2 / 0','2 / 0','2 / 0',' 2 '})                  % reskaling x-axis
    print('-dpng',sprintf('v%g.png',v));
    
    if v == vfinal % If loop reaches the correct velocity (estimated with constant velocity scan)
        mig(1:nt,1:ns) = sum(Kirchhoffdepth,3); % summing CO-Gather
        
        % Plot of the migration result
        fx=figure(v+1);
        set(fx, 'Position', [0 0 1280 1024] );
        imagesc(((1:ns)-1)*dcmp,Skala(:,1),mig(:,:),[-max(max(abs(mig))) max(max(abs(mig)))])
        %title('Tiefenmigration','Fontsize',24)
        xlabel('CMP [m]','Fontsize',24)
        ylabel('Depth [m]','Fontsize',24)
        colormap([ones(101,1),(0:.01:1)',(0:.01:1)';(1:-.01:0)',(1:-.01:0)',ones(101,1)])
        colorbar
        set(gca,'Fontsize',24)
        print('-dpng',sprintf('v%g.png',v));
        
        % Plot trace 51 normalized
        figure
        plot(((1:nt)-1)*dt,filtdata(:,51,1)/max(filtdata(:,51,1)),'r')
        hold on
        plot(((1:nt)-1)*dt,mig(:,51)/max(mig(:,51)),'k')
        ylabel('Normalisierte Amplitude','Fontsize',24)
        xlabel('Zeit [s]','Fontsize',24)
        legend('SNR Input','SNR Migriert','Location','best')
        set(gca,'Fontsize',24)
        print('-dpng','SNRnorm.png');
        
        % Plot trace 51 not normalized
        figure
        plot(((1:nt)-1)*dt,filtdata(:,51,1),'r')
        hold on
        plot(((1:nt)-1)*dt,mig(:,51),'k')
        ylabel('Amplitude','Fontsize',24)
        xlabel('Zeit [s]','Fontsize',24)
        legend('SNR Input','SNR Migriert','Location','best')
        set(gca,'Fontsize',24)
        print('-dpng','SNRreal.png');
        
        % Calculation of Signal-to-Noise-Ratio
        SNRin = log(max(max(filtdata(:,:,1)))/mean(mean(abs(filtdata(100:200,:)))));
        SNRout = log(max(max(mig(:,:)))/mean(mean(abs(mig(100:200,:)))));
        
        fprintf('Verbesserung der Signal-to-Noise ratio von %f2 auf %f2\n',SNRin,SNRout)
        
        % Fileoutput of datamatrices
        dlmwrite('mig.dat',mig)
        dlmwrite('COGatherh0.dat',Kirchhoffdepth(:,:,1));
        dlmwrite('COGatherh250.dat',Kirchhoffdepth(:,:,2));
        dlmwrite('COGatherh500.dat',Kirchhoffdepth(:,:,3));
        dlmwrite('COGatherh750.dat',Kirchhoffdepth(:,:,4));
        dlmwrite('COGatherh1000.dat',Kirchhoffdepth(:,:,5));
        
    end
    
end 

%% Frequency analysis of migrated and summed data

%Use FFT and faxis vectors from above
migdata = fft(mean(mig,2),NFFT)/nt; % FFT of migrated dataset

% Plot of the normalized frequency spektrum
fy = figure;
plot(faxis,abs(migdata(1:length(faxis)))/max(abs(migdata(1:length(faxis)))));
xlabel('Frequency [Hz]','Fontsize',24);
ylabel('Normalized Amplitude','Fontsize',24);
set(gca,'Fontsize',24)
set(fy, 'Position', [0 0 1280 1024] );    % Size of the new frame
axis ([0 75 0 1])
print('-dpng','freq_mig.png');                % Outputfile of figure


%Input signal normalized
figure
plot(((1:nt)-1)*dt,filtdata(:,51,1)/max(filtdata(:,51,1)),'r')
hold on
plot(((1:nt)-1)*dt,data(:,51,1)/max(data(:,51,1)),'k')
ylabel('Normalisierte Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('Filtered data','Original data','Location','best')
set(gca,'Fontsize',24)
print('-dpng','wavelet.png');

%Input signal not normalized
figure
plot(((1:nt)-1)*dt,filtdata(:,51,1),'r')
hold on
plot(((1:nt)-1)*dt,data(:,51,1),'k')
ylabel('Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('Filtered data','Original data','Location','best')
set(gca,'Fontsize',24)
print('-dpng','wavelet_unnorm.png');
