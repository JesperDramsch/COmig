 %{

COGkirch.m - Constant offset Kirchhoff migration.
Copyright (C) 2013 Jesper S Dramsch, Matthias Schneider, 
                   Dela Spickermann, Jan Walda

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

%% Open file
% Original data
filename = 'data-7/SEIS-orig';

fid = fopen(filename,'r');
data = reshape(fread(fid, [nt nh*ns],'single'),nt,ns,nh);
fclose(fid);

% With sqrt{-i omega} filtered data
filenamefilt = 'data-7/SEIS-filt';

fidfilt = fopen(filenamefilt,'r');
filtdata = reshape(fread(fidfilt, [nt nh*ns],'single'),nt,ns,nh);
fclose(fidfilt);

%% Input data plots

% Compare max amplitudes of each CMP in original data
figure
plot(((1:ns)-1)*dcmp,max(abs(data(:,:,1)),[],1),'k')
hold on
plot(((1:ns)-1)*dcmp,max(abs(data(:,:,2)),[],1),'b')
hold on
plot(((1:ns)-1)*dcmp,max(abs(data(:,:,3)),[],1),'r')
hold on
plot(((1:ns)-1)*dcmp,max(abs(data(:,:,4)),[],1),'g')
hold on
plot(((1:ns)-1)*dcmp,max(abs(data(:,:,5)),[],1),'y')
xlabel('CMP','Fontsize',24)
ylabel('Maximum amplitude','Fontsize',24)  
legend('h = 0','h = 250','h = 500','h = 750','h = 1000','Location','best')
set(gca,'Fontsize',24)
print('-dpng',sprintf('output/ampsorig.png'));

% Compare max amplitudes of each CMP in original data
figure
plot(((1:ns)-1)*dcmp,max(abs(filtdata(:,:,1)),[],1),'k')
hold on
plot(((1:ns)-1)*dcmp,max(abs(filtdata(:,:,2)),[],1),'b')
hold on
plot(((1:ns)-1)*dcmp,max(abs(filtdata(:,:,3)),[],1),'r')
hold on
plot(((1:ns)-1)*dcmp,max(abs(filtdata(:,:,4)),[],1),'g')
hold on
plot(((1:ns)-1)*dcmp,max(abs(filtdata(:,:,5)),[],1),'y')
xlabel('CMP','Fontsize',24)
ylabel('Maximum amplitude','Fontsize',24)  
legend('h = 0','h = 250','h = 500','h = 750','h = 1000','Location','best')
set(gca,'Fontsize',24)
print('-dpng',sprintf('output/ampsfilt.png'));

% Input signal normalized
figure
plot(((1:nt)-1)*dt,filtdata(:,51,1)/max(filtdata(:,51,1)),'r')
hold on
plot(((1:nt)-1)*dt,data(:,51,1)/max(data(:,51,1)),'k')
ylabel('Normalisierte Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('Filtered data','Original data','Location','best')
set(gca,'Fontsize',24)
print('-dpng','output/waveletnorm.png');

% Input signal not normalized
figure
plot(((1:nt)-1)*dt,filtdata(:,51,1),'r')
hold on
plot(((1:nt)-1)*dt,data(:,51,1),'k')
ylabel('Normalisierte Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('Filtered data','Original data','Location','best')
set(gca,'Fontsize',24)
print('-dpng','output/waveletorig.png');
        
%% Frequency analysis
% NFFT = 2^nextpow2(nt);                    % calc next 2^n for FFT
% fdata = fft(mean(data(:,:,1),2),NFFT)/nt; % FFT of extended dataset
% faxis = Fs/2*linspace(0,1,NFFT/2+1);      % Skaling of the x-axis
% 
% % Plot of the normalized frequency spektrum
% fx = figure(1);
% plot(faxis,abs(fdata(1:length(faxis)))/max(abs(fdata(1:length(faxis)))));
% xlabel('Frequency','Fontsize',24);
% ylabel('Normalized Amplitude','Fontsize',24);
% %title('Frequency analysis','Fontsize',24)
% set(gca,'Fontsize',24)
% set(fx, 'Position', [0 0 1280 1024] );    % Size of the new frame
% axis ([0 75 0 1])
% print('-dpng','output/freq.png');         % Outputfile of figure

%% discretizing and initializing
h=0:dh:hmax;                       % half offset
t_orig=0:dt:((nt-1)*dt);
i_v=0;
h_aper = round(.5*aper/dcmp);      % half aperture
% x-sampling = cmp-sampling

%% Loop over velocities
for v = vmin:dv:vmax;
    i_v = i_v+1;

    t_depth=t_orig*v*0.5;              % TWT-time to depth conversion
    zmax = max(t_depth);               % Max depth [m]    
    z=0:dz:zmax;                       % Depthsampling
    COG(1:length(z),1:ns,1:nh)=0;      % (depth, CMP, halfoffset)
    mig(1:length(z),1:ns)=0;           % stacking result
    filt_interp(1:length(z),1:ns,1:nh)=0;

    %% Loop over offsets
    for i_h = 1:nh
        disp('half offset'); disp((i_h-1)*dh);
        % to get the feeling it still runs when using 
        % interpolation at zdiff
        
        %% Interpolation from t to z domain
        % must be calculated before because aperture may need data in front
        % of the current cmp
        for i_cmp = 1:ns            
            filt_interp(:,i_cmp,i_h) = interp1(t_depth,...
            filtdata(:,i_cmp,i_h),z,'spline');
        end
        
        %% Loop over CMPs
        for i_cmp = 1:ns
            % Aperture limits
            bound_l = max(floor(i_cmp-h_aper), 1);
            bound_r = min(floor(i_cmp+h_aper), ns);
            
            % Control if everything runs smoothly
            disp(['     CMP ||' ' left boundary ||'...
                ' right boundary ||' ' half aperture ||' ' velocity']);
            disp([(i_cmp-1)*dcmp (bound_l-1)*dcmp...
                (bound_r-1)*dcmp h_aper*dcmp v]);
            %% Loop over contributing samples (Aperture)
            for i_aper=bound_l:bound_r
                
                %% Loop over Depth
                for i_z=1:length(z)
                    
                    % Compute diffraction hyperbola, /2 because data is TWT
                    zdiff = 0.5*( sqrt(z(i_z)^2 ...
                        + ((i_cmp-i_aper)*dcmp-h(i_h)).^2) ...
                        + sqrt(z(i_z)^2 ...
                        + ((i_cmp-i_aper)*dcmp+h(i_h)).^2) );
                    
                    % Exit if diffraction ist out of data
                    if(zdiff > (max(z)))
                        break;
                    end
                    
                    %% Interpolate zdiff
                    % ! only if with interpolation at zdiff
                    if(interpolate==1)
                        res_interp = interp1(z,...
                            filt_interp(:,i_aper,i_h),zdiff,'spline');
                    end
                    %% Compute amplitude correction
                    cosphi = z(i_z)/sqrt(z(i_z)^2+h(i_h)^2);
                    weight = cosphi/sqrt(zdiff*v);
                    
                    %% Sum up along diffraction
                    % ! with interpolation at zdiff
                    if(interpolate==1)
                        COG(i_z,i_cmp,i_h) = COG(i_z,i_cmp,i_h) ...
                            + res_interp * weight;
                    elseif(interpolate==0)
                        i_zdiff = floor(1.5+zdiff/dz); 
                        % +0.5 so it get rounded correctly and + 1 so its 
                        % start with index 1, +1+0.5 = +1.5
                        
                        % ! without interpolation at zdiff
                        COG(i_z,i_cmp,i_h) = COG(i_z,i_cmp,i_h) ...
                            + filt_interp(i_zdiff,i_aper,i_h) * weight;
                    else
                        disp('Error, no valid method!');
                    end
                end
            end
        end
        COG(1,:,i_h) = 0;         % NaN avoiding
    end
    
    %% CO-Gather for each velocity
    figure
    imagesc(((1:ns*nh)-1)*dcmp,z,COG(:,:),...
        [-max(max(abs(COG(:,:,1)))) +max(max(abs(COG(:,:,1))))])
    %title('Tiefenmigration','Fontsize',24)
    xlabel('CMP','Fontsize',24)
    ylabel('Depth','Fontsize',24)
    set(gca,'Fontsize',24)
    colormap([ones(101,1),(0:.01:1)',(0:.01:1)';(1:-.01:0)',(1:-.01:0)',...
        ones(101,1)])    % polarized plot
    colorbar
    set(gca,'Fontsize',24)
    set(gca,'XTickLabel',{' 0 ','2 / 0','2 / 0','2 / 0','2 / 0',' 2 '})                  % reskaling x-axis
    print('-dpng',sprintf('output/COGv%g.png',v));    
    
    % Compare max amplitudes of each CMP in COG
    figure
    plot(((1:ns)-1)*dcmp,max(abs(COG(:,:,1)),[],1),'k')
    hold on
    plot(((1:ns)-1)*dcmp,max(abs(COG(:,:,2)),[],1),'b')
    hold on
    plot(((1:ns)-1)*dcmp,max(abs(COG(:,:,3)),[],1),'r')
    hold on
    plot(((1:ns)-1)*dcmp,max(abs(COG(:,:,4)),[],1),'g')
    hold on
    plot(((1:ns)-1)*dcmp,max(abs(COG(:,:,5)),[],1),'y')
    xlabel('CMP','Fontsize',24)
    ylabel('Maximum amplitude','Fontsize',24)  
    legend('h = 0','h = 250','h = 500','h = 750','h = 1000','Location','best')
    set(gca,'Fontsize',24)
    print('-dpng',sprintf('output/ampsv%g.png',v));
        
    if v == vfinal % If loop reaches the correct velocity 
                   % (estimated with constant velocity scan)
        mig(1:length(z),1:ns) = sum(COG,3); % summing CO-Gather
        
        %% Plot of the migration result
        figure
        imagesc(((1:ns)-1)*dcmp,z,mig(:,:),...
            [-max(max(abs(mig))) +max(max(abs(mig)))])
        %title('Tiefenmigration','Fontsize',24)
        xlabel('CMP','Fontsize',24)
        ylabel('Depth','Fontsize',24)
        colormap([ones(101,1),(0:.01:1)',(0:.01:1)';(1:-.01:0)',...
            (1:-.01:0)',ones(101,1)])
        colorbar
        set(gca,'Fontsize',24)
        print('-dpng',sprintf('output/vm%g.png',v));
        
        % max amplitude of each CMP in migrated section
        figure
        plot(((1:ns)-1)*dcmp,max(abs(mig(:,:)),[],1),'k')
        xlabel('CMP','Fontsize',24)
        ylabel('Maximum amplitude','Fontsize',24)    
        set(gca,'Fontsize',24)
        print('-dpng',sprintf('output/migamp.png'));
        
        %%% Für alte traceplots müsste zurück interpoliert werden, da
        %%% length(z) != length(t) und length(z) hängt von v ab, siehe
        %%% t_depth=t_orig*v*0.5;  
        %%% zmax = max(t_depth);     % zmax nimmt mit steigendem v zu, aber   
        %%% z=0:dz:zmax;             % dz bleibt gleich !
        
        % Calculation of Signal-to-Noise-Ratio
        SNRin = log(max(max(abs(filtdata(:,:,1))))/...
            mean(mean(abs(filtdata(100:200,:)))));
        SNRorig = log(max(max(abs(data(:,:,1))))/...
            mean(mean(abs(data(100:200,:)))));
        SNRout = log(max(max(abs(mig(:,:))))/...
            mean(mean(abs(mig(100:200,:)))));
        
        fprintf('Signal-to-Noise ratio von %f2 (filtered) bzw. %f2 (original) auf %f2\n'...
            ,SNRin,SNRorig,SNRout)
        
        % Fileoutput of datamatrices
        dlmwrite('output/mig.dat',mig)
        dlmwrite('output/COGatherh0.dat',COG(:,:,1));
        dlmwrite('output/COGatherh250.dat',COG(:,:,2));
        dlmwrite('output/COGatherh500.dat',COG(:,:,3));
        dlmwrite('output/COGatherh750.dat',COG(:,:,4));
        dlmwrite('output/COGatherh1000.dat',COG(:,:,5));

        %%% Für alte Amplitudenvergleiche müsste zurück interpoliert werden, da
        %%% length(z) != length(t) und length(z) hängt von v ab, siehe
        %%% t_depth=t_orig*v*0.5;  
        %%% zmax = max(t_depth);     % zmax nimmt mit steigendem v zu, aber   
        %%% z=0:dz:zmax;             % dz bleibt gleich ! 
        
    end
end

tElapsed = toc(tStart);               % calculation time
