function mig_graphs(data,dcmp,dt,nh,flag)

%% Farben definieren
farbe = rand(nh,3);
farbe(:,1)=sort(farbe(:,1));
farbe(:,2)=sort(farbe(:,2),'Descend');

%% Input data plots
% Compare max amplitudes of each CMP in original data
figure
hold on
for k=1:nh
plot(((1:ns)-1)*dcmp,max(abs(data(:,:,k)),[],1),'Color',farbe(k,:))
end
xlabel('CMP','Fontsize',24)
ylabel('Maximum amplitude','Fontsize',24)
legend('h = 0','h = 250','h = 500','h = 750','h = 1000','Location','best')
set(gca,'Fontsize',24)
print('-dpng',sprintf('output/ampsorig.png'));


% Compare max amplitudes of each CMP in original data
figure
hold on
for k=1:nh
plot(((1:ns)-1)*dcmp,max(abs(filtdata(:,:,k)),[],1),'Color',farbe(k,:))
end
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

%% Frequency Plot
fx = figure(1);
plot(faxis,abs(fdata(1:length(faxis)))/max(abs(fdata(1:length(faxis)))));
xlabel('Frequency','Fontsize',24);
ylabel('Normalized Amplitude','Fontsize',24);
%title('Frequency analysis','Fontsize',24)
set(gca,'Fontsize',24)
set(fx, 'Position', [0 0 1280 1024] );    % Size of the new frame
axis ([0 75 0 1])
print('-dpng','output/freq.png');         % Outputfile of figure

%Migration result
fy = figure;
plot(faxis,abs(migdata(1:length(faxis)))/max(abs(migdata(1:length(faxis)))));
xlabel('Frequency [Hz]','Fontsize',24);
ylabel('Normalized Amplitude','Fontsize',24);
set(gca,'Fontsize',24)
set(fy, 'Position', [0 0 1280 1024] );    % Size of the new frame
axis ([0 75 0 1])
print('-dpng','freq_mig.png');                % Outputfile of figure

% Frequ comp
figure
plot(faxis,abs(fdata(1:length(faxis)))/max(abs(fdata(1:length(faxis)))),'k');
hold on
plot(faxis,abs(migdata(1:length(faxis)))/max(abs(migdata(1:length(faxis)))),'r');
xlabel('Frequency [Hz]','Fontsize',24);
ylabel('Normalized Amplitude','Fontsize',24);
legend('Original Data','Migrated Data')
set(gca,'Fontsize',24)
set(fx, 'Position', [0 0 1280 1024] );    % Size of the new frame
axis ([0 75 0 1])
print('-dpng','freq_comp.png');                % Outputfile of figure


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
hold on
for k=1:nh
plot(((1:ns)-1)*dcmp,max(abs(COG(:,:,k)),[],1),'Color',farbe(k,:))
end
xlabel('CMP','Fontsize',24)
ylabel('Maximum amplitude','Fontsize',24)
legend('h = 0','h = 250','h = 500','h = 750','h = 1000','Location','best')
set(gca,'Fontsize',24)
print('-dpng',sprintf('output/ampsv%g.png',v));

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

%% SNR
% Plot trace 51 normalized
figure
plot(((1:nt)-1)*dt,filtdata(:,51,1)/max(filtdata(:,51,1)),'r')
hold on
plot(((1:nt)-1)*dt,mig(:,51)/max(mig(:,51)),'k')
ylabel('Normalisierte Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('SNR Input','SNR Migriert','Location','NorthWest')
set(gca,'Fontsize',24)
print('-dpng','SNRnorm.png');

% Plot trace 51 not normalized
figure
plot(((1:nt)-1)*dt,filtdata(:,51,1),'r')
hold on
plot(((1:nt)-1)*dt,mig(:,51),'k')
ylabel('Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('SNR Input','SNR Migriert','Location','NorthWest')
set(gca,'Fontsize',24)
print('-dpng','SNRreal.png');


%% Ergebnis Plots
%Input signal normalized
figure
plot(((1:nt)-1)*dt,filtdata(:,51,1)/max(filtdata(:,51,1)),'r')
hold on
plot(((1:nt)-1)*dt,data(:,51,1)/max(data(:,51,1)),'k')
ylabel('Normalisierte Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('Filtered data','Original data','Location','NorthWest')
set(gca,'Fontsize',24)
print('-dpng','inoutNorm.png');

%Input signal not normalized
figure
plot(((1:nt)-1)*dt,filtdata(:,51,1),'r')
hold on
plot(((1:nt)-1)*dt,data(:,51,1),'k')
ylabel('Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('Filtered data','Original data','Location','NorthWest')
set(gca,'Fontsize',24)
print('-dpng','inout.png');

%Comparison results normalized
figure
plot(((1:nt)-1)*dt,mig(:,51,1)/max(mig(:,51,1)),'r')
hold on
plot(((1:nt)-1)*dt,data(:,51,1)/max(data(:,51,1)),'k')
ylabel('Normalisierte Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('Migration result','Original data','Location','NorthWest')
set(gca,'Fontsize',24)
print('-dpng','waveletNorm.png');

%Comparison results not normalized
figure
plot(((1:nt)-1)*dt,mig(:,51,1),'r')
hold on
plot(((1:nt)-1)*dt,data(:,51,1),'k')
ylabel('Amplitude','Fontsize',24)
xlabel('Zeit [s]','Fontsize',24)
legend('Migration result','Original data','Location','NorthWest')
set(gca,'Fontsize',24)
print('-dpng','wavelet.png');

% Plot of the migration result
fx=figure(v+1);
set(fx, 'Position', [0 0 1280 1024] );
imagesc(((1:ns)-1)*dcmp/1000,Skala(:,1)/1000,mig(:,:),[-max(max(abs(mig))) max(max(abs(mig)))])
%title('Tiefenmigration','Fontsize',24)
xlabel('CMP [km]','Fontsize',24)
ylabel('Depth [km]','Fontsize',24)
colormap([ones(101,1),(0:.01:1)',(0:.01:1)';(1:-.01:0)',(1:-.01:0)',ones(101,1)])
colorbar
set(gca,'Fontsize',24)
print('-dpng',sprintf('sum_v%g.png',v));