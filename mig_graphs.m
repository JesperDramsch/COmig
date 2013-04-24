function mig_graphs(switcher,varargin)
%{

mig_graphs.m - Plot function for Kirchhoff migration routine.
Copyright (C) 2013 Jesper S Dramsch, Matthias Schneider, Dela Spickermann, Jan Walda

 Use:
 SingleLine for Lineplot, i.e. Frequency plots
 CompLine for Two Lineplots, i.e. Trace comparison
 OffsetLine for one Lineplot per Offset, i.e. Amplitudechecks of different Offsets
 COG for Constant Offset Plots
 PolarPlot for seismic Polarity Plot, i.e. Summation of Migration result
 For templates see individual Plot option


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

switch switcher
    case 'SingleLine'
        % Template: mig_graphs('SingleLine',data,axis_vector,xlabel,ylabel,fileoutput)
        figure
        plot(varargin{2},varargin{1},'k');
        xlabel(varargin{3},'Fontsize',24);
        ylabel(varargin{4},'Fontsize',24);
        
    case 'CompLine'
        % Template: mig_graphs('CompLine',name1,input1,name2,input2,axis_vector,xaxis_vector,yaxis_vector,output)
        figure
        hold on
        plot(varargin{5},varargin{2},'k')
        plot(varargin{5},varargin{4},'r')
        ylabel(varargin{7},'Fontsize',24)
        xlabel(varargin{6},'Fontsize',24)
        legend(varargin{1},varargin{3},'Location','NorthEast')
        
    case 'OffsetLine'
        % Template: mig_graphs('OffsetLine',data,axis_vector,name)
        data = varargin{1};
        [nt,ns,nh] = size(data);
        farbe = rand(nh,3);
        farbe(:,1)=sort(farbe(:,1));
        farbe(:,2)=sort(farbe(:,2),'Descend');
        
        figure
        hold on
        for k=1:nh
            plot(varargin{2},max(abs(data(:,:,k)),[],1),'Color',farbe(k,:))
        end
        xlabel('CMP','Fontsize',24)
        ylabel('Maximum amplitude','Fontsize',24)
        legend('h = 0','h = 250','h = 500','h = 750','h = 1000','Location','best')
        
    case 'COG'
        % Template: mig_graphs('COG',data,CMPaxis_vector,yaxis_vector,ylabel,name)
        figure
        imagesc(varargin{2},varargin{3},varargin{1},[-max(max(max(abs(varargin{1})))) +max(max(max(abs(varargin{1}))))])
        xlabel('CMP [km]','Fontsize',24)
        ylabel(varargin{4},'Fontsize',24)
        colormap([ones(101,1),(0:.01:1)',(0:.01:1)';(1:-.01:0)',(1:-.01:0)',ones(101,1)])    % polarized plot
        colorbar
        set(gca,'XTick',(0:5)*2e3)                  
        set(gca,'XTickLabel',{' 0 ','2 / 0','2 / 0','2 / 0','2 / 0',' 2 '})                  % reskaling x-axis
        
        
    case 'PolarPlot'
        % Template: mig_graphs('PolarPlot',data,yaxis_vector,xaxis_vector,xlabel,ylabel,name)
        figure
        imagesc(varargin{2},varargin{3},varargin{1},[-max(max(max(abs(varargin{1})))) +max(max(max(abs(varargin{1}))))])
        %title('Tiefenmigration','Fontsize',24)
        xlabel(varargin{4},'Fontsize',24)
        ylabel(varargin{5},'Fontsize',24)
        colormap([ones(101,1),(0:.01:1)',(0:.01:1)';(1:-.01:0)',(1:-.01:0)',ones(101,1)])
        colorbar
        
        
    otherwise
        warndlg('Please choose valid Graph option')
end
%%% Für alte traceplots müsste zurück interpoliert werden, da
%%% length(z) != length(t) und length(z) hängt von v ab, siehe
%%% t_depth=t_orig*v*0.5;
%%% zmax = max(t_depth);     % zmax nimmt mit steigendem v zu, aber
%%% z=0:dz:zmax;             % dz bleibt gleich !

set(gca,'Fontsize',24)
print('-dpng',sprintf('%s/output/%s.png',pwd,varargin{end}));
