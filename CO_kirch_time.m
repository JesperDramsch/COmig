%{

CO_kirch_time.m - Common offset Kirchhoff migration.
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

function [COG, z] = CO_kirch_time(data, v, h, dt, dcmp, aper_half, flag_interp)

%% Open arrays and variables
[nt,ns] = size(data);
t_orig=(0:dt:((nt-1)*dt))';
t_max = max(t_orig);
COG(1:nt,1:ns) = 0;

%% Loop over CMPs
for i_cmp=1:ns % Indices of neighbouring CMPs
    % Aperture limits
    bound_l = max(floor(i_cmp-aper_half), 1);
    bound_r = min(floor(i_cmp+aper_half), ns);
    
    % Possible control if everything runs smoothly
    %{
    disp(['     CMP ||' ' left boundary ||'...
        ' right boundary ||' ' half aperture ||' ' velocity']);
    disp([(i_cmp-1)*dcmp (bound_l-1)*dcmp...
        (bound_r-1)*dcmp aper_half*dcmp v]);
    %}
    fprintf('||\tCMP \t||\tleft boundary \t||\tright boundary \t||\thalf aperture \t||\tvelocity\t||\n||\t%4g \t||\t\t%4g\t\t||\t\t%4g\t\t||\t\t%4g\t\t||\t%4g\t\t||\n',(i_cmp-1)*dcmp,(bound_l-1)*dcmp,...
        (bound_r-1)*dcmp,aper_half*dcmp,v)
    %% Loop over contributing samples (Aperture)
    for i_aper=bound_l:bound_r
        % Compute diffraction hyperbola, /2 because data is not TWT but depth
        t_diff = 0.5* ( sqrt((t_orig*v).^2 ...
            + ((i_cmp-i_aper)*dcmp-h).^2) ...
            + sqrt((t_orig*v).^2 ...
            + ((i_cmp-i_aper)*dcmp+h).^2) )/v;
        
        % Exit if diffraction ist out of data
        t_flag = (t_diff - t_max <= 0);
        
        %% Compute amplitude correction
        weight = 4*(t_diff*v)./(v^2.*t_orig);
        % based on Zhang Y. (2000)
        
        %% flag_interp zdiff
        % ! only if with interpolation at tdiff
        if(flag_interp==1)
            res_interp = interp1(t_orig,data(:,i_aper),t_diff,'spline');
            
            % Sinc approach leaves "Sinc-reverb" Looks funny, try it :-)
            %res_interp = sinc(t_diff(:,ones(size(t_orig))) - t_orig(:,ones(size(t_diff)))')*data(:,i_aper);

            %% Sum up along diffraction
            % ! with interpolation at tdiff
            COG(:,i_cmp) = COG(:,i_cmp) ...
                + res_interp .* weight .* t_flag;
        elseif(flag_interp==0)
            i_tdiff = ((floor(1.5+t_diff./dt)-1).*t_flag)+1;
            % +0.5 so it get rounded correctly and + 1 so its
            % start with index 1, +1+0.5 = +1.5
            % ! without interpolation at zdiff
            COG(:,i_cmp) = COG(:,i_cmp) ...
                + data(i_tdiff,i_aper) .* weight .* t_flag;
        else
            disp('Error, no valid method!');
        end
    end
end
% TWT, therefore * 0.5
COG(1,:) = 0;         % NaN avoiding

z = sqrt((h/(v)).^2+((0:nt-1)'*dt).^2)*v*0.5; % Depth skaling

return
