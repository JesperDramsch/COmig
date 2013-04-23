%{

CO_kirch_depth.m - Common offset Kirchhoff depth-migration.
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

function [COG] = CO_kirch_depth(data, v, h, dt, dz, dcmp, aper_half, flag_interp)


[nt,ns] = size(data);
t_orig=0:dt:((nt-1)*dt);
t_depth=t_orig*v*0.5;                   % TWT-time to depth conversion
z_max = max(t_depth);                    % Max depth [m]
z=0:dz:z_max;                            % Depthsampling
z_len = length(z);
COG(1:z_len,1:ns) = 0;                      % (depth, CMP)
cosphi = z./sqrt(z.^2+h^2);



%% Interpolation from t to z domain
% must be calculated before because aperture may need data in front
% of the current cmp
for i_cmp = 1:ns
    filt_interp(:,i_cmp) = interp1(t_depth,...
        data(:,i_cmp),z,'spline');
end

%% Loop over CMPs
for i_cmp = 1:ns
    % Aperture limits
    bound_l = max(floor(i_cmp-aper_half), 1);
    bound_r = min(floor(i_cmp+aper_half), ns);
    
    % Control if everything runs smoothly
    disp(['     CMP ||' ' left boundary ||'...
        ' right boundary ||' ' half aperture ||' ' velocity']);
    disp([(i_cmp-1)*dcmp (bound_l-1)*dcmp...
        (bound_r-1)*dcmp aper_half*dcmp v]);
    %% Loop over contributing samples (Aperture)
    for i_aper=bound_l:bound_r
        % Compute diffraction hyperbola, /2 because data is not TWT but depth
        z_diff = 0.5*( sqrt(z.^2 ...
            + ((i_cmp-i_aper)*dcmp-h).^2) ...
            + sqrt(z.^2 ...
            + ((i_cmp-i_aper)*dcmp+h).^2) );
        
        %% Compute amplitude correction
        weight = cosphi./sqrt(z_diff.*v);
        
        %% Loop over Depth
        for i_z=1:z_len
            
            % Exit if diffraction ist out of data
            if(z_diff(i_z) > z_max)
                break;
            end
            
            
            %% flag_interp zdiff
            % ! only if with interpolation at zdiff
            if(flag_interp==1)
                res_interp = interp1(z,...
                    filt_interp(:,i_aper),z_diff(i_z),'spline');
                
                %% Sum up along diffraction
                % ! with interpolation at zdiff
                COG(i_z,i_cmp) = COG(i_z,i_cmp) ...
                    + res_interp * weight(i_z);
            elseif(flag_interp==0)
                i_zdiff = floor(1.5+z_diff(i_z)/dz);
                % +0.5 so it get rounded correctly and + 1 so its
                % start with index 1, +1+0.5 = +1.5
                
                % ! without interpolation at zdiff
                COG(i_z,i_cmp) = COG(i_z,i_cmp) ...
                    + filt_interp(i_zdiff,i_aper) * weight(i_z);
            else
                disp('Error, no valid method!');
            end
        end
    end
end
COG(1,:) = 0;         % NaN avoiding
return
