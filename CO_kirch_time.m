                                                                     
                                                                     
                                                                     
                                             
 %{

CO_kirch.m - Constant offset Kirchhoff migration.
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

function [Kirchhoff, Skala] = CO_kirch(data, v, h, dt, dcmp, half_aper)

[nt,nx] = size(data);
Kirchhoff = zeros(nt,nx);
%% Loop over CMPs
for i_cmp=1:nx % Indices of neighbouring CMPs
    cmp = dcmp*i_cmp; % CMP-distance of the indices
    
%% Loop over timesamples
    for i_t=1:nt
        
        Tiefe = sqrt((2*h/(v)).^2+((i_t-1)*dt).^2); % traveltime depth
        
        t = sqrt(Tiefe^2 + (cmp/(v))^2); % TWT
        it = floor(1.5 + t/dt); % TWT index
        
        if(it > nt) % leave loop if out of Gather
            break;
        end
        
        amp = 4*(Tiefe*v)/(v^2*t); % Weightfunction
        % based on Zhang Y. (2000)
        
        % Aperture limits
        bound_l = max(floor(i_cmp-1), 1);
        bound_r = min(floor(nx-i_cmp), nx);
        if(i_cmp>(1+half_aper) && i_cmp<(nx-half_aper))
            bound_l = max(floor(i_cmp-1), -half_aper);
            bound_r = min(floor(nx-i_cmp), +half_aper);
        end
        
%% Loop over aperture
        for i_aper=bound_l:bound_r
            Kirchhoff(i_t,i_aper)=Kirchhoff(i_t,i_aper)+data(it,i_cmp+i_aper)*amp;
        end

    end
end
% TWT, daher * 0.5
Skala = sqrt((h/(v)).^2+((0:nt-1)'*dt).^2)*v*0.5; % Depth skaling
return 
