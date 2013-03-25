function [Kirchhoff, Skala] = CO_kirch(data, v, phi, h, dt, dcmp)

%{

    CO_kirch.m - Constant offset Kirchhoff migration.
    Copyright (C) 2013  Jesper S Dramsch, Matthias Schneider, Jan Walda

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

[nt,nx] = size(data);
Kirchhoff = zeros(nt,nx);

%% Schleife ueber CMPs
for i_cmp=-nx:nx
    cmp = dcmp*i_cmp;
    
%% Schleife ueber Samples
    for i_t=1:nt
        
        Tiefe = sqrt((h/(v)).^2+((i_t-1)*dt).^2);            % Laufzeittiefe
        
        t = sqrt(Tiefe^2 + (cmp/(v))^2);                     % TWT
        it = floor(1.5 + t/dt);                                % TWT Index
        
        if(it > nt)                                            % Abbruchkriterium
            break;
        end
        
        amp    = cos(phi)*Tiefe*sqrt(nt*dt/(2*pi*v*t));      % Gewichtsfunktion
        %amp    =  cos(phi)* (i_t-1)*dt * sqrt( nt * dt / (2  * v^2 * pi) );

        bound_l = max( floor(1-i_cmp),  1);                    % Apertur
        bound_r = min(floor(nx-i_cmp), nx);
        
%% Schleife ueber Apertur        
        for i_aper=bound_l:bound_r
            Kirchhoff(i_t,i_aper)=Kirchhoff(i_t,i_aper)+data(it,i_aper+i_cmp)*amp;
        end
        
    end
end
Skala = sqrt((h/(v)).^2+((0:nt-1)'*dt).^2)*v;                                       % Tiefenskalierung
return
