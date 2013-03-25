function [Kirchhoff, Skala] = CO_kirch(data, v, phi, h, dt, dcmp)

[nt,nx] = size(data);
Kirchhoff = zeros(nt,nx);

%% Schleife ueber CMPs
for i_cmp=-nx:nx
    cmp = dcmp*i_cmp;
    
%% Schleife ueber Samples
    for i_t=1:nt
        
        Tiefe = sqrt((h/(2*v)).^2+((i_t-1)*dt).^2);            % Laufzeittiefe
        
        t = sqrt(Tiefe^2 + (cmp/(2*v))^2);                     % TWT
        it = floor(1.5 + t/dt);                                % TWT Index
        
        if(it > nt)                                            % Abbruchkriterium
            break;
        end
        
        amp    = cos(phi)*Tiefe*sqrt(nt*dt/(2*pi*v*t))/t;      % Gewichtsfunktion
        
        bound_l = max( floor(1-i_cmp),  1);                    % Apertur
        bound_r = min(floor(nx-i_cmp), nx);
        
%% Schleife ueber Apertur        
        for i_aper=bound_l:bound_r
            Kirchhoff(i_t,i_aper)=Kirchhoff(i_t,i_aper)+data(it,i_aper+i_cmp)*amp;
        end
        
    end
end
Skala = sqrt((h/(2*v)).^2+((0:nt-1)'*dt).^2)*v;                                       % Tiefenskalierung
return
