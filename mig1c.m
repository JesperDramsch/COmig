%{

    mig1c.m - Analysis of Velocity errors in migration routines
    Copyright (C) 2013  Jesper S Dramsch, Matthias Schneider, Dela Spickermann, Jan Walda

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




clear all
close all
clc

z0 = 2000;                          %Wahre angenommene Reflektortiefe      
v_verhaeltnis   = [0.9;1;1.1].^2;   %Vm/V0
hz0 = 0:.001:3;                     %h/z0

vvvektor = repmat([0.9;1;1.1].^2,size(hz0));
hz0vektor = repmat(hz0,size(v_verhaeltnis)).^2;
z=z0*sqrt(vvvektor+hz0vektor.*(vvvektor-1));

plot(hz0,real(z'),'LineWidth',2)
set(gca,'YDir','reverse','FontSize',24);
legend('V_M = 0.9\cdotV_0','V_M = V_0','V_M = 1.1\cdotV_0','FontSize',24,'Orientation','horizontal','Location','SouthOutside')
title('Aufgabe 1c: Geschwindigkeitsfehler bei horizontalem Reflektor','FontSize',24)
xlabel('h/z_0','FontSize',24)
ylabel('Tiefenlage des Reflektors [m]','FontSize',24)
