
%%% Program to Generate ftn58 Fromat from Lattice Tight Binding Model %%%
%%% Author: Hans                                                      %%%
clear all

%%-- Initial Conditions --%%
tt =  -1;  %2.7; %%% t1  %%%
phi = pi/2;
t2 = 0.1;  % 0.0, 0.1, 0.2
u  = 0.;  % 1

pp1 = t2*exp(1i*phi);  %% t2 %%%
pp2 = t2*exp(-1i*phi);

save titletmp.mat *
%%%--- Pauli Matrix ---%%%
sx = [0 1;1 0];
sy = [0 -1i;1i 0];
sz = [1 0;0 -1];

%%-- Generating Hopping Terms --%%
%%% Kinetic Term %%%
ftn58_KT = [1 2 tt  0 0 0; 1 2 tt 0 -1 0;
            1 2 tt -1 0 0; 2 1 tt 0  0 0;
            2 1 tt  0 1 0; 2 1 tt 1  0 0;
            1 1 pp2  0  1 0; 1 1  pp1  0 -1 0;
            1 1 pp2  1 -1 0; 1 1  pp1 -1  1 0; 
            1 1 pp2 -1  0 0; 1 1  pp1  1  0 0;
            2 2 pp1  0  1 0; 2 2  pp2  0 -1 0;
            2 2 pp1  1 -1 0; 2 2  pp2 -1  1 0
            2 2 pp1 -1  0 0; 2 2  pp2  1  0 0;];
         
%%% Onsite Term %%%
ftn58_On = [1 1 +u 0 0 0; 2 2 -u 0 0 0];

ftn58 = [ftn58_KT;ftn58_On];
nbond = length(ftn58);
ftn58 = [(1:nbond)' ftn58];
ftn58 = [[2 nbond 0 0 0 0 0 ];ftn58];

save ftn58.mat ftn58
%%

% eval(['run sparse2ftn58full.m']);
eval(['run TBHmftn.m']);
eval(['run band_ftn.m']);
