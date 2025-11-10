%%% Program to Generate ftn58 Fromat from Lattice Tight Binding Model %%%
%%% Author: Hans                                                      %%%
clear all

%%-- Initial Conditions --%%
r0 = 2.7;                       %%% Nearest Neighbor hopping 4

Rso = 0.12;                      %%% intrinsic SOC                   (rank 2) -2
R1  = 0.0;                      %%% Rashba SOC with vertical E (break inversion)  NN  0
R2 = 0.0;                       %%% Rashba SOC which rotate spin direction          NNN 2

u  =  0;%1.165;                      %%% onsite energy of AB sublattice  (rank 4) 0.2
Hz =  0.0;                      %%% Zeeman field                    (rank 1) 0.3
Hx =  0.0;
Hy =  0.0;

isSO = 1; %  0: w/o spin || 1: w/ spin Degree of Freedom  %

%%% strange phase %%%
% r0 = 2.7;
% u  =  0.0;
% Hz =  0.2;
% Rso = 0.2;
% R1  = 0.2;
% R2 = 0.0;
%%%

%%% Renormalize factor %%%
sq3 = sqrt(3);
%R1 = 1i*R1*u;
R1 = 1i*R1;
Rso = Rso/(3*sqrt(3)); 
R2 = -1i*R2*2/3;

%%%--- Pauli Matrix ---%%%
sx = [0 1;1 0];
sy = [0 -1i;1i 0];
sz = [1 0;0 -1];

%%-- Generating Hopping Terms --%%

%%% graphene NN %%%
ftn58_r0 = [1 2 r0   0  0 0;
            1 2 r0  -1  0 0;
            1 2 r0   0 -1 0;];   

%ftn58 = [ftn58_r0];
%nbond = length(ftn58);
%ftn58 = [(1:nbond)' ftn58];
%ftn58 = [[2 nbond 0 0 0 0 0 ];ftn58];

if isSO == 1

%%% graphene NN %%%
ftn58_r0_spin = ftn58_r0;
ftn58_r0_spin(:,1:2) = ftn58_r0_spin(:,1:2) + 2;
ftn58_r0 = [ftn58_r0;ftn58_r0_spin];

%%% Rashba SOC : s_y*sig_x - s_x*sig_y %%%
% ftn58_R1 =  [1 4 -R1/(1i)*u   0  0 0;2 3  R1/(1i)*u   0  0 0; % - +(sx)/-(sy)
%              1 4  R1/(1i)*u  -1  0 0;2 3 -R1/(1i)*u  -1  0 0; % + -(sx)/+(sy)
%              1 4  (1i)*R1*u   0 -1 0;2 3  (1i)*R1*u   0 -1 0; %%% 1st term
%              1 4  R1/(1i)*u   0  0 0;2 3 -R1/(1i)*u   0  0 0;
%              1 4 -R1/(1i)*u  -1  0 0;2 3  R1/(1i)*u  -1  0 0;
%              1 4 -(1i)*R1*u   0 -1 0;2 3  (1i)*R1*u   0 -1 0;]; %%% 2nd term

ftn58_R1 =  [1 4  R1*(1+1i*sq3)/2   0  0 0;
             1 4  R1*(1-1i*sq3)/2  -1  0 0;
             1 4  R1*(-1)           0 -1 0;
             2 3  R1*(1+1i*sq3)/2   0  0 0;
             2 3  R1*(1-1i*sq3)/2  -1  0 0;
             2 3  R1*(-1)           0 -1 0;];

% ftn58_R1 =  [1 4 -R1/(1i)   0  0 0;2 3  R1/(1i)   0  0 0; % - +(sx)/-(sy)
%              1 4  R1/(1i)  -1  0 0;2 3 -R1/(1i)  -1  0 0; % + -(sx)/+(sy)
%              1 4  (1i)*R1   0 -1 0;2 3  (1i)*R1   0 -1 0; %%% 1st term
%              1 4  R1/(1i)   0  0 0;2 3 -R1/(1i)   0  0 0;
%              1 4 -R1/(1i)  -1  0 0;2 3  R1/(1i)  -1  0 0;
%              1 4 -(1i)*R1   0 -1 0;2 3  (1i)*R1   0 -1 0;]; %%% 2nd term

ftn58_R2 =  [1 3   1i*R2             1  0 0;    %tR_a1 = ()
             1 3   (sq3+(1i))*R2/2   0  1 0;    %H(1,3)=tR_a1*e^(ik.a2);
             1 3   (sq3-(1i))*R2/2  -1  1 0;
             1 3  -1i*R2            -1  0 0;
             1 3   (-sq3-(1i))*R2/2  0 -1 0;
             1 3   (-sq3+(1i))*R2/2  1 -1 0;
             2 4  -1i*R2             1  0 0;
             2 4  -(sq3+(1i))*R2/2   0  1 0;
             2 4  -(sq3-(1i))*R2/2  -1  1 0;
             2 4   1i*R2            -1  0 0;
             2 4   (sq3+(1i))*R2/2   0 -1 0;
             2 4   (sq3-(1i))*R2/2   1 -1 0;];

% ftn58_R2 =  [1 3 -1i*R2       1 -1 0;
%              1 3 -1i*R2       0 -1 0;
%              1 3 -1i*R2      -1  0 0;
%              1 3 -1i*R2*(-1) -1  1 0;
%              1 3 -1i*R2*(-1)  0  1 0;
%              1 3 -1i*R2*(-1)  1  0 0;
%              2 4  1i*R2       1 -1 0;
%              2 4  1i*R2       0 -1 0;
%              2 4  1i*R2      -1  0 0;
%              2 4  1i*R2*(-1) -1  1 0;
%              2 4  1i*R2*(-1)  0  1 0;
%              2 4  1i*R2*(-1)  1  0 0;];

ftn58_Rso = [1 1   Rso*1i  1 -1 0;
             1 1  -Rso*1i  0 -1 0;
             1 1   Rso*1i -1  0 0;
             1 1  -Rso*1i -1  1 0;
             1 1   Rso*1i  0  1 0;
             1 1  -Rso*1i  1  0 0;
             2 2  -Rso*1i  1 -1 0;
             2 2   Rso*1i  0 -1 0;
             2 2  -Rso*1i -1  0 0;
             2 2   Rso*1i -1  1 0;
             2 2  -Rso*1i  0  1 0;
             2 2   Rso*1i  1  0 0;
             3 3  -Rso*1i  1 -1 0;
             3 3   Rso*1i  0 -1 0;
             3 3  -Rso*1i -1  0 0;
             3 3   Rso*1i -1  1 0;
             3 3  -Rso*1i  0  1 0;
             3 3   Rso*1i  1  0 0;
             4 4   Rso*1i  1 -1 0;
             4 4  -Rso*1i  0 -1 0;
             4 4   Rso*1i -1  0 0;
             4 4  -Rso*1i -1  1 0;
             4 4   Rso*1i  0  1 0;
             4 4  -Rso*1i  1  0 0;];

% ftn58_Rso = [1 1  -Rso*1i  1 -1 0;
%              1 1   Rso*1i  0 -1 0;
%              1 1  -Rso*1i -1  0 0;
%              1 1   Rso*1i -1  1 0;
%              1 1  -Rso*1i  0  1 0;
%              1 1   Rso*1i  1  0 0;
%              2 2   Rso*1i  1 -1 0;
%              2 2  -Rso*1i  0 -1 0;
%              2 2   Rso*1i -1  0 0;
%              2 2  -Rso*1i -1  1 0;
%              2 2   Rso*1i  0  1 0;
%              2 2  -Rso*1i  1  0 0;
%              3 3   -Rso*1i  1 -1 0;
%              3 3  Rso*1i  0 -1 0;
%              3 3   -Rso*1i -1  0 0;
%              3 3  Rso*1i -1  1 0;
%              3 3   -Rso*1i  0  1 0;
%              3 3  Rso*1i  1  0 0;
%              4 4  Rso*1i  1 -1 0;
%              4 4   -Rso*1i  0 -1 0;
%              4 4  Rso*1i -1  0 0;
%              4 4   -Rso*1i -1  1 0;
%              4 4  Rso*1i  0  1 0;
%              4 4   -Rso*1i  1  0 0;];



% ftn58_VZ = [1 1 VZ 1 0 0; 1 1 VZ -1 0 0; 1 1 VZ 0 1 0; 1 1 VZ 0 -1 0;
%             2 2 VZ 1 0 0; 2 2 VZ -1 0 0; 2 2 VZ 0 1 0; 2 2 VZ 0 -1 0;
%             3 3 -VZ 1 0 0; 3 3 -VZ -1 0 0; 3 3 -VZ 0 1 0; 3 3 -VZ 0 -1 0;
%             4 4 -VZ 1 0 0; 4 4 -VZ -1 0 0; 4 4 -VZ 0 1 0; 4 4 -VZ 0 -1 0;];

%%% Zeeman Field Term %%%
ftn58_Mz = [1 1  Hz 0 0 0; 2 2  Hz 0 0 0;
            3 3 -Hz 0 0 0; 4 4 -Hz 0 0 0;];
% ftn58_Mz = [1 3  Hz 0 0 0; 2 4  Hz 0 0 0;
%             3 1  Hz 0 0 0; 4 2  Hz 0 0 0;];

ftn58_Mx = [1 3  Hx 0 0 0; 2 4  Hx 0 0 0;
            3 1  Hx 0 0 0; 4 2  Hx 0 0 0;];

ftn58_My = [1 3  -1i*Hy 0 0 0; 2 4  -1i*Hy 0 0 0;
            3 1   1i*Hy 0 0 0; 4 2   1i*Hy 0 0 0;];

%%% Onsite Term %%%
ftn58_On = [1 1 u 0 0 0; 2 2 -u 0 0 0;
            3 3 u 0 0 0; 4 4 -u 0 0 0;];

ftn58 = [ftn58_r0;ftn58_R1;ftn58_R2;ftn58_Rso;ftn58_Mz;ftn58_On;ftn58_Mx;ftn58_My];
% ftn58 = [ftn58_r0;ftn58_Rso;ftn58_Rso2;ftn58_Mz;ftn58_On];
%ftn58 = [ftn58_r0;ftn58_R1;ftn58_R2;ftn58_Rso;ftn58_Mz;ftn58_On;ftn58_VZ];
nbond = length(ftn58);
ftn58 = [(1:nbond)' ftn58];
ftn58 = [[4 nbond 0 0 0 0 0 ];ftn58];
end
             
save ftn58.mat ftn58
%%

eval(['run TBHmftn.m'])
%eval(['run sparse2ftn58full.m'])
eval(['run band_ftn_pap.m'])
% eval(['run band_ftn_pap.m'])
%eval(['run ./Berry-curvature-plot/Chern_main'])
