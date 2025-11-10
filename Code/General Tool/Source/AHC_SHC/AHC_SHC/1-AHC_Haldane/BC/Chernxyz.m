%%% Functon for calculating Spin Hall conductivity via Kubo formula (3D) %%%
%%% -------------------------------------------------------------------- %%%
%%% Define spin Berry curvature Omega^k_(i,j), where i and j stands for  %%%
%%% the direction of spin current and electric field. And k is for spin  %%%
%%% direction (k=x,y,z).
%%% Example: Sig^{Sk}_{Ji,Ej}
% Original Author : Chuang-Han Hsu
%%% -------------------------------------------------------------------- %%%
function [Berry,Ek] = Chernxyz(kpoints,ftn58sparse)

%% --- Initial setup for QSHC calcultaion --- %%%
BR   = ftn58sparse.BR;
abc  = ftn58sparse.abc; 
Norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
tt   = ftn58sparse.tt;
dd   = ftn58sparse.dd;
gam  = 1e-6;

T  = [BR(:,1)*abc(1) BR(:,2)*abc(2) BR(:,3)*abc(3)]; 
DD = dd*T;
DX = DD(:,1);
DY = DD(:,2);
DZ = DD(:,3);

%% --- Actual Procedure --- %%%
Berry   = zeros(Norb,1);  
LM      = zeros(Norb,1);
Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,Norb,Norb);
H0      = full(Hsparse);
H0      = (H0 + H0')/2;
[V,D]   = eig(H0);
Ek      = diag(D);   
            
vx = 1i*full(sparse(ii,jj,DX.*exp(1i*dd*kpoints').*tt,Norb,Norb));    % dH/dkx
vy = 1i*full(sparse(ii,jj,DY.*exp(1i*dd*kpoints').*tt,Norb,Norb));    % dH/dky
vz = 1i*full(sparse(ii,jj,DZ.*exp(1i*dd*kpoints').*tt,Norb,Norb));    % dH/dkz
vx = (vx + vx')/2;
vy = (vy + vy')/2;
vz = (vz + vz')/2;

%%%  Matrix element for <u|v|u>  %%%
V_x = V'*vx*V; 
V_x = (V_x + V_x')/2;  % matrix elements for velocity for E field
V_y = V'*vy*V; 
V_y = (V_y + V_y')/2;  % matrix elements for velocity for E field
V_z = V'*vz*V; 
V_z = (V_z + V_z')/2;  % matrix elements for velocity for E field

ch_xx = zeros(Norb,Norb);
ch_xy = zeros(Norb,Norb); 
ch_xz = zeros(Norb,Norb); 
ch_yx = zeros(Norb,Norb);
ch_yy = zeros(Norb,Norb); 
ch_yz = zeros(Norb,Norb); 
ch_zx = zeros(Norb,Norb);
ch_zy = zeros(Norb,Norb); 
ch_zz = zeros(Norb,Norb); 
 
% lm_yz = zeros(Norb,Norb);
% lm_zx = zeros(Norb,Norb);
% lm_xy = zeros(Norb,Norb);
            
for nb1 = 1:Norb  % no contribution for nb2=nb1 & time-reversal(nb1)
    for nb2 = nb1+1:Norb     % half-part of the matrix
        diff_energy = Ek(nb1) - Ek(nb2); 

        commom_term_Vxp = V_x(nb2,nb1)/(diff_energy^2+1i*gam);
        commom_term_Vyp = V_y(nb2,nb1)/(diff_energy^2+1i*gam);
        commom_term_Vzp = V_z(nb2,nb1)/(diff_energy^2+1i*gam);
        commom_term_Vxn = V_x(nb2,nb1)/(diff_energy^2-1i*gam);
        commom_term_Vyn = V_y(nb2,nb1)/(diff_energy^2-1i*gam);
        commom_term_Vzn = V_z(nb2,nb1)/(diff_energy^2-1i*gam);
             
%         if abs(diff_energy) <= 1e-10
%             commom_term_Vx = 0.0;
%             commom_term_Vy = 0.0;
%             commom_term_Vz = 0.0;
%         else
%             commom_term_Vx = V_x(nb2,nb1)*diff_energy^-2;
%             commom_term_Vy = V_y(nb2,nb1)*diff_energy^-2;
%             commom_term_Vz = V_z(nb2,nb1)*diff_energy^-2;
%         end

        ch_xx(nb1,nb2) = imag(V_x(nb1,nb2)*commom_term_Vxp + V_x(nb1,nb2)*commom_term_Vxn);  
        ch_xy(nb1,nb2) = imag(V_x(nb1,nb2)*commom_term_Vyp + V_x(nb1,nb2)*commom_term_Vyn);
        ch_xz(nb1,nb2) = imag(V_x(nb1,nb2)*commom_term_Vzp + V_x(nb1,nb2)*commom_term_Vzn);  
        ch_yx(nb1,nb2) = imag(V_y(nb1,nb2)*commom_term_Vxp + V_y(nb1,nb2)*commom_term_Vxn);
        ch_yy(nb1,nb2) = imag(V_y(nb1,nb2)*commom_term_Vyp + V_y(nb1,nb2)*commom_term_Vyn);  
        ch_yz(nb1,nb2) = imag(V_y(nb1,nb2)*commom_term_Vzp + V_y(nb1,nb2)*commom_term_Vzn);
        ch_zx(nb1,nb2) = imag(V_z(nb1,nb2)*commom_term_Vxp + V_z(nb1,nb2)*commom_term_Vxn);  
        ch_zy(nb1,nb2) = imag(V_z(nb1,nb2)*commom_term_Vyp + V_z(nb1,nb2)*commom_term_Vyn);
        ch_zz(nb1,nb2) = imag(V_z(nb1,nb2)*commom_term_Vzp + V_z(nb1,nb2)*commom_term_Vzn);  
        
%         ch_yx(nb1,nb2) = imag(V_y(nb1,nb2)*commom_term_Vx);          
%         ch_zy(nb1,nb2) = imag(V_z(nb1,nb2)*commom_term_Vy);  
%         ch_xz(nb1,nb2) = imag(V_x(nb1,nb2)*commom_term_Vz);   
    end  
end

ch_xx = ch_xx - ch_xx.';
ch_xy = ch_xy - ch_xy.';
ch_xz = ch_xz - ch_xz.';
ch_yx = ch_yx - ch_yx.';
ch_yy = ch_yy - ch_yy.';
ch_yz = ch_yz - ch_yz.';
ch_zx = ch_zx - ch_zx.';
ch_zy = ch_zy - ch_zy.';
ch_zz = ch_zz - ch_zz.';


Berry(:,1) = sum(ch_xx,2);Berry(:,4) = sum(ch_yx,2);Berry(:,7) = sum(ch_zx,2);
Berry(:,2) = sum(ch_xy,2);Berry(:,5) = sum(ch_yy,2);Berry(:,8) = sum(ch_zy,2);
Berry(:,3) = sum(ch_xz,2);Berry(:,6) = sum(ch_yz,2);Berry(:,9) = sum(ch_zz,2);

