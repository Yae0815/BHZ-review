%%% Program for calculate Spin Hall conductivity via Kubo formula (3D) %%%
%%% ------------------------------------------------------------------ %%%
%%% Define Berry curvature Omega^k_(i,j), where i and j stands for     %%%
%%% the direction of current and electric field.                       %%%                              %%%
%%% Example: \sigma_{Ji,Ej}                                            %%%
%%% ------------------------------------------------------------------ %%%
clear all

cquanta = 1;            % in units of (hbar/e) Ohm^-1 cm^-1
kB      = 8.6173324E-5; % in (eV*K^-1)

%% Initialization %%%
wcal = ReadInput('input.txt');
Nk1  = wcal.Nk1;
Nk2  = wcal.Nk2;
% allow 3D meshes (fallback to Nk2 if Nk3 missing)
if isfield(wcal,'Nk3') && ~isempty(wcal.Nk3)
    Nk3 = wcal.Nk3;
else
    Nk3 = Nk2; % minimal-change default for cubic grids
end
Temp = wcal.Temp;
Ef   = wcal.Ef;
nE   = wcal.nE; 
np   = wcal.ncpu;

g1   = wcal.g1(2:end-1);
sid  = find(isspace(g1));
G1   = [str2double(g1(1:sid(1)-1)) str2double(g1(sid(1)+1:end))];
 
g2   = wcal.g2(2:end-1);
sid  = find(isspace(g2));
G2   = [str2double(g2(1:sid(1)-1)) str2double(g2(sid(1)+1:end))];

% g3 may be absent in older inputs; default to full reduced span
if isfield(wcal,'g3') && ~isempty(wcal.g3)
    g3   = wcal.g3(2:end-1);
    sid  = find(isspace(g3));
    G3   = [str2double(g3(1:sid(1)-1)) str2double(g3(sid(1)+1:end))];
else
    G3   = [0 1];
end
 
em   = wcal.Emu(2:end-1);
sid  = find(isspace(em));
Emu  = [str2double(em(1:sid(1)-1)) str2double(em(sid(1)+1:end))];

%% --- Load TB-model --- %%s
instruct = load(wcal.file);
if isfield(instruct,'ftn58sparse')
    ftn58sparse = instruct.ftn58sparse;
elseif isfield(instruct,'SPftn58sparse')
    ftn58sparse = instruct.SPftn58sparse;
else
    ftn58sparse = instruct.Sftn58sparse;
    ftn58sparse.dd = [ftn58sparse.dd,zeros(size(ftn58sparse.dd,1),1)];
    ftn58sprase.BR = ftn58sparse.BR2D;
end
%% Transform to Wannier represetation
ftn58sparse = ftn58sparse_transform_representation(ftn58sparse,1);

%% K-mesh setup (3D BZ for cubic system) %%%
k1      = linspace(G1(1),G1(2),Nk1);
k2      = linspace(G2(1),G2(2),Nk2);
k3      = linspace(G3(1),G3(2),Nk3);
[K1,K2,K3] = ndgrid(k1,k2,k3);
kpts    = [K1(:), K2(:), K3(:)];
nks     = size(kpts,1);

%% In CARTESIAN coordinates
BR      = ftn58sparse.BR;
abc     = ftn58sparse.abc;
T       = [BR(:,1)*abc(1) BR(:,2)*abc(2) BR(:,3)*abc(3)];
Ta      = norm(T(1,:));
Tb      = norm(T(2,:));
Tc      = norm(T(3,:));

%% --- Open the parallel pool --- %%%
% c = parcluster('local');
% c.NumWorkers = np;
% parpool(c, c.NumWorkers);

%% --- Spin Curvature calculation --- %%%
norb     = ftn58sparse.norb;
SpinCurv = zeros(nks,norb,27);
Ek       = zeros(nks,norb);
fprintf('Start calculating curvature tensors ... \n')
tic
parfor nk=1:nks
        kvec             = kpts(nk,:)*2*pi;
        [Stmp,Etmp]      = SpinCurvature(kvec,ftn58sparse);
        SpinCurv(nk,:,:) = Stmp;
        Ek(nk,:)         = Etmp;               
end
toc 
save SpinCurvatureTexture.mat SpinCurv Ek nks ftn58sparse K1 K2 K3 Nk1 Nk2 Nk3 Ta Tb Tc G1 G2 G3 -v7.3

%% --- nonzero Conductivity tensor Calculation --- %%%
Volume    = abs(dot(T(1,:),cross(T(2,:),T(3,:))));                      % in unit of (Angstrom)^3
mod_dV    = ((G1(2)-G1(1))*(G2(2)-G2(1))*(G3(2)-G3(1))*(2*pi)^3)*(1/Volume)*(1/nks);  % integral measure of the 3D BZ
mu        = linspace(Emu(1),Emu(2),nE);                                 % the energy window
sigma     = zeros(length(mu),27);
tic
fprintf('Start scanning chemical potentials ... \n')
for imu=1:length(mu)
    FD   = (exp((Ek-(Ef+mu(imu)))/(kB*Temp)) + 1).^-1;            % Fermi-Dirac distribution   

    % sum over occupied orbitals and integrate over k-space
    sigma(imu,:) = mod_dV*cquanta*(1/(2*pi)^2)*sum(FD.*SpinCurv,[1 2]);
    fprintf('%d\n',imu)
end
toc
save Conductivity.mat mu sigma Temp Nk1 Nk2 Nk3 -v7.3
%%
% poolobj = gcp('nocreate');
% delete(poolobj);

%% -------------------------- %%
%% table of tensor components %%
% 1: xxx
% 2: xxy
% 3: xxz
% 4: yxx (jyx,vx; jyx = {sx,vy}/2)
% 5: yxy
% 6: yxz
% 7: zxx
% 8: zxy
% 9: zxz
% 10: xyx
% 11: xyy
% 12: xyz
% 13: yyx
% 14: yyy
% 15: yyz
% 16: zyx
% 17: zyy
% 18: zyz
% 19: xzx
% 20: xzy
% 21: xzz
% 22: yzx
% 23: yzy
% 24: yzz
% 25: zzx
% 26: zzy
% 27: zzz
%% -------------------------- %%
