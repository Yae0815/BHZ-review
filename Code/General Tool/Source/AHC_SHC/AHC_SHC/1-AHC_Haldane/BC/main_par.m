clc
clear

tic;
%%   load TBmodel ftn58sparse
load ftn58sparse_soc.mat
Nk=80;
EF=-3.7016;  %Fermi-level for Bi-H model
% nb=ftn58(1,1);
% ii=ftn58(2:end,2);
% jj=ftn58(2:end,3);
% tt=ftn58(2:end,4);
% dd=ftn58(2:end,5:7);
nb=ftn58sparse.norb;
ii=ftn58sparse.ij(:,1);
jj=ftn58sparse.ij(:,2);
tt=ftn58sparse.tt;
dd=ftn58sparse.dd;
%%  setup Kmesh line mode G-M-K-G
kxp=[0  0.5 1/3 0]*2*pi; 
kyp=[0  0.0 1/3 0]*2*pi; 
kzp=[0  0.0 0.0 0]*2*pi; 
kx=[];ky=[];kz=[];
for ikp=1:length(kxp)-1
    kx=cat(2,kx,linspace(kxp(ikp),kxp(ikp+1),Nk+1));kx(end)=[];
    ky=cat(2,ky,linspace(kyp(ikp),kyp(ikp+1),Nk+1));ky(end)=[];
    kz=cat(2,kz,linspace(kzp(ikp),kzp(ikp+1),Nk+1));kz(end)=[];
end
kpath=[kx(:),ky(:),kz(:)];
%  load kpath.mat

kpts=kpath;
nk=length(kpts);
Ek=zeros(nk,nb);

for ikk=1:nk
    fprintf('total k=%d now=%d\n',nk,ikk);
    kpoints=kpts(ikk,:);
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,nb,nb);
    H0      = full(Hsparse);
    H0      = (H0 + H0')/2;
    [V,D]   = eig(H0);
    Etmp    = diag(D);
    Ek(ikk,:)=Etmp;
end


toc;
Computation_time=num2str(toc); 

save bandEk.mat Ek

% return
%%  plot
figure,hold on,axis square,box on
set(gca,'linewidth',2,'FontSize',16)
ylabel('Energy (eV)','interpreter','latex')
for iil=1:ikp-1
    plot([Nk*iil Nk*iil],[-2 2],'k--')
end
axis([0 Nk*ikp -2 2])
set(gca,'XTick',[1 80 160 240],'XTickLabel',{'\Gamma','M','K','\Gamma'})
plot(Ek-EF,'k','LineWidth',2)



