%%%          Band Structure Plot for ftn58sparse         %%%
%%% The spin polariztion is implemented in this function %%%
%%% 3/6/2016 Hans                                        %%%
%%% ---------------------------------------------------- %%%
clear all

E_range = [-15 15];
Ef      = 0;%-16;
isSP    = 0;
is2D    = 0;

%% Initial info. %%
load ftn58sparse.mat
norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;
BR   = ftn58sparse.BR;
Sz   = [1 0;0 -1];

%% Kpoints %%%
nk = 100;

n = 0.5;
p1 = [linspace(0.0,-1/3,nk+1)' linspace(0.0,-2/3,nk+1)' linspace(0.0,0.0,nk+1)'];
p2 = [linspace(-1/3,0.0,nk+1)' linspace(-2/3,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
p3 = [linspace(0.0,1/3,nk+1)' linspace(0.0,2/3,nk+1)' linspace(0.0,0.0,nk+1)'];
p4 = [linspace(1/3,0.0,nk+1)' linspace(2/3,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
%p7 = [linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)' linspace(-0.5,0.0,nk+1)'];
%p8 = [linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.5,nk+1)'];
list    = [1 101 201 301 400];
%kpoints = [p1(1:nk,:);p2(1:nk,:);p3(1:nk,:);p4(1:nk,:);p5(1:nk,:);p6(1:nk,:);p7(1:nk,:);p8(1:nk,:)]*2*pi;
kpoints = [p1(1:nk,:);p2(1:nk,:);p3(1:nk,:);p4(1:nk,:)]*2*pi;
%kpoints = [p1(1:nk,:);p2(1:nk,:)]*2*pi;

if is2D==1
    N1 = 251; N2 = 251;
    p1 = linspace(-0.0,0.2,N1); 
    p2 = linspace(-0.0,0.2,N2);  
    [K1,K2] = meshgrid(p1,p2);

    kk = 1;
    for n1=1:N1
        for n2=1:N2
            kpoints(kk,:) = [p1(n1) p2(n2) 0]*2*pi;
            kk = kk + 1;
        end
    end

    nks     = length(kpoints);
end
    

%% Eigenvalue and Eigenvector calculations %%%
tic
eigvec = cell(size(kpoints,1),1);
nks    = length(kpoints);
Ek     = zeros(nks,norb);
SPz    = zeros(nks,norb);
for ik=1:nks
    time_init=tic;
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints(ik,:)').*tt,norb,norb);
    HH      = full(Hsparse);
    HH      = (HH+HH')/2 ;
    [vec, Etemp] = eig(HH);
     
    Ham{ik,1}    = HH;
    eigvec{ik,1} = vec;

%     A(ik*norb-norb+1:ik*norb,1) = ik;
%     A(ik*norb-norb+1:ik*norb,2) = diag(Etemp);
    
    Ek(ik,:)  = diag(Etemp);
    SSPz      = vec'*kron(Sz,eye(norb/2))*vec;
    SSPz      = (SSPz+SSPz')/2;
    SPz(ik,:) = diag(SSPz); 
    
    if mod(ik,1e3)==0
        fprintf('%3i/%i: %.3fs [%.2fm]\n',ik,nks,toc,toc(time_init)/60);
    end
    
end
toc

if is2D==1
    Ekk = reshape(Ek,[N1,N2,4]);

    %% 2D-Band
    figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
           'papertype','usletter','numbertitle','off','name','kx_kz',...
           'PaperPositionMode','manual','paperorientation','portrait',...
           'color','w');

    for iband = [2 3]
    %     surface(K1,K2,Ekk(:,:,iband)-Ef,'LineStyle','-','LineWidth',0.01,...
    %             'EdgeColor','none');
        surf(K1,K2,Ekk(:,:,iband)-Ef,'FaceAlpha','interp',...
            'FaceColor','blue',...    
            'AlphaDataMapping','scaled',...
            'AlphaData',gradient(abs(Ekk(:,:,iband)-Ef))*1e2,...        
            'EdgeColor','none');
        hold on
    end

    lighting gouraud;
    box on;
    ax = gca;
    ax.XGrid= 'on';
    ax.YGrid= 'on';
    ax.ZGrid= 'on';
    ax.DataAspectRatio = [1 1 15];
    axis([K1(1) K1(end) K2(1) K2(end) -1 1]);

    set(gca,'FontSize',16,'TickLabelInterpreter','latex','FontWeight','bold');
    xlabel('\bf{k$_{x}$($2\pi/a$)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    ylabel('\bf{k$_{y}$($2\pi/a$)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    zlabel('\bf{Energy (eV)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    return
end

% [~, Emin] = min(abs(Ek(1,:)-Ef));
% disp(Emin);

%% Plotting %%%
figure
color1 = [0.8500 0.3250 0.0980];
if isSP==1
    ncolor   = 1e4;
    spmap_p  = [linspace(1,1,ncolor+1)' linspace(0,1,ncolor+1)' linspace(0,1,ncolor+1)'];
    spmap_n  = [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)'];
    spmap    = [spmap_p(1:ncolor,:);spmap_n];

    for ii=1:norb
        for ikk = 1:nks
            if SPz(ikk,ii) > 0
                scatter(ikk,Ek(ikk,ii)-Ef,SPz(ikk,ii),'r.','LineWidth',1);
            elseif SPz(ikk,ii) < 0
                scatter(ikk,Ek(ikk,ii)-Ef,-SPz(ikk,ii),'b.','LineWidth',1);
            end
        hold on
        end
    end
    colorbar('TickLabelInterpreter','LaTex','FontSize',18);
else
    %for ii=1:norb
        %plot(1:nks,Ek(:,1)-Ef,'r','LineWidth',1.5);
        %hold on
        %plot(1:nks,Ek(:,3)-Ef,'b','LineWidth',1.5);
        hold on
        %plot(1:nks,Ek-Ef,'k','LineWidth',1.5);
        plot(1:nks,Ek-Ef,'color','k','LineWidth',1.5);
    %end
end

for il = 2:size(list,2)-1
line('XData', [list(il) list(il)], 'YData', [E_range(1) E_range(2)], 'LineStyle', '-', ...
    'LineWidth', 0.1, 'Color','k');
end
line('XData', [0 nks], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');

axis([1 nks E_range(1) E_range(2)]);
ax = gca;
ax.TickDir    = 'in';
ax.FontSize   = 14;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick      = list(:);
ax.YTick      = E_range(1):10:E_range(2);
ax.XTickLabel = {'\bf{$\Gamma$}' '-K' '\bf{$\Gamma$}' 'K' '\bf{$\Gamma$}'};
% ax.XTickLabel = {'\bf{X}' '\bf{$\Gamma$}' '\bf{X}'};
% ax.XTickLabel = {};
ax.LineWidth  = 1.2;
ax.TickLabelInterpreter='latex';
ylabel('E(eV)','FontSize',18,'interpreter','LaTex');