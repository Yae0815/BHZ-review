load SpinCurvatureTexture.mat

%% Inputs
orb = 3;      % to plot the Curvaturetexture of the orb^th orbital.
ii  = 22;      % the component of the Spin Curvature tensor
W   = 1;     % the intensity of the Curvaturetexture
cax = [-4 4]; % the range of the colormap
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

%% --- Plot Figure --- %%%
ncolor = 1e4;
%red <=> blue
map      = [[linspace(1,1,ncolor)' linspace(0,1,ncolor)' linspace(0,1,ncolor)'];...
            [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)']];

figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');

SpinCurv_ii = SpinCurv(:,orb,ii);
SpinCurv_ii = reshape(SpinCurv_ii,[Nk1 Nk2]);

h = pcolor(K1*(2*pi*(1/Ta)),K2*(2*pi*(1/Tb)),SpinCurv_ii*W);,shading interp;
%rotate(h,[0 0 1],45);
colormap(map);
caxis(cax);
colorbar('Ticks',cax,'TickLabels',{num2str(cax(1)),num2str(cax(2))});
box on

title(['Spin Curvaturetexture of ',num2str(orb),' th orbital at kz=0'])

axis([G1(1)/Ta G1(2)/Ta G2(1)/Tb G2(2)/Tb]*2*pi);
xlabel('$\bf{k_{x} (\AA)^{-1}}$','FontSize',24,'interpreter','LaTex');
ylabel('$\bf{k_{y} (\AA)^{-1}}$','FontSize',24,'interpreter','LaTex');
ax = gca;
ax.FontSize   = 20;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
%ax.XTickLabel = { };
ax.LineWidth  = 2;
ax.TickLabelInterpreter='latex';