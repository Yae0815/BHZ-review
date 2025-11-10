load Conductivity.mat

FIBZ = 1; % the constant to recover the total conductivity from irreducible BZ
ii   = 22; % the ii^th tensor component

color1 = [0.8500 0.3250 0.0980];
%% -------------------------- %%
%% table of tensor components %%
% 1: xxx  
% 2: xxy
% 3: xxz
% 4: yxx  (jyx,vx; jyx = {sx,vy}/2)
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
%% --- Conductivity -- %%%
% figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
%        'papertype','usletter','numbertitle','off',...
%        'PaperPositionMode','manual','paperorientation','portrait',...
%        'color','w');
% clf
figure
hold on
%plot(mu,FIBZ*sigma(:,ii),'k','LineStyle','-','Linewidth',1.0);
plot(mu,FIBZ*sigma(:,ii),'k','LineStyle','-','Linewidth',2.0);

%line('XData',[mu(1) mu(end)], 'YData', [0 0], 'LineStyle', '--', ...
%    'LineWidth', 0.5, 'Color','k');
%title(['Conductivity at T=',num2str(Temp),' with mesh ',num2str([Nk1 Nk2])])

% ls = legend('\bf{$\sigma_{xxy}$}');
% set(ls,'interpreter','LaTex','FontSize',14);
set(gca,'Ticklength',[0.05 0]);
set(gca,'Fontsize',14);
xlabel('$\mu$ (eV)','interpreter','LaTex');
ylabel('\bf{$\sigma^z_{xy}$} \bf{$(e^2/h)$}','interpreter','LaTex');
ax = gca;
ax.TickLabelInterpreter='latex';
ax.FontSize = 14;
ylim([-1.5 1.5]);
yline(1,'r--','LineWidth',1.6)
view(90,-90);
box on;