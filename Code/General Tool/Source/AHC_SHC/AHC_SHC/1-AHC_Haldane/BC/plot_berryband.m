clear all

load Berry_kpath.mat

EF=0;

berry1=Berry(:,:,1);
berry2=Berry(:,:,2);
berry3=Berry(:,:,3);
berry4=Berry(:,:,4);
berry5=Berry(:,:,5);
berry6=Berry(:,:,6);
berry7=Berry(:,:,7);
berry8=Berry(:,:,8);
berry9=Berry(:,:,9);

nk=length(Ek(:,1));
nb=length(Ek(1,:));

A=zeros(nk*nb,2);

for ikk=1:nk
    A((ikk-1)*nb+1:ikk*nb,1)=ikk;
    A((ikk-1)*nb+1:ikk*nb,2)=Ek(ikk,:)';
end

%return
%%
% Plot  band%
figure;
hold on
nknb=size(A(:,1));
scatter(A(:,1),A(:,2)-EF,50,reshape(berry2',nknb(1),1),'filled');
plot(Ek-EF,'k','LineWidth',2);
colorbar('vertical')
% title('berry2');
box on
axis square
colormap(colormap_rb_white_1)
caxis([-10 10]);
ax = gca;
ax.FontSize   = 14;
list    = [1 50 100 150 200];
ax.XTick      = list(:);
ax.XTickLabel = {'\Gamma' 'K' 'X' 'K^\prime' '\Gamma'};
ylim([-4 4]);
% Plot  band%
%{
figure;
hold on
nknb=size(A(:,1));
scatter(A(:,1),A(:,2)-EF,30,reshape(berry4',nknb(1),1),'filled');
plot(Ek-EF,'k');
colorbar('vertical')
title('berry4');
box on
axis square
colormap(colormap_rb_white_1)
%}
