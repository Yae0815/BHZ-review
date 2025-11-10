clear
%function [ftn58soc_full]=sparse2ftn58full(ftn58sparse)
load ftn58sparse.mat
% ftn58sparse=Sftn58sparse;
% clear Sftn58sparse

ftn58soc_full=zeros(length(ftn58sparse.ij)+1,7);
ftn58soc_full(1,1)=ftn58sparse.norb;
ftn58soc_full(1,2)=length(ftn58sparse.ij);
ftn58soc_full(2:end,1)=[1:length(ftn58sparse.ij)]';
ftn58soc_full(2:end,2:3)=real(ftn58sparse.ij);
ftn58soc_full(2:end,4)=ftn58sparse.tt;
ftn58soc_full(2:end,5:7)=real(ftn58sparse.dd);
% ftn58soc_full(2:end,6)=0;
% ftn58soc_full(2:end,6)=0;
%return

save ftn58soc_full.mat -v7.3 ftn58soc_full
