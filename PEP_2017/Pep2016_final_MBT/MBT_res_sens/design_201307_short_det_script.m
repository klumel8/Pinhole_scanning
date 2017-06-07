clc; clear all;

tic
L=40; 
R_target=6.0;
do_seq=0;
Ri=3.2;
d_factor=1.0;

[Sens,d,nph,S,rph,alpha_yz,axisph]=design_201307_short_det(L,R_target,do_seq,Ri,d_factor);
toc

figure(1);imagesc(squeeze(S(21,:,:)));axis image
% squeeze gooit overbodige dimenties weg (1x5x8 matrix --> 5x8 matrix)
% imagesc geeft kleur aan waarden van laag tot hoog
% axis image maakt rechthoekige pixels vierkant