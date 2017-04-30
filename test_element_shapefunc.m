
clc
clear all
close all

ex=[850 1000 1000];
ey=[900 700 800];

alpha=1;
ep=1;
T_inf=20;
% 
%  n=1
%  [Ke, fe]=robin_heat_tria(ex,ey,ep,alpha,T_inf,n)
%  [Ke2, fe2]=robin_heat_tria2(ex,ey,ep,alpha,T_inf,n)

n=2
[Ke, fe]=robin_heat_tria(ex,ey,ep,alpha,T_inf,n)
[Ke2, fe2]=robin_heat_tria2(ex,ey,ep,alpha,T_inf,n)

 n=3
 [Ke, fe]=robin_heat_tria(ex,ey,ep,alpha,T_inf,n)
 [Ke2, fe2]=robin_heat_tria2(ex,ey,ep,alpha,T_inf,n)





