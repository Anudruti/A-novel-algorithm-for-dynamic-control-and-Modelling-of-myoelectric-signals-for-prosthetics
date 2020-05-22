clc
clear all;
close all;
num=[0 1];
den=[3.51 1.3 0];
syg=tf(num,den);
bode(syg);
[Gm,Pm,Wgc,Wpc]=margin(syg)
grid on;