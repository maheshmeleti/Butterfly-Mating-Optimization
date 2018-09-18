clc
clear all
close all
x=1:1000;
y=sin(2*pi*0.005*x.*0.75);
plot(x,y);
figure(2);
y=cos(2*pi*0.005*x);
plot(x,y);