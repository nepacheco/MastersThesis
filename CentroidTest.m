clc;clear;close all
%% Running script to find neutral bending plan of off axis cutting segment
ro = 1;
ri = 0.9;
g = 0.85*2*ro;

phio = 2*acos((g-ro)/ro);
phii = 2*acos((g-ro)/ri);

ybaro = 4*ro*sin(phio/2)^3/(3*(phio - sin(phio)));
ybari = 4*ri*sin(phii/2)^3/(3*(phii - sin(phii)));
Ao = ro^2*(phio - sin(phio))/2;
Ai = ri^2*(phii - sin(phii))/2;
A = Ao - Ai;
ybar = (Ao*ybaro - Ai*ybari)/A;

figure()
xo = linspace(-sin(phio/2)*ro,sin(phio/2)*ro,100);
xi = linspace(-sin(phii/2)*ri,sin(phii/2)*ri,100);
xinter1 = linspace(xo(1),xi(1),2);
xinter2 = linspace(xi(end),xo(end),2);
yo = sqrt(ro^2 - xo.^2);
yi = sqrt(ri^2 - xi.^2);
yinter1 = linspace(yo(1),yi(1),2);
yinter2 = linspace(yi(end),yo(end),2);
xbar = xo;
hold on
title("Comparing On-Axis and Off-Axis Neutral Bending Planes")
plot(xo,yo,'b')
plot(xi,yi,'b')
plot(xinter1,yinter1,'b',xinter2,yinter2,'b');
plot(xbar,ybar*ones(100,1),'b');
axis equal
grid on;
%% Running a script to find neutral bending plane of on axis cutting segment
phi = phii;
ro = 1;
ri = 0.9;
ybaro = 2*ro*sin(phi/2)/(3*phi/2);
ybari = 2*ri*sin(phi/2)/(3*phi/2);
Ao = phi/2*ro^2;
Ai = phi/2*ri^2;
A = Ao - Ai;
ybar = (Ao*ybaro - Ai*ybari)/A;

xo = linspace(-sin(phi/2)*ro,sin(phi/2)*ro,100);
xi = linspace(-sin(phi/2)*ri,sin(phi/2)*ri,100);
xinter1 = linspace(xo(1),xi(1),2);
xinter2 = linspace(xi(end),xo(end),2);
yo = sqrt(ro^2 - xo.^2);
yi = sqrt(ri^2 - xi.^2);
yinter1 = linspace(yo(1),yi(1),2);
yinter2 = linspace(yi(end),yo(end),2);
xbar = xo;
plot(xo,yo,'r')
plot(xi,yi,'r')
plot(xinter1,yinter1,'r',xinter2,yinter2,'r');
plot(xbar,ybar*ones(100,1),'r');
xlabel("mm")
ylabel("mm");
axis equal
grid on;