function [x1,x2,x3]=GenerateSphereSurfaceForVisualisation(nth,nphi,a)

th=linspace(0,pi,nth);
phi=linspace(0,2*pi,nphi);

[th phi]=meshgrid(th,phi);

x1=a*cos(phi).*sin(th);
x2=a*sin(phi).*sin(th);
x3=a*cos(th);

