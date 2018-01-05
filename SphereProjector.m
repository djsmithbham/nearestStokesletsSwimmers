function [xp,yp,zp]=SphereProjector(x,y,z,a)

r=sqrt(x.^2+y.^2+z.^2);
th=acos(z./r);
phi=acos((x+1.0e-9)./sqrt((x+1.0e-9).^2+y.^2)).*(2*real(y>0)-1)+2*pi*real(y<0);
xp=a*sin(th).*cos(phi);
yp=a*sin(th).*sin(phi);
zp=a*cos(th);

