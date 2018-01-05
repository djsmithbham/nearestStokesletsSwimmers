function [x,z,dx,dz]=CiliumFromFourier(s,t,params)

% generates a cilium shape from Fourier representation

% params.M is the number of coefficients in the polynomial for cilium shape
% (for cubic M=4)

sigma=1;
ax0=params.ax0;
az0=params.az0;
ax=params.ax;
az=params.az;
bx=params.bx;
bz=params.bz;
M=params.M;
AAx=zeros(M,1);
AAz=zeros(M,1);
for m=1:M
    AAx(m)=FourierPolyEval(ax0{m},ax{m},bx{m},sigma,t);
    AAz(m)=FourierPolyEval(az0{m},az{m},bz{m},sigma,t);
end
x=polyval(AAx,s);
z=polyval(AAz,s);
dx=dpolyval(AAx,s);
dz=dpolyval(AAz,s);