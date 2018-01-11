function x=CreateSquareSurface(nx,ny)

% creates a unit square surface with nx, ny points per edge

sx=linspace(-0.5,0.5,nx);
sy=linspace(-0.5,0.5,ny);

[x1,x2]=ndgrid(sx,sy);
x3=0*x1;

x=[x1(:);x2(:);x3(:)];