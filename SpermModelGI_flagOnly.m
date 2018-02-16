function [x,v,X]=SpermModelGI_flagOnly(t,model)

% generates discretisation of a model sperm, griddedInterpolant flagellum,
% scalene ellipsoid head
%
% t - time
% model.ns - number of points along flagellum
% model.nh - head discretisation parameter. Total points are 6 nh^2
% model.a1 etc - head semi-axes
% model.F  - F{1} is x-interpolant, F{2} is y-interpolant

%-------------------------------------------------------------------
% coarse grid - position and velocity

% flagellum
s=linspace(0,1,model.ns);
% finite differences for velocity... can't be too precise as x,y are
% outputs of fsolve
dt=0.001;
tt=[t-dt/2  t  t+dt/2];

[sg,tg]=ndgrid(s,tt); % data goes in model.F
[xg,yg]=CalcxyFromPlanarInterp(sg,tg,model.F);

x1=xg(:,2);
x2=yg(:,2);
x3=0*xg(:,2);

x =[x1;x2;x3];

v =[0;0;0];

%-------------------------------------------------------------------------
% fine grid - position only
%
% flagellum
S=linspace(0,1,model.Ns);

[Sg,tg]=ndgrid(S,tt(:,2));
[X1,X2]=CalcxyFromPlanarInterp(Sg,tg,model.F);
X3=0*X1;

X1=model.a1+X1;

X =[X1;X2;X3];
