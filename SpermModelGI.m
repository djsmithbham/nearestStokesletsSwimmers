function [x,v,X]=SpermModelGI(t,model)

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

% generate head position
xh=GenerateSphereDiscr(model.nh,1);
[xh1,xh2,xh3]=ExtractComponents(xh);
xh1=model.a1*xh1;
xh2=model.a2*xh2;
xh3=model.a3*xh3;

% head is stationary in body frame
vh1=0*xh1;
vh2=0*xh2;
vh3=0*xh3;

% flagellum
s=linspace(0,1,model.ns);
% finite differences for velocity... can't be too precise as x,y are
% outputs of fsolve
dt=0.001;
tt=[t-dt/2  t  t+dt/2];

[sg,tg]=ndgrid(s,tt); % data goes in model.F
[xg,yg]=CalcxyFromPlanarInterp(sg,tg,model.F);

xt1=xg(:,2);xt2=yg(:,2);xt3=0*xg(:,2);

vt1=(xg(:,3)-xg(:,1))/dt;
vt2=(yg(:,3)-yg(:,1))/dt;
vt3=0*xg(:,2);

xt1=model.a1+xt1;

x1=[xh1;xt1];
x2=[xh2;xt2];
x3=[xh3;xt3];
x =[x1;x2;x3];

v1=[vh1;vt1];
v2=[vh2;vt2];
v3=[vh3;vt3];
v =[v1;v2;v3];

%-------------------------------------------------------------------------
% fine grid - position only
%
% generate head position
Xh=GenerateSphereDiscr(model.Nh,1);
[Xh1,Xh2,Xh3]=ExtractComponents(Xh);
Xh1=model.a1*Xh1;
Xh2=model.a2*Xh2;
Xh3=model.a3*Xh3;

% flagellum
S=linspace(0,1,model.Ns);

[Sg,tg]=ndgrid(S,tt(:,2));
[Xt1,Xt2]=CalcxyFromPlanarInterp(Sg,tg,model.F);
Xt3=0*Xt1;

Xt1=model.a1+Xt1;

X1=[Xh1;Xt1];
X2=[Xh2;Xt2];
X3=[Xh3;Xt3];
X =[X1;X2;X3];
