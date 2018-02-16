function [x,v,X]=SpermModelGI_headOnly(t,model)

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
x1=model.a1*xh1;
x2=model.a2*xh2;
x3=model.a3*xh3;
x =[x1;x2;x3];

v =[0;0;0];

%-------------------------------------------------------------------------
% fine grid - position only
%
% generate head position
Xh=GenerateSphereDiscr(model.Nh,1);
[Xh1,Xh2,Xh3]=ExtractComponents(Xh);
X1=model.a1*Xh1;
X2=model.a2*Xh2;
X3=model.a3*Xh3;

X =[X1;X2;X3];
