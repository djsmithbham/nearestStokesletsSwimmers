function [x,v,X]=ChlamyModel1_headonly(t,model) %#ok<INUSL>
%
% generates discretisation of a model chlamy scalene ellipsoid head only
%
% Input:
%
% t - time
% model.ns - number of points along each flagellum
% model.nh - head discretisation parameter. Total points are 6 nh^2
% model.a1 etc - head semi-axes
% model.ang - angle for flagella insertion (0 = both emerging from x2 axis)
% model.F  - F{1} is x-interpolant, F{2} is y-interpolant
%
% Output:
%
% x - Coarse grid - head and flagella positions
% v - Coarse grid - head and flagella velocities
% X - Fine grid - head and flagellar positions only

%-------------------------------------------------------------------
% coarse grid - position and velocity

% generate head position
xh=GenerateSphereDiscr(model.nh,1);
[xh1,xh2,xh3]=ExtractComponents(xh);
xh1=model.a1*xh1;
xh2=model.a2*xh2;
xh3=model.a3*xh3;
x=[xh1;xh2;xh3];

% head is stationary in body frame
vh1=0*xh1;
vh2=0*xh2;
vh3=0*xh3;
v=[vh1;vh2;vh3];

%-------------------------------------------------------------------------
% fine grid - position only
%
% generate head position
Xh=GenerateSphereDiscr(model.Nh,1);
[Xh1,Xh2,Xh3]=ExtractComponents(Xh);
Xh1=model.a1*Xh1;
Xh2=model.a2*Xh2;
Xh3=model.a3*Xh3;
X=[Xh1;Xh2;Xh3];
