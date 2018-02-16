function [x,v,X]=ChlamyModel1_flagonly(t,model)

% generates discretisation of a model chlamy, griddedInterpolant flagella,
% flagella are synchronised, scalene ellipsoid head
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
xh=[xh1;xh2;xh3];

% head is stationary in body frame
vh1=0*xh1;
vh2=0*xh2;
vh3=0*xh3;
vh=[vh1;vh2;vh3];

% flagellum
s=linspace(0,1,model.ns);
% finite differences for velocity... can't be too precise as x,y are
% outputs of fsolve
dt=0.001;
tt=[t-dt/2  t  t+dt/2];

[sg,tg]=ndgrid(s,tt); % data goes in model.F
[xg,yg]=CalcxyFromPlanarInterp(sg,tg,model.F);

% duplicate flagella, flip the trans flagellum
% xt1 xt2 xt3 - trans
% xc1 xc2 xc3 - cis
xt1=-xg(:,2); xt2=-yg(:,2); xt3=0*xg(:,2); xt=[xt1;xt2;xt3];
xc1= xg(:,2); xc2=-yg(:,2); xc3=0*xg(:,2); xc=[xc1;xc2;xc3];

% rotate the flagella
X0=[0;0;0];axx=3;
xt=RotatePoints(xt,X0, model.ang,axx);
xc=RotatePoints(xc,X0,-model.ang,axx);

% translate to surface of chlamy
xt=TranslatePoints(xt,[-model.a1*sin(model.ang),model.a2*cos(model.ang),0]);
xc=TranslatePoints(xc,[ model.a1*sin(model.ang),model.a2*cos(model.ang),0]);

% merge flagella positions
x=MergeVectorGrids(xt,xc);

% calculate velocities (reversing trans flagellum in v1 component)
vt1=(-xg(:,3)+xg(:,1))/dt; vt2=-(yg(:,3)-yg(:,1))/dt; vt3=0*xg(:,2); vt=[vt1;vt2;vt3];
vc1=( xg(:,3)-xg(:,1))/dt; vc2=-(yg(:,3)-yg(:,1))/dt; vc3=0*xg(:,2); vc=[vc1;vc2;vc3];

% rotate the velocities
vt=RotatePoints(vt,X0, model.ang,axx);
vc=RotatePoints(vc,X0,-model.ang,axx);

% merge flagella velocities
v=MergeVectorGrids(vt,vc);

%-------------------------------------------------------------------------
% fine grid - position only
%
% generate head position
Xh=GenerateSphereDiscr(model.Nh,1);
[Xh1,Xh2,Xh3]=ExtractComponents(Xh);
Xh1=model.a1*Xh1;
Xh2=model.a2*Xh2;
Xh3=model.a3*Xh3;
Xh=[Xh1;Xh2;Xh3];

% flagellum
S=linspace(0,1,model.Ns);

% duplicate flagella, flip the trans flagellum
% Xt1 Xt2 Xt3 - trans
% Xc1 Xc2 Xc3 - cis
[Sg,tg]=ndgrid(S,tt(:,2));
[Xg1,Xg2]=CalcxyFromPlanarInterp(Sg,tg,model.F);
Xt1=-Xg1; Xt2=-Xg2; Xt3=0*Xg1; Xt=[Xt1;Xt2;Xt3];
Xc1= Xg1; Xc2=-Xg2; Xc3=0*Xg1; Xc=[Xc1;Xc2;Xc3];

% rotate the flagella
X0=[0;0;0];axx=3;
Xt=RotatePoints(Xt,X0, model.ang,axx);
Xc=RotatePoints(Xc,X0,-model.ang,axx);

% translate to surface of chlamy
Xt=TranslatePoints(Xt,[-model.a1*sin(model.ang),model.a2*cos(model.ang),0]);
Xc=TranslatePoints(Xc,[ model.a1*sin(model.ang),model.a2*cos(model.ang),0]);

% merge flagella positions
X=MergeVectorGrids(Xt,Xc);

