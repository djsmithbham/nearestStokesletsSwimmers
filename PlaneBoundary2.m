function [x,X]=PlaneBoundary2(model)

% creates two opposed plane boundaries
% dimensions model.Lx by model.Ly
% number of points model.nx by model.ny (coarse)
%                  model.Nx by model.Ny (fine)
% model.O - origin of lower plane
% model.h - distance to upper plane

x=CreateSquareSurface(model.nx,model.ny);
X=CreateSquareSurface(model.Nx,model.Ny);

[x1,x2,x3]=ExtractComponents(x);
x1=x1*model.Lx;
x2=x2*model.Ly;
x=[x1;x2;x3];
x=TranslatePoints(x,model.O);
x=MergeVectorGrids(x,TranslatePoints(x,[0,0,model.h]));

[X1,X2,X3]=ExtractComponents(X);
X1=X1*model.Lx;
X2=X2*model.Ly;
X=[X1;X2;X3];
X=TranslatePoints(X,model.O);
X=MergeVectorGrids(X,TranslatePoints(X,[0,0,model.h]));