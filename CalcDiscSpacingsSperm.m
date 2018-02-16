function [hf,hq,f_hf,f_hq,h_hf,h_hq,f_min,q_min] = CalcDiscSpacingsSperm(zz,tt,swimmer,boundary)

% Head + flagella
x0=zz(1,1:3);
b1=zz(1,4:6);
b2=zz(1,7:9);
b3=cross(b1,b2);
B=[b1(:) b2(:) b3(:)];

[xi,~,Xi]=swimmer.fn(tt(1),swimmer.model);
xx0=ApplyRotationMatrix(B,xi); % x-x0

xs=TranslatePoints(xx0,x0);
Xs=ApplyRotationMatrix(B,Xi);
Xs=TranslatePoints(Xs,x0);

if isempty(boundary)
    xb=[];  
    Xb=[];
else
    [xb,Xb]=boundary.fn(boundary.model);
end

x=MergeVectorGrids(xs,xb);
X=MergeVectorGrids(Xs,Xb);

hf = CalcDiscr_h(x,0.2);
hq = CalcDiscr_h(X,0.2);

% Head only
x0=zz(1,1:3);
b1=zz(1,4:6);
b2=zz(1,7:9);
b3=cross(b1,b2);
B=[b1(:) b2(:) b3(:)];

[xi,~,Xi]=SpermModelGI_headOnly(tt(1),swimmer.model);
xx0=ApplyRotationMatrix(B,xi); % x-x0

xs=TranslatePoints(xx0,x0);
Xs=ApplyRotationMatrix(B,Xi);
Xs=TranslatePoints(Xs,x0);

if isempty(boundary)
    xb=[];
    Xb=[];
else
    [xb,Xb]=boundary.fn(boundary.model);
end

xHead=MergeVectorGrids(xs,xb);
XHead=MergeVectorGrids(Xs,Xb);

h_hf = CalcDiscr_h(xHead,0.2);
h_hq = CalcDiscr_h(XHead,0.2);

% Flagellum only
x0=zz(1,1:3);
b1=zz(1,4:6);
b2=zz(1,7:9);
b3=cross(b1,b2);
B=[b1(:) b2(:) b3(:)];

[xi,~,Xi]=SpermModelGI_flagOnly(tt(1),swimmer.model);
xx0=ApplyRotationMatrix(B,xi); % x-x0

xs=TranslatePoints(xx0,x0);
Xs=ApplyRotationMatrix(B,Xi);
Xs=TranslatePoints(Xs,x0);

if isempty(boundary)
    xb=[];
    Xb=[];
else
    [xb,Xb]=boundary.fn(boundary.model);
end

xFlag=MergeVectorGrids(xs,xb);
XFlag=MergeVectorGrids(Xs,Xb);

f_hf = CalcDiscr_h(xFlag,0.2);
f_hq = CalcDiscr_h(XFlag,0.2);

f_min = CalcDiscr_hMin(xHead,xFlag,0.2);
q_min = CalcDiscr_hMin(XHead,XFlag,0.2);

end