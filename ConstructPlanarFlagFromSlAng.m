function [x,v]=ConstructPlanarFlagFromSlAng(s,t,ThFun,dThFun,p)

% constructs discretisation of a planar-beating flagellum from sliding
% angle model
%
% ThFun(s,t,p) - function with model of sliding angle
% dThFun(s,t,p) - time derivative of sliding angle
% p - parameters structure for ThFun
% s - vector of N arclength values
% t - scalar time value
%
% output: x, contains 3N points, all x1 coords first, then x2 coords,...
%         v - time derivative of x

N=length(s);
x=zeros(3*N,1);

dx1=zeros(N-1,1);
dx2=zeros(N-1,1);
dv1=zeros(N-1,1);
dv2=zeros(N-1,1);

% split into integrals between each arclength
for j=1:N-1
    dx1(j,1)=integral(@(w)  cos(ThFun(w,t,p)),s(j),s(j+1));
    dx2(j,1)=integral(@(w)  sin(ThFun(w,t,p)),s(j),s(j+1));
    dv1(j,1)=integral(@(w) -sin(ThFun(w,t,p)).*dThFun(w,t,p),s(j),s(j+1));
    dv2(j,1)=integral(@(w)  cos(ThFun(w,t,p)).*dThFun(w,t,p),s(j),s(j+1));
end

% sum to get integrals 0 to s
x1=[zeros(1,N-1); tril(ones(N-1))]*dx1;
x2=[zeros(1,N-1); tril(ones(N-1))]*dx2;
v1=[zeros(1,N-1); tril(ones(N-1))]*dv1;
v2=[zeros(1,N-1); tril(ones(N-1))]*dv2;

% planar beat
x3=zeros(N,1);
v3=zeros(N,1);

x=[x1;x2;x3];
v=[v1;v2;v3];


