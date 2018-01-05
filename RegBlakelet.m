function S=RegBlakelet(x,X,eps)
% regularised Blakelet, i.e. Ainley et al. 2008 solution (see Smith 2009 for sign error correction)
% x is a row vec of field points: 3*collocation points, M
% X is a row vec of source points: 3*quadrature points, N
% the matrix output is expected to be multipled by a N x M matrix of quadrature weights 

% primary regularised stokeslet

x=x(:);
X=X(:);

M=length(x);
N=length(X);

r=x*ones(1,N)-ones(M,1)*X';

X1=r(      1:M/3,       1:N/3); %x differences
X2=r(M/3+1:2*M/3, N/3+1:2*N/3); %y differences
X3=r(  2*M/3+1:M,   2*N/3+1:N); %z differences

rsq=X1.^2+X2.^2+X3.^2;
ireps3=1./(sqrt((rsq+eps^2)).^3);
iso=kron(eye(3),(rsq+2.0*eps^2).*ireps3);
Ireps3=kron(ones(3,3),ireps3);
dyad=[X1.*X1 X1.*X2 X1.*X3; X2.*X1 X2.*X2 X2.*X3; X3.*X1 X3.*X2 X3.*X3].*(Ireps3);

S=iso+dyad;

clear rsq ireps3 Ireps3 iso dyad

% regularised image system of Ainley

Y=X;
Y(2*N/3+1:N)=-Y(2*N/3+1:N); % reflect sources in x3=0 plane

R=x*ones(1,N)-ones(M,1)*Y'; % R=x-Y

X3=R(  2*M/3+1:M,   2*N/3+1:N); % recalculate X3

Rsq=X1.^2+X2.^2+X3.^2;
iReps=1./sqrt(Rsq+eps^2);
iReps3=iReps.^3;
iso=kron(eye(3),(Rsq+2.0*eps^2).*iReps3);
IReps3=kron(ones(3,3),iReps3);
dyad=[X1.*X1 X1.*X2 X1.*X3; X2.*X1 X2.*X2 X2.*X3; X3.*X1 X3.*X2 X3.*X3].*(IReps3);

% regularised stokeslet image

SI=-(iso+dyad);
S=S+SI; % running total (to save memory)
clear SI;

% Higher order terms...

%  Delta tensor
Delta=diag([1 1 -1]);

%  height of source points
h=repmat(X(2*N/3+1:N)',[M 3]);

% blob term

iReps5=iReps.^5;
phi=3*eps^2*iReps5;
BT=-2*h.^2.*kron(Delta,phi);
S=S+BT;
clear BT;

% potential source dipole


IReps5=kron(ones(3,3),iReps5);
PD=2*h.^2.*(kron(Delta,iReps3)-3*IReps5.*[X1.*X1 X1.*X2 -X1.*X3; X2.*X1 X2.*X2 -X2.*X3; X3.*X1 X3.*X2 -X3.*X3]);
S=S+PD;
clear PD;

% regularised stokes dipole

SD=2*h.*( [zeros(2*M/3,N);[X1.*(Rsq+4*eps^2).*iReps5 X2.*(Rsq+4*eps^2).*iReps5 -X3.*(Rsq+4*eps^2).*iReps5]] - kron(Delta,X3.*iReps3) +    [zeros(M,2*N/3) [X1.*iReps3; X2.*iReps3; X3.*iReps3]] + 3*IReps5.*kron(ones(3),X3).*[X1.*X1 X1.*X2 -X1.*X3; X2.*X1 X2.*X2 -X2.*X3; X3.*X1 X3.*X2 -X3.*X3]);
S=S+SD;
clear SD;

% rotlet difference term

RD=-(6.0*h*eps^2.*IReps5).*([zeros(2*M/3,N);[X1 X2 X3]]-kron(eye(3),X3));
S=S+RD;
clear RD;

% final step: 1/8pi scaling
S=S*(1.0/(8.0*pi));

