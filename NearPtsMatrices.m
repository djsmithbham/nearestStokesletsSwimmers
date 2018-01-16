function [xmat,Xmat,matNear,idx]=NearPtsMatrices(x,X,rho)

%idx = rangesearch(X,Y,r)

% x is a 3N matrix of force points
% X is a 3Q matrix of quadrature points
% rho is a threshold distance
% xmat and Xmat are 3N by 3Q sparse matrices
% which are populated only where x and X are within
% a distance rho

% true zeros will be used to omit values so set zero coordinates to very small numerical value instead

x(x==0)=2e-15;
X(X==0)=2e-15;

x=x(:);
X=X(:);

N=length(x)/3;
Q=length(X)/3;

xmat=sparse(3*N,3*Q);
Xmat=sparse(3*N,3*Q);

x1=x(1:N);
x2=x(N+1:2*N);
x3=x(2*N+1:3*N);

X1=X(1:Q);
X2=X(Q+1:2*Q);
X3=X(2*Q+1:3*Q);

if rho>0
    idx=rangesearch([X1 X2 X3],[x1 x2 x3],rho);
else
    idx=cell(N); % empty cell array
end

cidx=cell2mat(idx');
ridx=zeros(1,length(cidx));
s=ones(1,length(cidx));
counter=1;
for m=1:N
    ridx(counter:counter+length(idx{m})-1)=m*ones(1,length(idx{m}));
    counter=counter+length(idx{m});
end
matNear=sparse(ridx,cidx,s,N,Q);

Xmatr{1}=repmat(X1',N,1).*matNear;
Xmatr{2}=repmat(X2',N,1).*matNear;
Xmatr{3}=repmat(X3',N,1).*matNear;
Xmat=blkdiag(Xmatr{:});

xmatr{1}=repmat(x1,1,Q).*matNear;
xmatr{2}=repmat(x2,1,Q).*matNear;
xmatr{3}=repmat(x3,1,Q).*matNear;
xmat=blkdiag(xmatr{:});

xd=xmat-Xmat;
xd1=xd(1:N,1:Q);
xd2=xd(N+1:2*N,Q+1:2*Q);
xd3=xd(2*N+1:3*N,2*Q+1:3*Q);
