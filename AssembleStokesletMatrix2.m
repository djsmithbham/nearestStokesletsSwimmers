function [AS,NN]=AssembleStokesletMatrix2(x,X,ep,domain,rho,blockSize)

N=length(x)/3;
Q=length(X)/3;

[~,~,matNear,idx]=NearPtsMatrices(x,X,rho);  

% walltime(1)=wt(1);
% walltime(2)=wt(2);
% walltime(3)=wt(3);

% x1=x(1:N);
% x2=x(N+1:2*N);
% x3=x(2*N+1:3*N);
% 
% X1=X(1:Q);
% X2=X(Q+1:2*Q);
% X3=X(2*Q+1:3*Q);

% calculate index of near quadrature points (within distance rho - symbol psi)
% if rho>0
%     idx=rangesearch([X1 X2 X3],[x1 x2 x3],rho);
% else
%     idx=cell(N); % empty cell array
% end

% calculate nearest neighbour matrix (symbol nu)
NN=NearestNeighbourMatrix(X,x,blockSize);
    
% assemble near stokeslets matrix
ASnear=RegStokeslet2(x,X,idx,ep)*NN;

%
% blocking... haven't figured out how to do this yet in the present context
% calculate number of quadrature nodes that can be used for each block
% blockNodes=floor(blockSize*2^27/(9*N));
% ASnear=zeros(3*N,3*N);
% for iMin=1:blockNodes:Q
%     iMax=min(iMin+blockNodes-1,Q);
%     iRange=[iMin:iMax Q+iMin:Q+iMax 2*Q+iMin:2*Q+iMax];
%     switch domain
%            case 'i'
%                ASnear=ASnear+RegStokeslet2(x,X(iRange),idx,ep)*NN(iRange,:);
%            case 'h'
%                ASnear=ASnear+RegBlakelet2( x,X(iRange),idx,ep)*NN(iRange,:);
%     end
% end

% assemble far stokeslets matrix
matFar=sum(NN,1)-kron(ones(3),matNear)*NN; % 3N x 3N matrix with weights for far stokeslets  
AScoarse=RegStokeslet(x,x,ep).*matFar;
  
AS=ASnear+AScoarse;
