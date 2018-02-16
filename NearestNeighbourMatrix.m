function NClosest=NearestNeighbourMatrix(xQuad,xTrac,varargin)

	% Given 3Q x 1 vector of quadrature points xQuad
	% and a 3N x 1 vector of traction points xTrac
	% this function creates a 3Q x 3N sparse matrix
	% with a single '1' in each column to indicate the closest
	% traction point to each quadrature point.
	%
	% Vectors should be supplied with all x1 coordinates listed first
	% then all x2 coordinates, then all x3 coordinates. 
    
    % if varargin is nonempty, then it should contain blockSize
    
    Q=length(xQuad)/3;
	N=length(xTrac)/3;
    
    if ~isempty(varargin)
        blockSize=varargin{1};
        % calculate number of quadrature nodes that can be done in each 
        %   block
        blockNodes=floor(blockSize*2^27/(9*N));
    else
        blockNodes=Q;
    end
    
    xQ1=xQuad(1:Q);
    xQ2=xQuad(Q+1:2*Q);
    xQ3=xQuad(2*Q+1:3*Q);
    
    xT1=xTrac(1:N);
    xT2=xTrac(N+1:2*N);
    xT3=xTrac(2*N+1:3*N);
        
    nMin=zeros(Q,1);
    
%  	X1=xQ1*ones(1,N)-ones(Q,1)*xT1';
%	X2=xQ2*ones(1,N)-ones(Q,1)*xT2';
%	X3=xQ3*ones(1,N)-ones(Q,1)*xT3';
%	distsq=X1.^2+X2.^2+X3.^2;
%	[dummy,nMin]=min(distsq,[],2);
        
    for iMin=1:blockNodes:Q
        iMax=min(iMin+blockNodes-1,Q);
        blockCurr=iMax-iMin+1;
        X1=xQ1(iMin:iMax)*ones(1,N)-ones(blockCurr,1)*xT1';
        X2=xQ2(iMin:iMax)*ones(1,N)-ones(blockCurr,1)*xT2';
        X3=xQ3(iMin:iMax)*ones(1,N)-ones(blockCurr,1)*xT3';
        distsq=X1.^2+X2.^2+X3.^2;
        [~,nMin(iMin:iMax)]=min(distsq,[],2);
    end    
	NClosest=sparse(Q,N); % creates sparse all-zero matrix
	NClosest((1:Q)'+Q*(nMin-1))=1;
    NClosest=kron(speye(3),NClosest);
