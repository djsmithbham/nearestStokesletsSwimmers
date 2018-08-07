function [A,NN]=AssembleStokesletMatrix(xCollNodes,xQuadNodes,xForceNodes,eps,domain,blockSize,varargin)

	% blockSize is in GB
	% stokeslet values are calculated in blocks to avoid memory overflow
	% each block is multiplied by the nearest neighbour matrix 
	% the result is summed to yield the stokeslet matrix
    
    if ~isempty(varargin)
        % supplied nearest-neighbour matrix
        NN=varargin{1};
    else
        % find closest force node to each quadrature node
        NN=NearestNeighbourMatrix(xQuadNodes,xForceNodes,blockSize);
    end
    
	M=length(xCollNodes)/3;
	N=length(xForceNodes)/3;
	Q=length(xQuadNodes)/3;

	% calculate number of quadrature nodes that can be used for each block
	blockNodes=floor(blockSize*2^27/(9*M));

	%assemble stokeslet matrix 
	A=zeros(3*M,3*N);
	for iMin=1:blockNodes:Q
		iMax=min(iMin+blockNodes-1,Q);
		iRange=[iMin:iMax Q+iMin:Q+iMax 2*Q+iMin:2*Q+iMax];
                switch domain
                       case 'i'
                           A=A+RegStokeslet(xCollNodes,xQuadNodes(iRange),eps)*NN(iRange,:);
                       case 'h'
  		       	           A=A+RegBlakelet(xCollNodes,xQuadNodes(iRange),eps)*NN(iRange,:);
                end
	end

