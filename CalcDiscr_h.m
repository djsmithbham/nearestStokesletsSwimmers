function [h] = CalcDiscr_h(x,blockSize)
%CalcDiscr_h This function calculates the maximum over all points in a discretisation
%           x          of the distance to the nearest neighbour point
%           blockSize  max size of block of points to work with, measured
%                      in GB - to prevent memory overrun

    N=length(x)/3;

    x1=x(1:N);
    x2=x(N+1:2*N);
    x3=x(2*N+1:3*N);
    
% calculate number of nodes that can be used for each block
	blockNodes=floor(blockSize*2^27/(9*N));
    
    hMin=zeros(N,1);    
    for iMin=1:blockNodes:N
        iMax=min(iMin+blockNodes-1,N);
        blockCurr=iMax-iMin+1;        
    	X1=x1(iMin:iMax)*ones(1,N)-ones(blockCurr,1)*x1';
        X2=x2(iMin:iMax)*ones(1,N)-ones(blockCurr,1)*x2';
        X3=x3(iMin:iMax)*ones(1,N)-ones(blockCurr,1)*x3';
        distsq=X1.^2+X2.^2+X3.^2+[zeros(blockCurr,iMin-1) 100*eye(blockCurr) zeros(blockCurr,N-blockCurr-iMin+1)];
        [hMin(iMin:iMax),~]=min(distsq,[],2); 
    end
    h=sqrt(max(hMin));
