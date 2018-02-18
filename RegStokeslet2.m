function S=RegStokeslet2(x,X,idx,ep)
	% x is a 3N vector of field points
    % X is a 3Q vector of source points
    % idx{m} is an index of which source points to evaluate for x(m)
    % ep is regularisation parameter
	% outputs a possibly sparse array of regularised stokeslets between field and source points
	% blocks are [Sxx, Sxy, Sxz; Syx, Syy, Syz; Szx, Szy, Szz] where Sxx is N x Q etc. 

    x=x(:);
    X=X(:);
    
	N=length(x)/3;
	Q=length(X)/3;
%    inz=find(xmat(1:N,1:Q));
    
    %xmat(151)
    %xmat(151,1)
    
%     r1=sparse(N,Q);
%     r2=sparse(N,Q);
%     r3=sparse(N,Q);
%     rsq=sparse(N,Q);
%     irep3=sparse(N,Q);
%     iso=sparse(N,Q);
    
    counter=0;
    for m=1:N
        clear r1 r2 r3 rsq irep3 iso
        
        nnz=length(idx{m});
        
        r1=x(m)-X(idx{m});
        r2=x(N+m)-X(Q+idx{m});
        r3=x(2*N+m)-X(2*Q+idx{m});
        rsq=r1.^2+r2.^2+r3.^2;
        irep3=(rsq+ep^2).^(-3/2);
        iso=(rsq+2.0*ep^2);
        
        v11(counter+1:counter+nnz)=(iso+r1.*r1).*irep3/8.0/pi;
        v12(counter+1:counter+nnz)=(    r1.*r2).*irep3/8.0/pi;
        v13(counter+1:counter+nnz)=(    r1.*r3).*irep3/8.0/pi;
        v21(counter+1:counter+nnz)=(    r2.*r1).*irep3/8.0/pi;
        v22(counter+1:counter+nnz)=(iso+r2.*r2).*irep3/8.0/pi;
        v23(counter+1:counter+nnz)=(    r2.*r3).*irep3/8.0/pi;
        v31(counter+1:counter+nnz)=(    r3.*r1).*irep3/8.0/pi;
        v32(counter+1:counter+nnz)=(    r3.*r2).*irep3/8.0/pi;
        v33(counter+1:counter+nnz)=(iso+r3.*r3).*irep3/8.0/pi;
 
        irow(counter+1:counter+nnz)=m;
        icol(counter+1:counter+nnz)=idx{m};
        counter=counter+nnz;
    end    
    
    S11=sparse(irow,icol,v11,N,Q);
    S12=sparse(irow,icol,v12,N,Q);
    S13=sparse(irow,icol,v13,N,Q);
    S21=sparse(irow,icol,v21,N,Q);
    S22=sparse(irow,icol,v22,N,Q);
    S23=sparse(irow,icol,v23,N,Q);
    S31=sparse(irow,icol,v31,N,Q);
    S32=sparse(irow,icol,v32,N,Q);
    S33=sparse(irow,icol,v33,N,Q);
        
    S=[S11 S12 S13; S21 S22 S23; S31 S32 S33];
    
%     machEps=xmat(1:N,1:Q)*eps;
%     epMat=(~xmat(1:N,1:Q)==0)*ep;
% 	r1=             xmat(1:N,1:Q) -             Xmat(1:N,1:Q) + machEps;    
% 	r2=     xmat(N+1:2*N,Q+1:2*Q) -     Xmat(N+1:2*N,Q+1:2*Q) + machEps;  
% 	r3= xmat(2*N+1:3*N,2*Q+1:3*Q) - Xmat(2*N+1:3*N,2*Q+1:3*Q) + machEps;
%     rsq=r1.^2+r2.^2+r3.^2;
% %    irep3=sparse(N,Q);
%     irep3=(rsq+epMat).^(-3/2);
%     
%     irep3(151)
%     
%     rsq(151)
%     epMat(151)
%     
% %    iso=sparse(N,Q);
%     iso=(rsq+2.0*epMat).*irep3;
%     
    
