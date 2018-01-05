function S=RegStokeslet(x,X,eps)
	% x is a vector of field points:  3*M
	% X is a vector of source points: 3*Q
	% eps is regularisation parameter
	% outputs an array of regularised stokeslets between field and source points
	% blocks are [Sxx, Sxy, Sxz; Syx, Syy, Syz; Szx, Szy, Szz] where Sxx is M x Q etc. 

	x=x(:);
	X=X(:);
	M=length(x)/3;
	Q=length(X)/3;
	r1=      x(1:M)*ones(1,Q)-ones(M,1)*      X(1:Q)';
	r2=  x(M+1:2*M)*ones(1,Q)-ones(M,1)*  X(Q+1:2*Q)';
	r3=x(2*M+1:3*M)*ones(1,Q)-ones(M,1)*X(2*Q+1:3*Q)';
	rsq=r1.^2+r2.^2+r3.^2;
	ireps3=1./(sqrt((rsq+eps^2)).^3);
	isotropic=kron(eye(3),(rsq+2.0*eps^2).*ireps3);
	dyadic=[r1.*r1 r1.*r2 r1.*r3; r2.*r1 r2.*r2 r2.*r3; r3.*r1 r3.*r2 r3.*r3].*kron(ones(3,3),ireps3);

	S=(1.0/(8.0*pi))*(isotropic+dyadic);
