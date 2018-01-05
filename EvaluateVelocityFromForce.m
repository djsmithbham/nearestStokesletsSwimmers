function u=EvaluateVelocityFromForce(xf,X,x,f,epsilon,domain,blockSize)

% input:
% xf - field points
% X  - quadrature grid
% x  - force grid
% f  - forces
%
% output:
% u  - velocity field at xf

[A,~]=AssembleStokesletMatrix(xf,X,x,epsilon,domain,blockSize);

u=A*f;
