function [t,z]=SolveMultiSwimmingTrajectoryAndForces(x00,b10,b20,tRange,swimmer,boundary,epsilon,domain,blockSize,varargin)

% Solves problem of multi swimmer trajectory and accompanying time-varying forces
% for velocity field reconstruction
%
% input:
%   x00        Initial position vector
%   b10, b20   Initial basis vectors (third basis vector is calculated via
%              cross product)
%   tRange     time range over which to calculate trajectory
%   swimmer    Structure describing how to calculate swimmer shape and
%              kinematics
%   boundary   Structure describing boundary
%   epsilon    regularisation parameter
%   domain     switch for use of image systems
%              'i' = use infinite fluid (usual)
%              'h' = use half space images (Ainley/Blakelet)
%   blockSize  memory size for stokeslet matrix block assembly in GB.
%              0.2 is a safe choice.
%   varargin   if varargin{1}='f' then includes force components

Nsw=length(swimmer);
if Nsw==1
    swimmertemp=swimmer;
    clear swimmer;
    swimmer{1}=swimmertemp;
    x00temp=x00;
    clear x00;
    x00{1}=x00temp;
    b10temp=b10;
    clear b10;
    b10{1}=b10temp;
    b20temp=b20;
    clear b20;
    b20{1}=b20temp;
end
xis=[];
for n=1:Nsw
    [xistemp,~,~] = swimmer{n}.fn(0,swimmer{n}.model);
    xis=[xis; xistemp];
end
[xib,~]   = boundary.fn(boundary.model);
DOF=length(xis)+length(xib);

z0=zeros(9*Nsw+DOF,1);
for n=1:Nsw
    z0(n)      = x00{n}(1);
    z0(n+Nsw)  = x00{n}(2);
    z0(n+2*Nsw)= x00{n}(3);
    z0(3*Nsw+n)= b10{n}(1);
    z0(4*Nsw+n)= b10{n}(2);
    z0(5*Nsw+n)= b10{n}(3);
    z0(6*Nsw+n)= b20{n}(1);
    z0(7*Nsw+n)= b20{n}(2);
    z0(8*Nsw+n)= b20{n}(3);
end

%z0=[x00;b10;b20;zeros(DOF,1)];

[t,z]=ode45(@(t,z) SolveMultiSwimmingProblemWithBoundary(z,swimmer,boundary,t,epsilon,domain,blockSize,varargin),tRange,z0);
