function [t,z]=SolveSwimmingTrajectoryAndForces(x00,b10,b20,tRange,swimmer,boundary,epsilon,domain,blockSize,varargin)

% Solves problem of swimmer trajectory and accompanying time-varying forces
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

[xis,~,~] = swimmer.fn(0,swimmer.model);
if isempty(boundary)
    xib=[];
else
    [xib,~]   = boundary.fn(boundary.model);
end

DOF=length(xis)+length(xib);
z0=[x00;b10;b20;zeros(DOF,1)];
[t,z]=ode45(@(t,z) SolveSwimmingProblemWithBoundary(z,swimmer,boundary,t,epsilon,domain,blockSize,varargin),tRange,z0);
