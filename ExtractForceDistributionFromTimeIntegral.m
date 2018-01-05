function f=ExtractForceDistributionFromTimeIntegral(t,H,t0)

% extracts force distribution at 3N points from the time integral
% produced by ode45

% input:    t,  ode45 time values
%           H,  time-integrated force distribution produced by ode45
%           t0, scalar time value at which force is required
% output:   f,  force distribution at 3N points

dt=1e-5;

DOF=size(H,2);
[T,NP]=ndgrid(t,1:DOF);

T0=repmat(t0,1,DOF);

F=griddedInterpolant(T,NP,H);%,T2,1:DOF
f=(F(T0+dt,1:DOF)-F(T0-dt,1:DOF))/(2*dt);
