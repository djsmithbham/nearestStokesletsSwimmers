function [t,z]=odesolver(fun,tRange,z0)

dt=0.01;
N=ceil((tRange(end)-tRange(1))/dt);

t=linspace(tRange(1),tRange(end),N);

z(1,:)=z0;
for n=1:N-1
    z(n+1,:)=z(n,:)+fun(t(n),z(n,:))';
end
