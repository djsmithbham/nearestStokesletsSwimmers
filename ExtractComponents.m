function [x1,x2,x3]=ExtractComponents(x)

  N=length(x)/3;
  x1=x(1:N);
  x2=x(N+1:2*N);
  x3=x(2*N+1:3*N);
