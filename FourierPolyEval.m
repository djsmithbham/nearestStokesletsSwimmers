function y=FourierPolyEval(a0,a,b,sigma,x)

% evaluates Fourier series in sine/cosine form

y=a0;

Nf=length(a);

for n=1:Nf
    y=y+a(n)*cos(n*sigma*x)+b(n)*sin(n*sigma*x);
end

