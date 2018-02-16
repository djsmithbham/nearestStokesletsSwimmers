function [h] = CalcDiscr_hMin(x,y,blockSize)
%CalcDiscr_h This function calculates the minimum point distance between
%points on each module (e.g. head + flagellum)

    Nx=length(x)/3;

    x1=x(1:Nx);
    x2=x(Nx+1:2*Nx);
    x3=x(2*Nx+1:3*Nx);
    
    Ny = length(y)/3;
    y1=y(1:Ny);
    y2=y(Ny+1:2*Ny);
    y3=y(2*Ny+1:3*Ny);
    
    D1 = x1 - y1';
    D2 = x2 - y2';
    D3 = x3 - y3';
    
    distsq = D1.^2 + D2.^2 + D3.^2;
    h = sqrt(min(distsq(:)));
