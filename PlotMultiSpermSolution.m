function PlotMultiSpermSolution(fig,swimmer,plotSwitch,t,z,m,varargin)

if ~isempty(varargin)
    mksz=varargin{1};
else
    mksz=1;
end

Nsw=length(swimmer);

if Nsw==1
    swimmertemp=swimmer;
    clear swimmer;
    swimmer{1}=swimmertemp;
end

for n=1:Nsw
    x0{n}=z(m,       n:Nsw:n+2*Nsw);             %3*(n-1)+1:3*n);
    b1{n}=z(m, 3*Nsw+n:Nsw:n+5*Nsw);       %:3*Nsw+3*n);
    b2{n}=z(m, 6*Nsw+n:Nsw:n+8*Nsw);       %+3*(n-1)+1:6*Nsw+3*n);
    b3{n}=cross(b1{n},b2{n});
    B{n}=[b1{n}(:) b2{n}(:) b3{n}(:)];
    O{n}=x0{n}; %[x0{n}(1) x0{n}(2) x0{n}(3)];
    n
    t
    [xi{n},vi{n},Xi{n}]=swimmer{n}.fn(t(m),swimmer{n}.model);
    x{n}=ApplyRotationMatrix(B{n},xi{n});
    X{n}=ApplyRotationMatrix(B{n},Xi{n});
    x{n}=TranslatePoints(x{n},O{n});
    X{n}=TranslatePoints(X{n},O{n});
    [x1{n},x2{n},x3{n}]=ExtractComponents(x{n});
    [X1{n},X2{n},X3{n}]=ExtractComponents(X{n});
end

figure(fig);hold on;
for n=1:Nsw
    switch plotSwitch
        case 'f' %force points
            hp=plot3(x1{n},x2{n},x3{n},'b.');set(hp,'markersize',mksz);
        case 'q' %quadrature points
            hp=plot3(X1{n},X2{n},X3{n},'r.');set(hp,'markersize',mksz);
        case 'a' %force and quadrature points
            hp=plot3(x1{n},x2{n},x3{n},'b.');set(hp,'markersize',mksz);
            hp=plot3(X1{n},X2{n},X3{n},'r.');set(hp,'markersize',mksz);
    end
end
