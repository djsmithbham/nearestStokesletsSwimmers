function PlotMultiSwimmerVelocityField(Xg,Yg,Zg,t,z,m,swimmer,boundary,epsilon,domain,blockSize,varargin)

if ~isempty(varargin)
    scl=varargin{1};
else
    scl=1;
end

% plotting grid
xf=[Xg(:);Yg(:);Zg(:)];

% check how many swimmers and extract kinematics
Nsw=length(swimmer)

if Nsw==1
    swimmertemp=swimmer;
    clear swimmer;
    swimmer{1}=swimmertemp;
end
for n=1:Nsw
    x0{n}=z(m,       n:Nsw:n+2*Nsw);
    b1{n}=z(m, 3*Nsw+n:Nsw:n+5*Nsw);
    b2{n}=z(m, 6*Nsw+n:Nsw:n+8*Nsw);
    b3{n}=cross(b1{n},b2{n});
    B{n}=[b1{n}(:) b2{n}(:) b3{n}(:)];
    O{n}=x0{n};
    [xi{n},vi{n},Xi{n}]=swimmer{n}.fn(t(m),swimmer{n}.model);
    x{n}=ApplyRotationMatrix(B{n},xi{n});
    X{n}=ApplyRotationMatrix(B{n},Xi{n});
    x{n}=TranslatePoints(x{n},O{n});
    X{n}=TranslatePoints(X{n},O{n});
    [x1{n},x2{n},x3{n}]=ExtractComponents(x{n});
    [X1{n},X2{n},X3{n}]=ExtractComponents(X{n});
end

% merge with boundary to get force points
if isempty(boundary)
    xb=[];
    Xb=[];
else
    [xb,Xb]=boundary.fn(boundary.model);
end
xs=x{1};
Xs=X{1};
if Nsw>1
    for n=2:Nsw
        xs=MergeVectorGrids(xs,x{n});
        Xs=MergeVectorGrids(Xs,X{n});
    end
end
x=MergeVectorGrids(xs,xb);
X=MergeVectorGrids(Xs,Xb);

% extract force from integral
H=z(:,9*Nsw+1:end);
f=ExtractForceDistributionFromTimeIntegral(t,H,t(m)); 
f=f(:);

% evaluate and plot velocity field on grid
u=EvaluateVelocityFromForce(xf,X,x,f,epsilon,domain,blockSize);

[Ug,Vg,Wg]=ExtractComponents(u);
quiver(Xg(:),Yg(:),scl*Ug,scl*Vg,0);
