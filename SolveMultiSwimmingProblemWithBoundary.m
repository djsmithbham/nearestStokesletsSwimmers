function dz=SolveMultiSwimmingProblemWithBoundary(z,swimmer,boundary,t,epsilon,domain,blockSize,varargin)

    % solves force-free multiple swimmer problem for translational and angular velocity U, Om
    % in the presence of a stationary rigid boundary (e.g. finite plane wall)
    %
    % input:           swimmer{.} - array of structures describing how to construct swimmer - all same size
    %                  (internal variable Nsw is length of swimmer, i.e. number of swimmers, Ns is no. force points on swimmer)
    %        z       - position/orientation of swimmer
    %                  dz(1:3*Nsw)        = x0      - origin of swimmer
    %                  dz(3*Nsw+1:6*Nsw)  = b1      - first basis vector of swimmer frame
    %                  dz(6*Nsw+1:9*Nsw)  = b2      - second basis vector of swimmer frame
    %        t       - time (scalar)
    %        epsilon - regularisation parameter
    %        domain  - 'i' for no image systems, 'h' for Ainley/Blake halfspace
    %        blockSize - memory size to use for stokeslet matrix construction in GB. 0.2 works for any reasonable hardware
    %        varargin  - if varargin{1}='f' then include force components
    %
    % variables: U      - translational velocity
    %            Om     - angular velocity
    %            xb     - discretisation of stationary boundary - force points
    %            Xb     - discretisation of stationary boundary - quadrature points
    %            Is{n}  - index of first point of swimmer n
    %
    % output: dz        - rate of change of position-orientation (x,b1,b2)
    %                     first 3*Nsw compts are dx0/dt, 
    %                     next  3*Nsw are db1/dt,
    %                     next  3*Nsw are db2/dt
    %                     each are ordered by 1-components first, then 2-components, then 3-components
    %                   - if varargin{1}='f' then includes force components
   
    Nsw=length(swimmer);
    
    for n=1:Nsw
        x0{n}=z(n:Nsw:n+2*Nsw);             %3*(n-1)+1:3*n);
        b1{n}=z(3*Nsw+n:Nsw:n+5*Nsw);       %:3*Nsw+3*n);
        b2{n}=z(6*Nsw+n:Nsw:n+8*Nsw);       %+3*(n-1)+1:6*Nsw+3*n);
        b3{n}=cross(b1{n},b2{n});
        B{n}=[b1{n}(:) b2{n}(:) b3{n}(:)];
        [xi{n},vi{n},Xi{n}]=swimmer{n}.fn(t,swimmer{n}.model);
        xx0{n}=ApplyRotationMatrix(B{n},xi{n});
        xsca{n}=TranslatePoints(xx0{n},x0{n});
        vsca{n}=ApplyRotationMatrix(B{n},vi{n});
        Xsca{n}=ApplyRotationMatrix(B{n},Xi{n});
        Xsca{n}=TranslatePoints(Xsca{n},x0{n});
    end
    Is{1}=1;
    if Nsw>1
        for n=1:Nsw
            Is{n+1}=Is{n}+length(xsca{n})/3;
        end
    else
        Is{2}=Is{1}+length(xsca{1})/3;
    end
    
    [xb,Xb]=boundary.fn(boundary.model);
    vb=xb*0; % boundary is stationary
       
    xs=xsca{1};
    Xs=Xsca{1};
    vs=vsca{1};
    if Nsw>1
        for n=2:Nsw
            xs=MergeVectorGrids(xs,xsca{n});
            Xs=MergeVectorGrids(Xs,Xsca{n});
            vs=MergeVectorGrids(vs,vsca{n});
        end
    end
    
    x=MergeVectorGrids(xs,xb);
    X=MergeVectorGrids(Xs,Xb);
    v=MergeVectorGrids(vs,vb);
    
    % now formulate fluid dynamics problem... essentially a mobility problem 

    Ns=length(xs)/3;
    Qs=length(Xs)/3;
    Nb=length(xb)/3;
    Qb=length(Xb)/3;
    
    N=Ns+Nb;
    
    [AS,~]=AssembleStokesletMatrix(x,X,x,epsilon,domain,blockSize);
    NNss=NearestNeighbourMatrix(Xs,xs,blockSize);

%    AU =-kron(eye(3),[ones(Ns,1);zeros(Nb,1)]) ; 
    
    au = zeros(Ns+Nb,Nsw);
    
    for n=1:Nsw
        au(Is{n}:Is{n+1}-1,n) = -ones(Is{n+1}-Is{n},1);
    end
    AU = kron(eye(3),au);

%    AU =-kron(eye(3), ...
%              [kron(eye(Nsw),ones(Ns,1)); zeros(Nb,Nsw)]); % component of velocity due to translational velocity of swimmers; zero velocity of boundary
    
    af = zeros(Nsw,Ns+Nb);
    for n=1:Nsw
        af(n,Is{n}:Is{n+1}-1)=sum(NNss(1:Qs,Is{n}:Is{n+1}-1),1);
    end
    AF = kron(eye(3),af);
          
%    AF = kron(eye(3),[kron(eye(Nsw),sum(NNss(1:Qs,1:Ns),1)) zeros(Nsw,Nb)]); % force summation - only on swimmers
    
    %ze=0*x1; 
%    AOm=[ze -x3 x2; zeros(Nb,3); x3 ze -x1; zeros(Nb,3); -x2 x1 ze; zeros(Nb,3)];

%    Ze = zeros(Nsw*Ns+Nb,Nsw);
    % component of velocity due to rotation of swimmer about x0; zero velocity of boundary

    Ze  = zeros(Ns+Nb,Nsw);
    x1m = zeros(Ns+Nb,Nsw);
    x2m = zeros(Ns+Nb,Nsw);
    x3m = zeros(Ns+Nb,Nsw);
    for n=1:Nsw
        [x1,x2,x3]=ExtractComponents(xx0{n});
        x1m(Is{n}:Is{n+1}-1,n)=x1;   % x1(Is{n}:Is{n+1}-1);  
        x2m(Is{n}:Is{n+1}-1,n)=x2;
        x3m(Is{n}:Is{n+1}-1,n)=x3;
    end    
    AOm = [ Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze];
    
%     AOm=[ Ze                                  [kron(eye(Nsw),-x3); zeros(Nsw*Ns)]  [kron(eye(Nsw), x2); zeros(Nsw*Ns)]; ...
%           [kron(eye(Nsw), x3); zeros(Nsw*Ns)] Ze                                   [kron(eye(Nsw),-x1); zeros(Nsw*Ns)]; ...
%           [kron(eye(Nsw),-x2); zeros(Nsw*Ns)] [kron(eye(Nsw), x1); zeros(Nsw*Ns)]  Ze                                ];    
    
     
    % ze=0*x1;
%     AM=[ ze zeros(1,Nb) -x3 zeros(1,Nb)  x2 zeros(1,Nb); ...
%          x3 zeros(1,Nb)  ze zeros(1,Nb) -x1 zeros(1,Nb); ...
%         -x2 zeros(1,Nb)  x1 zeros(1,Nb)  ze zeros(1,Nb)];  
    
%    Ze = zeros(Nsw,Nsw*Ns+Nb);
   
    Ze  = zeros(Nsw,Ns+Nb);
    x1m = zeros(Nsw,Ns+Nb);
    x2m = zeros(Nsw,Ns+Nb);
    x3m = zeros(Nsw,Ns+Nb);
    [x1,x2,x3]=ExtractComponents(Xs'*NNss); % moment summation
    for n=1:Nsw
        x1m(n,Is{n}:Is{n+1}-1)=x1(Is{n}:Is{n+1}-1);
        x2m(n,Is{n}:Is{n+1}-1)=x2(Is{n}:Is{n+1}-1);
        x3m(n,Is{n}:Is{n+1}-1)=x3(Is{n}:Is{n+1}-1);
    end    
    AM = [Ze -x3m x2m; x3m Ze -x1m; -x2m x1m Ze];
    
%     AM=[ Ze                                 [kron(eye(Nsw),-x3) zeros(Nsw,Nb)]  [kron(eye(Nsw), x2) zeros(Nsw,Nb)]; ...
%          [kron(eye(Nsw), x3) zeros(Nsw,Nb)] Ze                                  [kron(eye(Nsw),-x1) zeros(Nsw,Nb)]; ...
%          [kron(eye(Nsw),-x2) zeros(Nsw,Nb)] [kron(eye(Nsw), x1) zeros(Nsw,Nb)]  Ze                               ];

    A=[AS AU AOm; AF zeros(3*Nsw,6*Nsw); AM zeros(3*Nsw,6*Nsw)];
    
    % rhs assembly
    rhs=[v;zeros(6*Nsw,1)];

    % solve and extract f, U, Omega
    sol=A\rhs;
    f = sol(1           : 3*N);
    U = sol(3*N+1       : 3*N+3*Nsw);
    Om= sol(3*N+3*Nsw+1 : 3*N+6*Nsw);

    dz(1      :3*Nsw ,1)  = U;
    for n=1:Nsw
        om{n}=Om(n:Nsw:n+2*Nsw);
        dz(3*Nsw+n:Nsw:n+5*Nsw,1)  = cross(om{n},b1{n});
        dz(6*Nsw+n:Nsw:n+8*Nsw,1)  = cross(om{n},b2{n});
    end
    
    if ~isempty(varargin)
        dz(9*Nsw+1:9*Nsw+3*N,1)=f;
    end
    
    
