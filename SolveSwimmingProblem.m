function dz=SolveSwimmingProblem(z,swimmer,t,epsilon,domain,blockSize)

    % solves force-free swimming problem for translational and angular velocity U, Om
    %
    % input: z       - position/orientation of swimmer
    %                  dz(1:3) = x0      - origin of swimmer
    %                  dz(4:6) = b1      - first basis vector of swimmer frame
    %                  dz(7:9) = b2      - second basis vector of swimmer frame
    %        swimmer - structure describing how to construct swimmer
    %        t       - time (scalar)
    %
    % variables: U      - translational velocity
    %            Om     - angular velocity
    %
    % output: dz        - rate of change of position-orientation (x,b1,b2)
    %                     first three compts are dx0/dt, 
    %                     next three are db1/dt,
    %                     next three are db2/dt

    [xi,vi,Xi]=swimmer.fn(t,swimmer.model);
    x0=z(1:3);
    b1=z(4:6);
    b2=z(7:9);
    b3=cross(b1,b2);
    B=[b1(:) b2(:) b3(:)];
    xx0=ApplyRotationMatrix(B,xi); % x-x0
    x=TranslatePoints(xx0,x0);
    v=ApplyRotationMatrix(B,vi);
    X=ApplyRotationMatrix(B,Xi);
    X=TranslatePoints(X,x0);

    % now formulate fluid dynamics problem... essentially a mobility problem 

    N=length(x)/3;
    Q=length(X)/3;

    % lhs assembly
    [AS,NN]=AssembleStokesletMatrix(x,X,x,epsilon,domain,blockSize);

    AF=[sum(NN(1:Q,:),1);sum(NN(Q+1:2*Q,:),1);sum(NN(2*Q+1:3*Q,:),1)]; % force summation
    AU=-kron(eye(3),ones(N,1)); % component of velocity due to translational velocity of frame

    [x1 x2 x3]=ExtractComponents(xx0);ze=0*x1; % component of velocity due to rotation of frame about x0
    AOm=[ze -x3 x2; x3 ze -x1; -x2 x1 ze];

    [x1 x2 x3]=ExtractComponents(X'*NN);ze=0*x1; % moment summation
    AM=[ze -x3 x2; x3 ze -x1; -x2 x1 ze];

    A=[AS AU AOm; AF zeros(3,6); AM zeros(3,6)];

    % rhs assembly
    rhs=[v;zeros(6,1)];

    % solve and extract f, U, Omega
    sol=A\rhs;
    f=sol(1:3*N);
    U=sol(3*N+1:3*N+3);
    Om=sol(3*N+4:3*N+6);

    dz(1:3,1)=U;
    dz(4:6,1)=cross(Om,b1);
    dz(7:9,1)=cross(Om,b2);
    