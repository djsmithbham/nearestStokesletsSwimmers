%CONVBOUNDARY Solves the sperm swimming problem for convergence testing
%
%CONVRESULTS = CONVBOUNDARY(F1,F2,F3,F4,NBEATS,BNDRY) solves the sperm
%swimming problem with the number of points in the boundary discretisation 
%specified by F1, F2, F3, F4, over NBEATS of the flagellum.
%
%   -----
%   Input
%   -----
%   F1      - No. of force points in boundary discretisation in x direction
%   F2      - No. of force points in boundary discretisation in y direction
%   F3      - No. of quad. points in boundary discretisation in x direction
%   F4      - No. of quad. points in boundary discretisation in y direction
%   NBEATS  - Number of beats of the flagellum
%
%   ------
%   Output
%   ------
%   CONVRESULTS.POSITION        - Swimmer position in time
%              .B1              - 1st basis vector for sperm frame of ref.
%              .B2              - 2nd basis vector for sperm frame of ref.
%              .FORCES          - Forces on each point at each time
%              .TIME            - Time in multiples of nBeats
%              .EPSILON         - Regularisation parameter 
%              .ns              - No. of force points on flagellum
%              .nh              - No. of force points on cell head = 6*F2^2
%              .nS              - No. of quadrature points on flagellum
%              .nH              - No. of quad. points on cell head = 6*F4^2
%              .hx              - Force disc. spacing in x direction
%              .hy              - Force disc. spacing in y direction
%              .hX            	- Quad. disc. spacing in x direction
%              .hY           	- Quad. disc. spacing in y direction
%              .SOLUTIONTIME    - CPU time to solve
%              .DISTRAVELED     - Straight line distance traveled 
%
% Author: M.T.Gallagher, all rights reserved 2018
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
function convResults = ConvBoundary(f1,f2,f3,f4,nBeats)

% Beats
%nBeats = 1;
tRange = [0, 2*pi*nBeats];

% Time
dt = 2*pi*0.05;

% Domain - 'i' infinite
%        - 'h' use half space images (Ainley/Blakelet)
domain = 'i';

% Memory allocation
blockSize = 0.2;

% Initialise swimmer ------------------------------------------------------
% swimmer no. 1
args.phase=pi/5;
args.k=2*pi;
nns=30;
s = linspace(0,1,nns);
t = tRange(1):dt:tRange(2);
[S,T]=ndgrid(s,t);
x00=[0;0;0.2];
B=RotationMatrix(0*pi/3,3);
b10=B(:,1);
b20=B(:,2);
xyWaveFn=@DKActSpermWave;
swimmer.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);
swimmer.fn=@SpermModelGI;

% semi-axes
swimmer.model.a1=2.0/45;
swimmer.model.a2=1.6/45;
swimmer.model.a3=1.0/45;
swimmer.model.ns=40;
swimmer.model.nh=4;
% Fine grid
swimmer.model.Ns=100;
swimmer.model.Nh=10;

% boundary ---------------------
boundary.fn=@PlaneBoundary2;
boundary.model.h=0.4;
boundary.model.Lx=3;
boundary.model.Ly=3;
boundary.model.O=[0 0 0];

%% Boundary discretisation convergence
epsilon=0.25/45;

% Discretisation parameters for boundary
% Coarse grid
boundary.model.nx=f1;
boundary.model.ny=f2;
% Fine grid
boundary.model.Nx=f3;
boundary.model.Ny=f4;

tic
[t,z]=SolveSwimmingTrajectoryAndForces(x00,b10,b20, ...
    tRange,swimmer,boundary,epsilon,domain, blockSize);
solTime = toc;
fprintf('time taken = %f\n',toc)

convResults.position = z(:,1:3);
convResults.b1 = z(:,4:6);
convResults.b2 = z(:,7:9);
convResults.forces = z(:,10:end);
convResults.time = t/2/pi;
convResults.epsilon = epsilon;
convResults.ns = swimmer.model.ns;
convResults.nh = swimmer.model.nh;
convResults.nS = swimmer.model.Ns;
convResults.nH = swimmer.model.Nh;
convResults.nx = boundary.model.nx;
convResults.ny = boundary.model.ny;
convResults.nX = boundary.model.Nx;
convResults.nY = boundary.model.Ny;
% Discretisations
convResults.hx = boundary.model.Lx/(boundary.model.nx-1);
convResults.hy = boundary.model.Ly/(boundary.model.ny-1);
convResults.hX = boundary.model.Lx/(boundary.model.Nx-1);
convResults.hY = boundary.model.Ly/(boundary.model.Ny-1);
convResults.solutionTime = solTime;
convResults.distTraveled = ...
    norm(convResults.position(end,:) - ...
    convResults.position(1,:));

disp(convResults)

end