%CONVSPERM Solves the sperm swimming problem for convergence testing
%
%CONVRESULTS = CONVSPERM(F1,F2,F3,F4,NBEATS,BNDRY) solves the sperm
%swimming problem with the number of points in the discretisation specified
%by F1, F2, F3, F4, over NBEATS of the flagellum with or without boundary
%BNDRY
%
%   -----
%   Input
%   -----
%   F1      - Number of force points on flagellum
%   F2      - Number of force points on cell head = 6*F2^2
%   F3      - Number of quadrature points on flagellum
%   F4      - Number of quadrature points on cell head = 6*F4^2
%   NBEATS  - Number of beats of the flagellum
%   BNDRY   - TRUE if boundary included
%           - FALSE otherwise
%   ------
%   Output
%   ------
%   CONVRESULTS.POSITION        - Swimmer position in time
%              .B1              - 1st basis vector for sperm frame of ref.
%              .B2              - 2nd basis vector for sperm frame of ref.
%              .FORCES          - Forces on each point at each time
%              .TIME            - Time in multiples of nBeats
%              .EPSILON         - Regularisation parameter 
%              .ns              - F1
%              .nh              - F2
%              .nS              - F3
%              .nH              - F4
%              .HF              - Max. dist. between force points
%              .HQ              - Max. dist. between quadrature points
%              .F_HF            - Max. dist. between flagellum force points
%              .F_HQ            - Max. dist. between flagellum quad. points
%              .H_HF            - Max. dist. between head force points
%              .H_HQ            - Max. dist. between head quadrature points
%              .F_MIN           - Min. dist. between force points
%              .Q_MIN           - Min. dist. between quadrature points
%              .SOLUTIONTIME    - CPU time to solve
%              .DISTRAVELED     - Straight line distance traveled 
%
% Author: M.T.Gallagher, all rights reserved 2018
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
function convResults = ConvSperm(f1,f2,f3,f4,nBeats,bndry)

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

% boundary ---------------------
if bndry
    boundary.fn=@PlaneBoundary2;
    boundary.model.h=0.4;
    boundary.model.nx=16;
    boundary.model.ny=15;
    boundary.model.Nx=32;
    boundary.model.Ny=30;
    boundary.model.Lx=3;
    boundary.model.Ly=3;
    boundary.model.O=[0 0 0];
else
    boundary = [];
end

%% Boundary discretisation convergence
epsilon=0.25/45;

% Discretisation parameters for sperm
% Coarse grid
swimmer.model.ns=round(f1);
swimmer.model.nh=round(f2);
% Fine grid
swimmer.model.Ns=round(f3);
swimmer.model.Nh=round(f4);

tic
[t,z]=SolveSwimmingTrajectoryAndForces(x00,b10,b20, ...
    tRange,swimmer,boundary,epsilon,domain, blockSize);
solTime = toc;
fprintf('time taken = %f\n',toc)

% Calculate hf and hq
[hf,hq,f_hf,f_hq,h_hf,h_hq,f_min,q_min] ...
    = CalcDiscSpacingsSperm(z,t,swimmer,boundary);

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
convResults.hf = max(hf);
convResults.hq = max(hq);
convResults.f_hf = f_hf;
convResults.f_hq = f_hq;
convResults.h_hf = h_hf;
convResults.h_hq = h_hq;
convResults.f_min = f_min;
convResults.q_min = q_min;
convResults.solutionTime = solTime;
convResults.distTraveled = ...
    norm(convResults.position(end,:) - ...
    convResults.position(1,:));
disp(convResults)

end