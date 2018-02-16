%CONVCHLAMY Solves the chlamy swimming problem for convergence testing
%
%CONVRESULTS = CONVCHLAMY(F1,F2,F3,F4,NBEATS,BNDRY) solves the chlamy
%swimming problem with the number of points in the discretisation specified
%by F1, F2, F3, F4, over NBEATS of the flagella
%
%   -----
%   Input
%   -----
%   F1      - Number of force points on flagella
%   F2      - Number of force points on cell head = 6*F2^2
%   F3      - Number of quadrature points on flagella
%   F4      - Number of quadrature poitns on cell head = 6*F4^2
%   NBEATS  - Number of beats of the flagella
%
%   ------
%   Output
%   ------
%   CONVRESULTS.POSITION        - Swimmer position in time
%              .B1              - 1st basis vector for chlamy frame of ref.
%              .B2              - 2nd basis vector for chlamy frame of ref.
%              .FORCES          - Forces on each point at each time
%              .TIME            - Time in multiples of nBeats
%              .EPSILON         - Regularisation parameter 
%              .ns              - F1
%              .nh              - F2
%              .nS              - F3
%              .nH              - F4
%              .HF              - Max. dist. between force points
%              .HQ              - Max. dist. between quadrature points
%              .F_HF            - Max. dist. between flagella force points
%              .F_HQ            - Max. dist. between flagella quad. points
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
function convResults = ConvChlamy(f1,f2,f3,f4,nBeats)

% Beats
%nBeats = 1;
tRange = [0, 2*pi*nBeats];

% Time
dt = 2*pi*0.05;
t = tRange(1) : dt : tRange(2);


% Initial position vector
x00 = [0; 0; 0];

% Initial rotation matrix
B = RotationMatrix(0*pi/3, 3);

% Initial basis vectors (3rd vector is x-prod)
b10 = B(:,1);
b20 = B(:,2);

% Generate gridded interpolant for cilium --------------------------------
% Arclength discretisation for prescribed waveform
nS = 30;
s = linspace(0, 1, nS);

[S,T] = ndgrid(s,t);

stWaveFn=@ChlamyFromModel;
swimmer.model.F=ConstructInterpolantFromSTForm(S,T,stWaveFn,[]);
%--------------------------------------------------------------------------

% Domain - 'i' infinite
%        - 'h' use half space images (Ainley/Blakelet)
domain = 'i';

% Memory allocation
blockSize = 0.2;

% Initialise swimmer ------------------------------------------------------
% Swimmer function
swimmer.fn = @ChlamyModel1;

% Discretisation parameters
% Coarse grid
swimmer.model.ns=40;
swimmer.model.nh=4;
% Fine grid
swimmer.model.Ns=100;
swimmer.model.Nh=10;

% Chlamy body semi-axes
swimmer.model.a1=0.5;
swimmer.model.a2=0.6;
swimmer.model.a3=0.6;

% chlamy flagellar angle
swimmer.model.ang=pi/6;

%% flagellar discretisation convergence
epsilon = 0.0001;

% Discretisation parameters
% Coarse grid
swimmer.model.ns=round(f1);
swimmer.model.nh=round(f2);
% Fine grid
swimmer.model.Ns=round(f3);
swimmer.model.Nh=round(f4);

tic
[t,z]=SolveSwimmingTrajectoryAndForces(x00,b10,b20, ...
    tRange,swimmer,[],epsilon,domain, blockSize);
solTime = toc;
fprintf('time taken = %f\n',toc)

% Calculate hf and hq
[hf,hq,f_hf,f_hq,h_hf,h_hq,f_min,q_min] ...
    = CalcDiscSpacingsChlamy(z,t,swimmer,[]);

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