%NNSWIMMERTEMPLATE Template for the Nearest-Neighbour swimmer code
%
%% Generate swimmer
% The swimmer should be a struct (or cell array of structs for the case of
% multiple swimmers) with the following form
%
%   swimmer.model   - Should contain all information to calculate the
%  (swimmer{ii}.model)  discretisation of the swimmer
%                       Typically this would include:
%                           - Number of force points for head
%                           - Number of quadrature points for head
%                           - Number of force points for flagellum
%                           - Number of quadrature points for flagellum
%                           - Dimensions of cell (e.g. head axis lengths,
%                               flagellum length, etc)
%                           - F - gridded interpolant for cell position at
%                               time t (in body-frame coordinates)
%
%   swimmer.fn      - Function for swimmer discretisation
%  (swimmer{ii}.fn)        Inputs should be:
%                               t     - time 
%                               model - swimmer.model
%                       Outputs should be:
%                               x     - Force points for swimmer
%                               v     - Velocities at force points
%                                       (calculated by e.g. finite
%                                       differences
%                               X     - Quadrature points for swimmer
%
%
% Swimmers also need the following associated information:
% x00 (x00{ii}) - The origin of the swimmer at t = 0
% b10 (b10{ii}) - The first component of the swimmer basis vector
% b20 (b20{ii}) - The second component of the swimmer basis vector

%% Generate Boundary
% The boundary can either be empty ( = [] ) or contain the following
%   information
%
% boundary.model    - Contain all information to calculate the
%                       discretisation of the boundary
%                       Typically this would include:
%                           - Number of force points in each direction
%                           - Number of quadrature points in each direction
%                           - Boundary lengths in each direction
%                           - Boundary separations in the case of multiple
%                               boundaries
%
% boundary.fn       - Function for boundary discretisation
%                       Inputs should be:
%                               boundary.model
%
%                       Outputs should be:
%                               x    - Force points on boundary
%                               X    - Quadrature points on boundary
%
%% Numerical parameters
% Here you need to specify:
%   epsilon    - The regularisation parameter
%   domain     - 'i' for infinite fluid
%                'h' to use half space images (Ainley / Blakelet)
%   tRange     - Time range over which to calculate trajectory
%   blockSize  - Memory size for stokeslet matrix block assembly in GB.
%                    0.2 is a safe choice.
%
%% Solve the swimming problem
%
% The swimming problem can then be solved using either:
%
%   [t,z]=SolveSwimmingTrajectoryAndForces(x00,b10,b20,tRange,swimmer, ...
%           boundary,epsilon,domain,blockSize)
%
% for a single swimmer, or
%
%   [t,z]=SolveMultiSwimmingTrajectoryAndForces(x00,b10,b20,tRange, ...
%           swimmer,boundary,epsilon,domain,blockSize)
%
%   
% The outputs are:
%   t - Time steps the solution is calculated at 
%   z - Solution with the following form:
%           z(t,1:3)    - Position of swimmer at time t in x_1, x_2, x_3
%           z(t,4:6)    - Swimmer first basis vector at time t
%           z(t,7:9)    - Swimmer second basis vector at time t
%           z(t,10:end) - Forces calculated at each force point for each
%                           point in time t