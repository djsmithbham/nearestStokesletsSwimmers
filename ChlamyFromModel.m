%CHLAMYFROMMODEL Generates a chlamy flagellum shape from the model
%   PSI(S,T) = PSI_0(S) + PSI_1(S)*SIN(WT + PHI(S))
%
%   -----
%   Input
%   -----
%   S - Arclength
%   T - Time
%
%   ------
%   Output
%   ------
%   (X,Z)   - Position of flagellum
%   (DX,DZ) - Derivatives in X,Z directions
%
%   ---------------
%   Child functions
%   ---------------   
%   NONE
%
% Author: M.T.Gallagher, all rights reserved 2018
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
function [x,z,dx,dz] = ChlamyFromModel(s,t,params) %#ok<INUSD> 
                                                      % Ignore params input

[ns,nt] = size(s);

% Generate points at higher resolution for integration
NS = 2*ns + 1;
S = linspace(min(s(:)),max(s(:)),NS);
S = repmat(S(:),1,nt);
T = repmat(t(1,:),NS,1);

% Model components
psi0 = linspace(0,-2.5,NS)';
psi0 = repmat(psi0,1,nt);

psi1 = 0.7 + 0.15*sin(2*pi * S);

phi = linspace(0,-2*pi,NS)';
phi = repmat(phi,1,nt);

psi = psi0 - psi1 .* cos(T + phi);

% Integrate
CX = zeros(ns,nt); CY = CX;
for ii = 1:nt
    % Step size
    ds = S(3,1)-S(1,1);
    
    % Coefficient matrix
    iO = 1:ns;     iO = [iO  , iO, iO  ];
    jO = 2*(1:ns); jO = [jO-1, jO, jO+1];
    aO = [1/6*ones(1,ns), 2/3*ones(1,ns), 1/6*ones(1,ns)].*ds;
    A = sparse(iO,jO,aO,ns,2*ns+1);
    
    intX = cos(psi(:,ii));
    intY = sin(psi(:,ii));
    
    % Integrate
    X = A*intX;
    Y = A*intY;
    
    % Sum
    X = cumsum(X);
    
    Y = cumsum(Y);
    
    CX(:,ii) = X;
    CY(:,ii) = Y;
end

CX = CX - CX(1,:);
CY = CY - CY(1,:);

% Rotate
RTh =pi/2 - atan2(CY(5,1),CX(5,1));
x = cos(RTh)*CX - sin(RTh)*CY;
z = sin(RTh)*CX + cos(RTh)*CY;
z = -z;

% Derivatives from finite differences
dx = nan;
dz = nan;






