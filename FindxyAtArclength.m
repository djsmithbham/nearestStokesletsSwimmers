function [x,y]=FindxyAtArclength(s,t,xyWaveFn,varargin)

% for waveforms not prescribed in terms of arclength, this function
% takes arclength coordinate input and outputs x,y coordinates at that
% arclength

if isempty(varargin) 
    Resi=@(xii) s-integral(@(xiii) xyWaveFn(xiii,t,'ds'),0,xii);
else
    Resi=@(xii) s-integral(@(xiii) xyWaveFn(xiii,t,'ds',varargin{1}),0,xii);
end

options=optimoptions(@fsolve,'Display','off');

xii0=s;


xii=fsolve(Resi,xii0,options);

if isempty(varargin)
    outputs=DKActSpermWave(xii,t,'xy');
else
    outputs=DKActSpermWave(xii,t,'xy',varargin{1});
end

x=outputs(:,1);
y=outputs(:,2);

