function F=ConstructInterpolantFromSTForm(s,t,stWaveFn,params)

% constructs a griddedInterpolant for flagellar waveforms of (s,t) form
%
% (s,t) are produced by ndgrid
% stWaveFn - function to give wave for vector arclength (s) and scalar time (t)
% params - structure with coefficients defining waveform

x=zeros(size(s,1),size(t,2));
y=zeros(size(s,1),size(t,2));

for k=1:size(t,2)
    [x(:,k),y(:,k),~,~]=stWaveFn(s(:,1),t(1,k),params);
end

F{1}=griddedInterpolant(s,t,x,'spline');
F{2}=griddedInterpolant(s,t,y,'spline');