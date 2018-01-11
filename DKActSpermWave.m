function outputs=DKActSpermWave(xii,t,outopt,varargin)

% generates 'activated'waveform used by Dresdner & Katz 1981

% if varargin is non-empty, wavenumber is read in from first entry

% ds is arclength metric

a=0.1087;
b=0.0543;

if ~isempty(varargin)
    k=varargin{1}.k;
    phase=varargin{1}.phase;
else
    k=2*pi;
    phase=0;
end

t=t+phase;

% set up in 'wave frame' first
yii=(a*xii+b).*sin(k*xii-t)-b*sin(-t);

% rotation angle into body frame (so flagellum is tangential to head)
th=atan(a*sin(-t)+k*b*cos(-t));

% rotate into body frame
xi= cos(th)*xii+sin(th)*yii;
yi=-sin(th)*xii+cos(th)*yii;

switch outopt
    case 'xy'
        outputs=[xi,yi];
    case 'ds'
        outputs=sqrt(1+(a*sin(k*xii-t)+k*(a*xii+b).*cos(k*xii-t)).^2);
end