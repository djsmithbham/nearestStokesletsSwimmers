%GENERATESWIMMINGFIGURESPERM Generates figure 2 from the manuscript
%
function GenerateSwimmingFigureSperm

% Font options for printing figures
fs=8;
fn='times';
wd=7.0;
ht=6.0;

%% Panel a and c  - image of sperm beat pattern
nbeats=1; 
tRange=[0 2*pi*nbeats];
nns=30;

dt=2*pi*0.05;

% Generate swimmer ----------------------------------
% waveform
xyWaveFn=@DKActSpermWave;
args.phase=0;
args.k=2*pi;
s = linspace(0,1,nns);
t = tRange(1):dt:tRange(2);
[S,T]=ndgrid(s,t);
x00{1}=[0;0;0.2];
B{1}=RotationMatrix(0*pi/3,3);
b10{1}=B{1}(:,1);
b20{1}=B{1}(:,2);
swimmer{1}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);
swimmer{1}.fn=@SpermModelGI;

% discretisation parameters - number of points
swimmer{1}.model.ns=40;
swimmer{1}.model.nh=4;
swimmer{1}.model.Ns=100;
swimmer{1}.model.Nh=10;
% semi-axes
swimmer{1}.model.a1=2.0/45;
swimmer{1}.model.a2=1.6/45;
swimmer{1}.model.a3=1.0/45;

% Panel a ----------------------------------
figure(1);clf;hold on;
nth=10;
nphi=20;
a=1;
[x3,x2,x1]=GenerateSphereSurfaceForVisualisation(nth,nphi,a);
x1=x1*swimmer{1}.model.a1;
x2=x2*swimmer{1}.model.a2;
x3=x3*swimmer{1}.model.a3;

surf(x1,x2,x3,0*x3);shading flat;light;

clrs = summer(length(t));
for nt=1:length(t)
    [xis,~,~]  = swimmer{1}.fn(t(nt),swimmer{1}.model);
    [x1,x2,x3] = ExtractComponents(xis);
    plot3(x1,x2,x3,'.','color',clrs(nt,:));
end
axis equal;axis off;
set(1,'paperunits','centimeters');
set(1,'papersize',[wd ht]);
set(1,'paperposition',[0 0 wd ht]);
print(1,'-dpng','-r600','figureSperm_a.png');


% Panel c ----------------------------------
figure(2);clf;hold on;

[x3,x2,x1]=GenerateSphereSurfaceForVisualisation(nth,nphi,a);
x1=x1*swimmer{1}.model.a1;
x2=x2*swimmer{1}.model.a2;
x3=x3*swimmer{1}.model.a3;

surf(x1,x2,x3,0*x3);shading flat;light;

clrs = spring(2*length(t)); clrs = clrs(round(length(t)/2)+1:end,:);
for nt=1:length(t)
    [~,~,Xis]  = swimmer{1}.fn(t(nt),swimmer{1}.model);
    [X1,X2,X3] = ExtractComponents(Xis);
    plot3(X1,X2,X3,'.','color',clrs(nt,:));
end
axis equal;axis off;
set(2,'paperunits','centimeters');
set(2,'papersize',[wd ht]);
set(2,'paperposition',[0 0 wd ht]);
print(2,'-dpng','-r600','figureSperm_c.png');

%% Solve swimming problem

nBeats = 5;
tRange=[0 2*pi*nBeats];
dt=2*pi*0.05;

%waveform
xyWaveFn=@DKActSpermWave;

% swimmer no. 1
args.phase=pi/5;
args.k=2*pi;
s = linspace(0,1,nns);
t = tRange(1):dt:tRange(2);
[S,T]=ndgrid(s,t);
x00{1}=[0;0;0.2];
B{1}=RotationMatrix(0*pi/3,3);
b10{1}=B{1}(:,1);
b20{1}=B{1}(:,2);
swimmer{1}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);
swimmer{1}.fn=@SpermModelGI;
% discretisation parameters - number of points
swimmer{1}.model.ns=40;
swimmer{1}.model.nh=4;
swimmer{1}.model.Ns=100;
swimmer{1}.model.Nh=10;
% semi-axes
swimmer{1}.model.a1=2.0/45;
swimmer{1}.model.a2=1.6/45;
swimmer{1}.model.a3=1.0/45;

% swimmer no. 2
% no args input sets phase = 0, k = 2*pi
x00{2}=[-1;0.5;0.2];
B{2}=RotationMatrix(pi/3,3);
b10{2}=B{2}(:,1);
b20{2}=B{2}(:,2);
swimmer{2}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn);
swimmer{2}.fn=@SpermModelGI;
% discretisation parameters - number of points
swimmer{2}.model.ns=40;
swimmer{2}.model.nh=4;
swimmer{2}.model.Ns=100;
swimmer{2}.model.Nh=10;
% semi-axes
swimmer{2}.model.a1=2.5/45;
swimmer{2}.model.a2=2.0/45;
swimmer{2}.model.a3=1.0/45;

% swimmer no. 3
args.phase=pi/8;
args.k=7*pi/3;
x00{3}=[0.5;-0.8;0.2];
B{3}=RotationMatrix(4*pi/3,3);
b10{3}=B{3}(:,1);
b20{3}=B{3}(:,2);
swimmer{3}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);
swimmer{3}.fn=@SpermModelGI;
% discretisation parameters - number of points
swimmer{3}.model.ns=40;
swimmer{3}.model.nh=4;
swimmer{3}.model.Ns=100;
swimmer{3}.model.Nh=10;
% semi-axes
swimmer{3}.model.a1=2.2/45;
swimmer{3}.model.a2=1.8/45;
swimmer{3}.model.a3=1.0/45;

% swimmer no. 4
args.phase=pi/8;
args.k=9*pi/3;
x00{4}=[1.0;1.1;0.2];
B{4}=RotationMatrix(2*pi/3,3);
b10{4}=B{4}(:,1);
b20{4}=B{4}(:,2);
swimmer{4}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);
swimmer{4}.fn=@SpermModelGI;
% discretisation parameters - number of points
swimmer{4}.model.ns=40;
swimmer{4}.model.nh=4;
swimmer{4}.model.Ns=100;
swimmer{4}.model.Nh=10;
% semi-axes
swimmer{4}.model.a1=2.2/45;
swimmer{4}.model.a2=1.8/45;
swimmer{4}.model.a3=1.0/45;

% swimmer no. 5
args.phase=pi/8;
args.k=8*pi/3;
s = linspace(0,1,nns);
t = tRange(1):2*pi*0.05:tRange(2);
[S,T]=ndgrid(s,t);
x00{5}=[-1.0;-0.6;0.2];
B{5}=RotationMatrix(2*pi/3,3);
b10{5}=B{5}(:,1);
b20{5}=B{5}(:,2);
swimmer{5}.model.F=ConstructInterpolantFromxyForm(S,T,xyWaveFn,args);
swimmer{5}.fn=@SpermModelGI;
% discretisation parameters - number of points
swimmer{5}.model.ns=40;
swimmer{5}.model.nh=4;
swimmer{5}.model.Ns=100;
swimmer{5}.model.Nh=10;
% semi-axes
swimmer{5}.model.a1=2.2/45;
swimmer{5}.model.a2=1.8/45;
swimmer{5}.model.a3=1.0/45;


% boundary
boundary.fn=@PlaneBoundary2;
boundary.model.h=0.4;
boundary.model.nx=16;
boundary.model.ny=15;
boundary.model.Nx=32;
boundary.model.Ny=30;
boundary.model.Lx=3;
boundary.model.Ly=3;
boundary.model.O=[0 0 0];

% numerical parameters
epsilon=0.25/45;
domain='i';
blockSize=0.2;

tic
fprintf('starting solver\n')

[t,z]=SolveMultiSwimmingTrajectoryAndForces(x00,b10,b20,tRange,swimmer,...
    boundary,epsilon,domain,blockSize);

solveTime = toc;
fprintf('CPU time taken = %f\n',solveTime)

save('multiResults.mat');

%% Panel e - trajectory over 5 beats

figure(3);clf;hold on;
Nsw=length(swimmer);
x0 = cell(1,Nsw);
for n=1:Nsw
    x0{n}=z(:,n:Nsw:n+2*Nsw);
    plot(x0{n}(:,1),x0{n}(:,2),'k');
end
hx=xlabel('\(x_1\) (flagellar lengths)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(3,'paperunits','centimeters');
set(3,'papersize',[wd 0.88*ht]);
set(3,'paperposition',[0 0.04*ht wd 0.04*ht+0.84*ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
axis equal;
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
box on;
set(gca,'tickdir','out');
print(3,'-dpdf','-r600','figureSperm_e.pdf');

%% Panel b, d, f - flow fields

% This can be generated independantly of solving by loading
% 'multiResults.mat'

Nx=30;
Ny=31;

xg=linspace(-boundary.model.Lx/2,boundary.model.Lx/2,Nx);
yg=linspace(-boundary.model.Ly/2,boundary.model.Ly/2,Ny);
[Xg,Yg]=ndgrid(xg,yg);Zg=0*Xg+boundary.model.h/2;



tt = linspace(t(1),t(end),630);
zz=interp1(t,z,tt,'pchip');

m=1;
figure(4);clf;set(4,'position',[1 29 1366 660]);
PlotMultiSwimmerVelocityField(Xg,Yg,Zg,tt,zz,m,swimmer,boundary,epsilon,domain,blockSize,2);hold on;
PlotMultiSpermSolution(4,swimmer,'q',tt,zz,m);hold off;
view([0 90]);xlim([-boundary.model.Lx/2,boundary.model.Lx/2]);ylim([-boundary.model.Ly/2,boundary.model.Ly/2]);axis equal;
xlim([-boundary.model.Lx/2,boundary.model.Lx/2]);ylim([-boundary.model.Ly/2,boundary.model.Ly/2]);
hx=xlabel('\(x_1\) (flagellar lengths)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(4,'paperunits','centimeters');
set(4,'papersize',[wd 0.88*ht]);
set(4,'paperposition',[0 0.04*ht wd 0.04*ht+0.84*ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
box on;
set(gca,'tickdir','out');
print(4,'-dpdf','-r600','figureSperm_b.pdf');

m=315;
figure(5);clf;set(5,'position',[1 29 1366 660]);
PlotMultiSwimmerVelocityField(Xg,Yg,Zg,tt,zz,m,swimmer,boundary,epsilon,domain,blockSize,2);hold on;
PlotMultiSpermSolution(5,swimmer,'q',tt,zz,m);hold off;
view([0 90]);xlim([-boundary.model.Lx/2,boundary.model.Lx/2]);ylim([-boundary.model.Ly/2,boundary.model.Ly/2]);axis equal;
xlim([-boundary.model.Lx/2,boundary.model.Lx/2]);ylim([-boundary.model.Ly/2,boundary.model.Ly/2]);
hx=xlabel('\(x_1\) (flagellar lengths)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(5,'paperunits','centimeters');
set(5,'papersize',[wd 0.88*ht]);
set(5,'paperposition',[0 0.04*ht wd 0.04*ht+0.84*ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
box on;
set(gca,'tickdir','out');
print(5,'-dpdf','-r600','figureSperm_d.pdf');

m=630;
figure(6);clf;set(6,'position',[1 29 1366 660]);
PlotMultiSwimmerVelocityField(Xg,Yg,Zg,tt,zz,m,swimmer,boundary,epsilon,domain,blockSize,2);hold on;
PlotMultiSpermSolution(6,swimmer,'q',tt,zz,m);hold off;
view([0 90]);xlim([-boundary.model.Lx/2,boundary.model.Lx/2]);ylim([-boundary.model.Ly/2,boundary.model.Ly/2]);axis equal;
xlim([-boundary.model.Lx/2,boundary.model.Lx/2]);ylim([-boundary.model.Ly/2,boundary.model.Ly/2]);
hx=xlabel('\(x_1\) (flagellar lengths)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(6,'paperunits','centimeters');
set(6,'papersize',[wd 0.88*ht]);
set(6,'paperposition',[0 0.04*ht wd 0.04*ht+0.84*ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
box on;
set(gca,'tickdir','out');
print(6,'-dpdf','-r600','figureSperm_f.pdf');

end