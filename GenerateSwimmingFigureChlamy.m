%GENERATESWIMMINGFIGURECHLAMY Generates figure 1 from the manuscript
%
function GenerateSwimmingFigureChlamy

% Font options for printing figures
fs=8;
fn='times';
wd=7.0;
ht=6.0;

%% Panel a and c  - image of chlamy beat pattern

% Generate swimmer ------------------------------------------------------------------------
nbeats=1;
tRange=[0 2*pi*nbeats];

% generate gridded interpolant ------------------------
s = linspace(0,1,30);
t = (tRange(1):2*pi*0.05:tRange(2));
[S,T]=ndgrid(s,t);
stWaveFn=@ChlamyFromModel;
swimmer.model.F=ConstructInterpolantFromSTForm(S,T,stWaveFn,[]);
%------------------------------------------------------

swimmer.fn=@ChlamyModel1;
% discretisation parameters - number of points
swimmer.model.ns=40;
swimmer.model.nh=4;
swimmer.model.Ns=400;
swimmer.model.Nh=10;
% chlamy body semi-axes
swimmer.model.a1=0.5;
swimmer.model.a2=0.6;
swimmer.model.a3=0.6;
% chlamy flagellar angle
swimmer.model.ang=pi/5;
% -----------------------------------------------------------------------------------------

% Panel a ----------------------------------
figure(1);clf;hold on;
nth=40;
nphi=80;
a=1;
[x3,x2,x1]=GenerateSphereSurfaceForVisualisation(nth,nphi,a);
x1=x1*swimmer.model.a1;
x2=x2*swimmer.model.a2;
x3=x3*swimmer.model.a3;
surf(x1,x2,x3,0*x3);shading flat;light;
clrs = summer(length(t));
for nt=1:length(t)
    [xis,~,~]  = swimmer.fn(t(nt),swimmer.model);
    [x1,x2,x3] = ExtractComponents(xis);
    plot3(x1,x2,x3,'.','color',clrs(nt,:));
end
plot3(x1(1:swimmer.model.nh.^2*6),x2(1:swimmer.model.nh.^2*6), ...
    x3(1:swimmer.model.nh.^2*6),'.','color',clrs(round(3*nt/4),:))
axis equal;axis off;
set(1,'paperunits','centimeters');
set(1,'papersize',[wd ht]);
set(1,'paperposition',[0 0 wd ht]);
set(gca,'fontname',fn);set(gca,'fontsize',fs);
print(1,'-dpng','-r600','figureChlamy_a.png');

% Panel c ----------------------------------
figure(2);clf;hold on;

[x3,x2,x1]=GenerateSphereSurfaceForVisualisation(nth,nphi,a);
x1=x1*swimmer.model.a1;
x2=x2*swimmer.model.a2;
x3=x3*swimmer.model.a3;

surf(x1,x2,x3,0*x3);shading flat;light;
clrs = spring(2*length(t));
clrs = clrs(round(length(t)/2)+1:end,:);
for nt=1:length(t)
    [~,~,Xis]  = swimmer.fn(t(nt),swimmer.model);
    [X1,X2,X3] = ExtractComponents(Xis);
    plot3(X1,X2,X3,'.','color',clrs(nt,:));
end
plot3(x1(1:swimmer.model.nh.^2*6),x2(1:swimmer.model.nh.^2*6), ...
    x3(1:swimmer.model.nh.^2*6),'.','color',clrs(round(3*nt/4),:))
axis equal;axis off;
set(2,'paperunits','centimeters');
set(2,'papersize',[wd ht]);
set(2,'paperposition',[0 0 wd ht]);
set(gca,'fontname',fn);set(gca,'fontsize',fs);
print(2,'-dpng','-r600','figureChlamy_c.png');


%% Solve swimming problem

nbeats=5;
tRange=[0 2*pi*nbeats];
x00=[0;0;0];

dt=2*pi*0.05;
B=RotationMatrix(0*pi/3,3);
b10=B(:,1);
b20=B(:,2);

% generate gridded interpolant ------------------------
s = linspace(0,1,30);
t = (tRange(1):2*pi*0.05:tRange(2));
[S,T]=ndgrid(s,t);
stWaveFn=@ChlamyFromModel;
swimmer.model.F=ConstructInterpolantFromSTForm(S,T,stWaveFn,[]);
%------------------------------------------------------

epsilon=0.001;
domain='i';
blockSize=0.2;

boundary=[];

tic
fprintf('starting solver\n')

[t,z]=SolveSwimmingTrajectoryAndForces(x00,b10,b20,tRange,swimmer,boundary,epsilon,domain,blockSize);
solveTime=toc;

fprintf('CPU time taken = %f\n',solveTime)

save('biflagellateResults.mat')

%% Panel e - trajectory over 5 beats

figure(3);clf;
plot(t/2/pi,z(:,2),'k');
hx=xlabel('\(t\) (beats)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(3,'paperunits','centimeters');
set(3,'papersize',[wd 0.88*ht]);
set(3,'paperposition',[0 0.04*ht wd 0.04*ht+0.84*ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
xlim([0 nbeats]);
box on;
set(gca,'tickdir','out');
print(3,'-dpdf','-r600','figureChlamy_e.pdf');

%% Panel b, d, f - flow fields

% This can be generated independantly of solving by loading
% 'biflagellateResults.mat'

Nx=15;
Ny=16;
xg=linspace(-1.5,1.5,Nx);
yg=linspace(-1.0,2.0,Ny);
[Xg,Yg]=ndgrid(xg,yg);Zg=0*Xg;

tt=(0:dt:2*pi)';
zz=interp1(t,z,tt,'pchip');
scl=2;
mksz=0.75;

fig=4;figure(fig);clf;hold on;m=floor(length(tt)/3);

PlotMultiSwimmerVelocityField(Xg,Yg,Zg,tt,zz,m,swimmer,boundary,epsilon,domain,blockSize,scl);
PlotMultiSpermSolution(fig,swimmer,'q',tt,zz,m,mksz);
axis equal;
hx=xlabel('\(x_1\) (flagellar lengths)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(fig,'paperunits','centimeters');
set(fig,'papersize',[wd ht]);
set(fig,'paperposition',[0 0 wd ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
box on;
set(gca,'tickdir','out');
print(fig,'-dpdf','-r600','figureChlamy_b.pdf');

fig=5;figure(fig);clf;hold on;m=floor(2*length(tt)/3);

PlotMultiSwimmerVelocityField(Xg,Yg,Zg,tt,zz,m,swimmer,boundary,epsilon,domain,blockSize,scl);
PlotMultiSpermSolution(fig,swimmer,'q',tt,zz,m,mksz);
axis equal;
hx=xlabel('\(x_1\) (flagellar lengths)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(fig,'paperunits','centimeters');
set(fig,'papersize',[wd ht]);
set(fig,'paperposition',[0 0 wd ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
box on;
set(gca,'tickdir','out');
print(fig,'-dpdf','-r600','figureChlamy_d.pdf');

fig=6;figure(fig);clf;hold on;m=floor(3*length(tt)/3);

PlotMultiSwimmerVelocityField(Xg,Yg,Zg,tt,zz,m,swimmer,boundary,epsilon,domain,blockSize,scl);
PlotMultiSpermSolution(fig,swimmer,'q',tt,zz,m,mksz);
axis equal;
hx=xlabel('\(x_1\) (flagellar lengths)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(fig,'paperunits','centimeters');
set(fig,'papersize',[wd ht]);
set(fig,'paperposition',[0 0 wd ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
box on;
set(gca,'tickdir','out');
print(fig,'-dpdf','-r600','figureChlamy_f.pdf');

end