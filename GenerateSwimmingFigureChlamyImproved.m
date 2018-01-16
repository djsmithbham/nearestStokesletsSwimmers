% CreateSwimmingPaperChlamyFigure.m

clear all;

fs=8;
fn='times';
wd=7.0;
ht=6.0;

% Panels a/c image of chlamy beat pattern

nbeats=1;
tRange=[0 2*pi*nbeats];
x00=[0;0;0];

dt=2*pi*0.05;
B=RotationMatrix(0*pi/3,3);
b10=B(:,1);
b20=B(:,2);

% generate gridded interpolant ------------------------
s = linspace(0,1,30);
t = [tRange(1):2*pi*0.05:tRange(2)];
[S,T]=ndgrid(s,t);
fourierParams=load('paramecium_fourier_coefs.mat');
stWaveFn=@CiliumFromFourier;
swimmer.model.F=ConstructInterpolantFromSTForm(S,T,stWaveFn,fourierParams);
%------------------------------------------------------

swimmer.fn=@ChlamyModel1;
% discretisation parameters - number of points
swimmer.model.ns=40;
swimmer.model.nh=4;
swimmer.model.Ns=100;
swimmer.model.Nh=10;
% chlamy body semi-axes
swimmer.model.a1=0.5;
swimmer.model.a2=0.6;
swimmer.model.a3=0.6;
% chlamy flagellar angle
swimmer.model.ang=pi/6;

figure(1);clf;hold on;
nth=40;
nphi=80;
a=1;
[x3,x2,x1]=GenerateSphereSurfaceForVisualisation(nth,nphi,a);
x1=x1*swimmer.model.a1;
x2=x2*swimmer.model.a2;
x3=x3*swimmer.model.a3;
surf(x1,x2,x3,0*x3);shading flat;light;
for nt=1:length(t)
    [xis,~,~]  = swimmer.fn(t(nt),swimmer.model);
    [x1,x2,x3] = ExtractComponents(xis);
    plot3(x1,x2,x3,'c.');
end
axis equal;axis off;
set(1,'paperunits','centimeters');
set(1,'papersize',[wd ht]);
set(1,'paperposition',[0 0 wd ht]);
%set(gca,'fontname',fn);set(gca,'fontsize',fs);
print(1,'-dpng','-r600','figureChlamy_a.png');

clf;hold on;
[x3,x2,x1]=GenerateSphereSurfaceForVisualisation(nth,nphi,a);
x1=x1*swimmer.model.a1;
x2=x2*swimmer.model.a2;
x3=x3*swimmer.model.a3;
surf(x1,x2,x3,0*x3);shading flat;light;
for nt=1:length(t)
    [~,~,Xis]  = swimmer.fn(t(nt),swimmer.model);
    [X1,X2,X3] = ExtractComponents(Xis);
    plot3(X1,X2,X3,'m.');
end
axis equal;axis off;
set(1,'paperunits','centimeters');
set(1,'papersize',[wd ht]);
set(1,'paperposition',[0 0 wd ht]);
%set(gca,'fontname',fn);set(gca,'fontsize',fs);
print(1,'-dpng','-r600','figureChlamy_c.png');


%%
% panel e - trajectory over 5 beats

nbeats=1;
tRange=[0 2*pi*nbeats];
x00=[0;0;0];

dt=2*pi*0.05;
B=RotationMatrix(0*pi/3,3);
b10=B(:,1);
b20=B(:,2);

% generate gridded interpolant ------------------------
s = linspace(0,1,30);
t = [tRange(1):2*pi*0.05:tRange(2)];
[S,T]=ndgrid(s,t);
fourierParams=load('paramecium_fourier_coefs.mat');
stWaveFn=@CiliumFromFourier;
swimmer.model.F=ConstructInterpolantFromSTForm(S,T,stWaveFn,fourierParams);
%------------------------------------------------------

%%

t=[tRange(1):dt:tRange(2)];Nt=length(t);

epsilon=0.001;
domain='i';
rho=0.5;
blockSize=0.2;

[xis,~,~] = swimmer.fn(0,swimmer.model);DOF=length(xis);z0=[x00;b10;b20;zeros(DOF,1)];boundary=[];
%z=zeros(Nt,length(z0));

%%

tic

% improved code matches old code very well
% for the chlamy case, the old code is faster (as it's a smaller problem
% probably)
%[t,z]=SolveSwimmingTrajectoryAndForces(        x00,b10,b20,tRange,swimmer,boundary,epsilon,domain,    blockSize);
[t,z]=SolveSwimmingTrajectoryAndForcesImproved(x00,b10,b20,tRange,swimmer,boundary,epsilon,domain,rho,blockSize);
toc

%%
figure(1);clf;
plot(t/2/pi,z(:,2),'k');
hx=xlabel('\(t\) (beats)','interpreter','latex');
hy=ylabel('\(x_2\) (flagellar lengths)','interpreter','latex');
set(1,'paperunits','centimeters');
set(1,'papersize',[wd 0.88*ht]);
set(1,'paperposition',[0 0.04*ht wd 0.04*ht+0.84*ht]);
set(gca,'fontsize',fs); set(gca,'fontname',fn);
set(hx,'fontsize',fs);  set(hx,'fontname',fn);
set(hy,'fontsize',fs);  set(hy,'fontname',fn);
xlim([0 nbeats]);
box on;
set(gca,'tickdir','out');
print(1,'-dpdf','-r600','figureChlamyode113eps0.001I_e.pdf');

%%

Nx=15;
Ny=16;
xg=linspace(-1.5,1.5,Nx);
yg=linspace(-1.0,2.0,Ny);
[Xg,Yg]=ndgrid(xg,yg);Zg=0*Xg;

tt=[0:dt:2*pi]';
zz=interp1(t,z,tt,'pchip');
scl=2;
mksz=0.75;

fig=1;figure(fig);clf;hold on;m=floor(length(tt)/3);

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

fig=2;figure(fig);clf;hold on;m=floor(2*length(tt)/3);

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

fig=3;figure(fig);clf;hold on;m=floor(3*length(tt)/3);

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

