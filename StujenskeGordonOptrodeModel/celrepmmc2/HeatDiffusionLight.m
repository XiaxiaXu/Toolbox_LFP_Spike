function [u_time,u_space,t,r,depth]=HeatDiffusionLight(frac_abs,params,u_step,t_on,t_off,t_max,Power,t_save,r_avg)
%[u_time,t,u_space]=HeatDiffusionLight(frac_abs,params,u_step,t_on,t_off,t_max,Power,t_save,r_avg)
%This function models temperature change (in degrees Celsius) due to light
%from an optical fiber in a block of brain tissue based on the Pennes
%Bio-heat Equation (Pennes, 1948). Must be run after MonteCarloLight.
%
%Inputs:
%1. FRAC_ABS - Fraction absorbed (output of MonteCarloLight)
%2. PARAMS - Model parameters (output of MonteCarloLight)
%3. U_STEP - Discretization step in mm for model (does not need to match
%light model, and in fact, it is better for the light model to have a
%smaller discretization factor. This discretization step should be larger
%because the heat equation solution is more computationally taxing. Default
%is .03 mm.
%4. T_ON - Times when the laser turns on (e.g. [0, 1, 2]) in s. Default is 0.
%5. T_OFF - Times when the laser turns off (e.g. [0.05 1.05 2.05]) in s.
%Default is t_max.
%6. T_MAX - Time to end simulation. Default is 120 or the largest value of
%t_on/t_off (whichever is larger).
%7. POWER - Laser power in mW. Default is 1 mW.
%8. T_SAVE - Times to save heat in u_space in s. Default is t_max.
%9. R_AVG - Radial distance in mm over which to average for u_time. Default
%is .25 mm.
%
%Outputs:
%1. U_TIME - matrix with temperature change (degrees Celsius) over time (s)
%at each depth (depth x time) relative to baseline. 0 indicates no change.
%2. U_SPACE - matrix with temperature values in the plane of the fiber
% [depth x radial distance x length(t_save)]
%3. T - times (s) corresponding to second dimension of u_time.
%4. R - radial distance (mm) corresponding to second dimension of u_space.
%5. DEPTH - depths (mm) corresponding to first dimension of u_time and
%u_space.
%
%Written by Joseph M. Stujenske
%Stujenske et al., Cell Reports

if nargin<2
    error('You must specify absorbed fraction and params.');
end
if ~exist('u_step','var') || isempty(u_step)
    u_step=.03; %discretization step in mm
end
if ~exist('Power','var') || isempty(Power)
    Power=1; %output power in mW
end
if ~exist('t_on','var') || isempty(t_on)
    t_on=0; %times when the laser is on
end
if ~exist('t_max','var') || isempty(t_max)
    t_max=max([120 max(t_on) max(t_off)]); %time to end simulation
end
if ~exist('t_off','var') || isempty(t_off)
    t_off=t_max; %time when the laser is off; default is on and then off
end
if ~exist('t_save','var') || isempty(t_save)
    t_save=t_max; %output power in mW
end
if ~exist('r_avg','var') || isempty(r_avg)
    r_avg=.25; %output power in mW
end
if r_avg<u_step/2
    r_avg=u_step/2;
    disp('Warning: Changing r_avg value to include central two voxels (minimum r)')
end

%Initialize Variables
dr=params.dstep; %discretization step for the light model
zmin=params.zrange(1);
zmax=params.zrange(3);
rmax=params.rmax; %mm; how far in the x-y plane to calculate heat
corrector=rem(rmax,u_step);
rmax=rmax+corrector;

%Constants from Literature (Elwassif et al., 2006; Janssen et al., 2005)
density=1.04*10^-6; %density of brain
spheat=3650*10^3; %specific heat of brain
kbrain=.527; %whole brain thermal diffusivity constant;
wb=8.5*10^-3; %blood perfusion rate of brain (units converted)
pblood=3.6*10^6; %specific heat of blood (units converted)
densblood=1.06*10^-6; %density of blood (units converted)
%(units converted)
Tinit=37; %Temperature of brain at rest
Ta=36.7; %Arterial temperature
qm=(Tinit-Ta)*wb*pblood*densblood; %value chosen to balance equation

%Time Step
deltat=(u_step^2)/(6*(kbrain/(spheat*density))); %s; this is the
%minimum value needed for numerical stability

%Initialize vectors and matrices
depth=(zmin:u_step:zmax)';
r=0:u_step:rmax;
I=zeros(length(r),length(depth));

%CALCULATING LIGHT INTENSITY FOR EVERY X AND Z:
lightmodelrs=0:dr:rmax;
for rep=1:length(r)
    [~,in]=min(abs(r(rep)-lightmodelrs));
    I(rep,:)=frac_abs(1:u_step/dr:end,in)*Power;
end

%Initialize constants for later calculations:
uchange=(I/(density*spheat))*deltat; %change in temperature due to light
u=ones(length(r),length(depth))*Tinit; %temperature at t = 0
pc=1/(spheat*density); %constant for later calculations
k=ones(length(r),length(depth))*kbrain; %This matrix can be made non-homogenous
%for your applications, but numerical stability will not be assured. For
%this reason, this functionality is not built-in.

%Initalize Variable to Save Data
u_time=zeros(length(depth),length(0:deltat:t_max));
u_time(1,:)=Tinit;
t=0:deltat:t_max;
stepper=1;
ri=repmat(r',1,length(depth))+u_step/2;

%Initialize Progress Output
u_space=zeros(length(depth),length(r),length(t_save));
fprintf('\nProgress:\n')
stringout=['t=',num2str(0),' out of ',num2str(t_max),' s'];
fprintf(stringout)
tplot=0;
tsavecount=1;

%%RUN SIMULTATION%%
for t_loop=deltat:deltat:t_max;
    stepper=stepper+1;
    [m,n] = size(u);
    
    %Finite Difference Method for Heat Equation
    vrr = ((k(1:m,1:n).*(ri+u_step/2)).*u([2:m m],1:n) -(k(1:m,1:n).*(ri+u_step/2)+k([1 1:m-1],1:n).*(ri-u_step/2)).*u(1:m,1:n) +...
        (k([1 1:m-1],1:n).*(ri-u_step/2)).*u([1 1:m-1],1:n))./(u_step.^2*ri);
    vzz = ((k(1:m,1:n)).*u(1:m,[2:n n]) -(k(1:m,1:n)+k(1:m,[1 1:n-1])).*u(1:m,1:n) +...
        (k(1:m,[1 1:n-1])).*u(1:m,[1 1:n-1]))./(u_step.^2);
    
    %Discretization of Laplacian: Deltav = v_xx + v_yy + v_zz
    deltau = (vrr+vzz)*deltat*pc;
    
    %Heat Change Due to Perfusion by Blood and Metabolic Heat
    deltaperfusion=((Ta-u)*pblood*densblood*wb+qm)*pc*deltat;
    
    %Check if light is on or off
    %Total Heat Change
    
    t_ondiff=t_on-t_loop;
    t_offdiff=t_off-t_loop;
    if (isempty(t_offdiff(t_offdiff<=0)) && ~isempty(max(t_ondiff(t_ondiff<=0)))) || (~isempty(max(t_ondiff(t_ondiff<=0))) && max(t_ondiff(t_ondiff<=0))>max(t_offdiff(t_offdiff<=0)))
        uchange1=(deltau+uchange+deltaperfusion);
    else
        uchange1=(deltau+deltaperfusion);
    end
    
    
    %Absorbing ends
    uchange1(end,:)=0;
    uchange1(:,1)=0;
    uchange1(:,end)=0;
    
    %Increase Temperature
    u=u+uchange1;
    
    %Save Data
    u_time(:,stepper)=(nansum(u((r+u_step/2)<=r_avg,:).*repmat(r(r+u_step/2<=r_avg)',1,size(u,2)))./(nansum(r(r+u_step/2<=r_avg),2)))-Tinit;
    if tsavecount<=length(t_save) && t_save(tsavecount)-t_loop<=deltat && t_save(tsavecount)-t_loop>0
        u_space(:,:,tsavecount)=u'-Tinit;
        tsavecount=tsavecount+1;
    end
    
    %Progress Update
    if t_loop>tplot
        removeout=repmat('\b',1,length(stringout));
        stringout=['t=',num2str(floor(t_loop)),' out of ',num2str(t_max),' s'];
        fprintf([removeout stringout])
        tplot=tplot+1; %Update Every Second
    end
    
end
%%END SIMULATION%%