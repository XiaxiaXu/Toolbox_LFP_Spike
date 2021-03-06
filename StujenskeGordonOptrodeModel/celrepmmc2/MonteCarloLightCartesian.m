function [frac_abs,frac_trans,params,x,depth]=MonteCarloLightCartesian(modelnumber,r,NA,varargin)
%[frac_abs,frac_trans,params,x,depth]=montecarlolight_function(modelnumber,r,NA,wavelength)
%This function outputs a matrix of light intensity values for light
%emitted out of a fiber optic, simulated in cartesian coordinates.
%This model is based on the work of Jacques and colleagues (see Jacques,
%1998; Jacques, 2003; Jacques, 2010)
%Inputs:
%1. MODELNUMBER - 1 = Johansson, 2010 (default); 2 = Yaroslavsky et al.,
%2002; Data from one of these papers is used to linearly interpolate the
%scattering and absorption coefficient and the anisotropy parameter for the
%wavelength of choice. Otherwise, set modelnumber = 3 for custom input
%parameters, and specify optical parameters as the last input with the
%syntax: MonteCarloLight(3,r,NA,[a s g])
%where a, s, and g are the absorption coefficient, scattering coefficient,
%and anisotropy parameter, respectively
%2. R - radius of optical fiber in mm.
%3. NA - numerical aperture of fiber (between 0 and 1).
%4. WAVELENGTH - wavelength (nm) of light to use (model 1,2) -OR-
%   OPTICAL_PARAMETERS - [a s g] (model 3)
%
%Optional Inputs:
%[]=montecarlolight_function(modelnumber,...,nphotonpackets,zrange,rmax,
%   step)
%5. NPHOTONPACKETS - 100,000 * nphotonpackets are launched (requires
%integer value; rounds up if non-integer value is given). Packets are
%launched in multiples of 100,000 for computational efficiency. 100 is
%default.
%6. ZRANGE - [A B C], where A is the minimum depth value, B is the location
%of the fiber tip, and C is the maximum depth value. [-4 0 6] is default.
%Out of bounds photon packets are absorbed and removed from the simulation.
%7. RMAX - maximum x and y distances away from the fiber. 6 mm is default.
%8. STEP - step to use for discretization of space (mm). .01 is default.
%
%Outputs:
%1. FRAC_ABS - absorbed fraction, A [1/mm^3] (depth x x-disance x y-distance)
%2. FRAC_TRANS - fractional transport, T [1/mm^2] (depth x x-disance x y-distance)
%3. PARAMS - structure with the following fields: dstep, zrange, rmax,
%            absorption, scattering, g, r, NA
%4. X - xy-distances in mm corresponding to second and third dimensions of
%frac_abs and frac_trans
%5. DEPTH - depth in mm, corresponding to first dimension of these matrices
%
%Equations:
%T = FRAC_ABS/absorption coefficient
%Fluence Rate=Power * FRAC_TRANS [mW/mm^2]
%[Note: this script is slower than MonteCarloLight, so it should only be
%used if modifying to specify an assymmetry about the fiber]
%
%Written by Joseph M. Stujenske
%Stujenske et al., Cell Reports

tic %Initialize timing
%Error checking on inputs and sorting out the variable arguments:
if nargin<4
    error('You must specify wavelength, model type, fiber diameter and numerical aperture.')
end
if ~ismember(modelnumber,1:3)
    error('Model Number can only be 1 (Johansson, 2010), 2 (Yaroslavsky et al., 2002), or 3 (custom).')
end
if ~(NA>=0 && NA<=1)
    error('NA must be between 0 and 1')
end
if length(varargin)<1
    if modelnumber==3
        error('For custom parameters, you must specify them as the last output in the form [a s g]')
    else
        error('For interpolation, you must specify a wavelength')
    end
end
if length(varargin)>=2
    nphotonpackets=varargin{2};
    if nphotonpackets<1
        nphotonpackets=1;
        disp('Warning: Number of Photon Packets Increased to 1*100000 (minimum)')
    end
end
if length(varargin)>=3
    zrange=varargin{3};
end
if length(varargin)>=4
    rmax=varargin{4};
end
if length(varargin)>=5
    dstep=varargin{5};
end
if ~exist('nphotonpackets','var') || isempty(nphotonpackets)
    nphotonpackets=100;
end
if ~exist('zrange','var') || isempty(zrange)
    zrange=[-4 0 6];
end
if ~exist('rmax','var') || isempty(rmax)
    rmax=6;
end
if ~exist('dstep','var') || isempty(dstep)
    dr=.01; %discretization parameter for space in mm
else
    dr=dstep;
end

%Input their custom input parameters or interpolate from the models:
if modelnumber==3
    asgvar=varargin{1};
    absorption=asgvar(1);
    scattering=asgvar(2);
    g=asgvar(3);
else
    wavelength=varargin{1};
    if ((wavelength<480 || wavelength >900) && modelnumber==1) || ((wavelength<450 || wavelength>1064) && modelnumber==2)
        disp('Warning: Extrapolating Outside the Range of Experimentally Recorded Data')
    end
    %Data from Johansson, 2010 (Model 1) and Yaroslavsky, 2002 (Model 2)
    absorptionmat{1}=[480 0.37;...
        560	0.26;...
        580	0.19;...
        640	0.05;...
        780	0.02;...
        900	0.02];
    absorptionmat{2}=[450 0.07;...
        510 0.04;...
        630 0.02;...
        670 0.02;...
        1064 0.05];
    scatteringmat{1}=[480 11;...
        580 9.7;...
        640 9.0;...
        700 8.2;...
        780 7.8;...
        900 6.6];
    scatteringmat{2}=[450 11.7;...
        510 10.6;...
        630 9.0;...
        670 8.4;...
        1064 5.7];
    gmat{1}=[480 .89;...
        580 .89;...
        640 .89;...
        700 .90 ;...
        780 .90;...
        900 .90];
    gmat{2}=[450 .88;...
        510 .88;...
        630 .89;...
        670 .91;...
        1064 .9];
    
    %Interpolate from models
    absorption=interp1(absorptionmat{modelnumber}(:,1),...
        absorptionmat{modelnumber}(:,2),wavelength,'linear','extrap'); % absorption coefficient in mm^-1
    scattering=interp1(scatteringmat{modelnumber}(:,1),...
        scatteringmat{modelnumber}(:,2),wavelength,'linear','extrap'); %scattering coefficient in mm^-1
    g=interp1(gmat{modelnumber}(:,1),gmat{modelnumber}(:,2),wavelength,...
        'linear','extrap'); %anisotropy factor (between 0 and 1)
    g(g<0)=0;
    g(g>1)=1;
end

%Initialize variables
nphotonstotal_touse=100000*nphotonpackets; %Total number of photon packets to launch
nphotons=100000; %The number of photons to have "in play" at any particular time.
nlaunched=nphotons;
dz=dr; %discretization parameter for space in mm

n=1.36; %refractive index of the brain
th_out=asind(NA/n); %acceptance angle for the fiber
zmin=zrange(1); zcenter=zrange(2); zmax=zrange(3); %define bounds in z



catcher=zeros(length(zmin:dz:zmax),length(-rmax:dr:rmax),...
    length(-rmax:dr:rmax)); %matrix to accept energy from photons

newplotthreshold=nphotons; %How often to give command line progress update

%Initialize Progress Update
fprintf('\nProgress:\n')
stringout=[num2str(nphotons./...
    nphotonstotal_touse),'%% of photons launched\nEstimated Time Remaining: Estimating'];
fprintf(stringout)

%Launch first batch of photons
[cors,dirs]=photoninit(nphotons,th_out,r,zcenter);

w=ones(1,nphotons); %weight matrix which indicates the percentage of
%photon energy left in each packet

while 1
    %choose step size for photon packets
    s=-log(rand(1,nphotons))/(absorption+scattering);
    
    %move photon packets in the x,y, and z direction by the step size s
    cors=cors+repmat(s,3,1).*dirs;
    
    %out of bounds photon packets are killed:
    outofbounds=sqrt(sum(cors(1:2,:).^2))>=rmax | cors(3,:)>=zmax | cors(3,:)<zmin;
    w(outofbounds)=0;
    
    %put out of bounds photon packets "inbounds" to avoid errors later on
    cors(:,w==0)=[zeros(2,sum(w==0));ones(1,sum(w==0))*zcenter];
    rcors=sqrt(sum(cors(1:2,:).^2));
    
    %photon packet weight is stored in the catcher matrix:
    zs=floor((cors(3,:)-zmin)/dz)+1;
    xs=floor((cors(1,:)+rmax)/dr)+1;
    ys=floor((cors(2,:)+rmax)/dr)+1;
    ins=sub2ind(size(catcher),zs,xs,ys);
    wacumm = accumarray(ins',w);
    ns=find(wacumm~=0);
    catcher(ns)=catcher(ns)+wacumm(ns)*(absorption/(absorption+scattering));
    
    %reduce weight of the photon packets by what was absorbed:
    w=w*(scattering/(absorption+scattering));
    
    %change direction of the photon packets:
    if g>0 && g<=1
        costh=((1+g.^2-((1-g.^2)./(1-g+2*g.*rand(1,nphotons))).^2)./(2.*g));
    elseif g==0
        costh=2*rand(1,nphotons)-1;
    else
        error('g must be between 0 and 1');
    end
    phi=2*pi*rand(1,nphotons);
    sinth=sqrt(1-costh.^2);
    temp=sqrt(1-dirs(3,:).^2);
    
    uxyz=repmat(sinth,3,1).*[(dirs(1:2,:).*repmat(dirs(3,:).*cos(phi),2,1)+...
        [-dirs(2,:);dirs(1,:)].*repmat(sin(phi),2,1))./repmat(temp,2,1);-cos(phi).*temp]+...
        dirs.*repmat(costh,3,1);
    %Above line is equivalent to:
    %             uxx=sinth.*(xdir.*zdir.*cos(phi)-ydir.*sin(phi))./temp+xdir.*costh;
    %             uyy=sinth.*(ydir.*zdir.*cos(phi)+xdir.*sin(phi))./temp+ydir.*costh;
    %             uzz=-sinth.*cos(phi).*temp+zdir.*costh;
    
    %photon packets very close to being vertical are treated as being
    %vertical:
    tofix=abs(dirs(1,:))<.0001 & abs(dirs(2,:))<.0001;
    if any(tofix)
        uxyz(:,tofix)=[repmat(sinth(tofix),2,1).*[cos(phi(tofix));...
            sin(phi(tofix))];costh(tofix).*sign(dirs(3,tofix))];
        %Above line is equilavent to:
        %             uxx(tofix)=sinth(tofix).*cos(phi(tofix));
        %             uyy(tofix)=sinth(tofix).*sin(phi(tofix));
        %             uzz(tofix)=costh(tofix).*sign(zdir(tofix));
    end
    dirs=uxyz;
    
    %correct dirs to be a unit vector if it has deviated
    %due to rounding errors:
    mag=sqrt(sum(dirs.^2));
    dirs=dirs./repmat(mag,3,1);
    
    %roulette to see if packets with low weight die while preserving
    %conservation of energy:
    chance=rand(1,nphotons);
    w(w<1e-4 &chance<=.1& w>0)=w(w<1e-4 &chance<=.1& w>0)./.1;
    todestroy=(w<1e-4 & chance>=.1 & w>0) | outofbounds;
    ntodestroy=sum(todestroy);
    
    %replace destroyed photon packets with new photon packets
    if ntodestroy>0
        if ntodestroy+nlaunched<=nphotonstotal_touse
            [cors(:,todestroy),dirs(:,todestroy)]=photoninit(ntodestroy,th_out,r,zcenter);
            w(todestroy)=1;
            nlaunched=nlaunched+ntodestroy;
        elseif nlaunched<nphotonstotal_touse
            which=find(todestroy);
            replaceins=(which(1:nphotonstotal_touse-nlaunched));
            [cors(:,replaceins),dirs(:,replaceins)]=photoninit(length(replaceins),th_out,r,zcenter);
            w(replaceins)=1;
            nlaunched=nlaunched+length(replaceins);
            w(which(nphotonstotal_touse-nlaunched+1:end))=0;
        else
            w(w<1e-4 &chance>=.1 & w>0)=0;
        end
    end
    
    %terminate the simulation if all of the packets have no weight:
    if all(w==0)
        break
    end
    
    %Progress output
    if nlaunched>=max(newplotthreshold,nphotons+1)
        tout=toc;
        
        %Estimate Time to Destroy the Rest of the Photon Packets
        timeremaining=tout.*((nphotonstotal_touse/(nlaunched-nphotons))-1); 
        
        newplotthreshold=newplotthreshold+nphotons;
        removeout=repmat('\b',1,length(stringout)-2);
        stringout=[num2str(nlaunched./nphotonstotal_touse*100),...
            '%% of photons launched\nEstimated Time Remaining: ',...
            num2str(timeremaining),' s'];
        fprintf([removeout stringout])
    end
end

%Clean Up - Generate Outputs
totaltime=toc;
fprintf(['\nSimulation Time: ',num2str(totaltime),' s','\n'])

%Save data!
frac_abs=catcher./(dr.^2*dz)./(nphotonstotal_touse); %Correct for cylindrical geometry and make fraction
frac_trans=frac_abs/absorption;
params.dstep=dr;
params.zrange=zrange;
params.rmax=rmax;
params.absorption=absorption;
params.scattering=scattering;
params.g=g;
params.r=r;
params.NA=NA;
x=-rmax:dr:rmax;
y=-rmax:dr:rmax;
depth=zmin:dz:zmax;
end

function [cors,dirs]=photoninit(num,th_out,r,zcenter)
%directions are encoded in spherical coordinates:
thinit=cosd(th_out)+rand(1,num)*(1-cosd(th_out)); %initial angle in the x-y plane
phiinit=-180+2*rand(1,num)*180; %inital angle in the z plane
%initialize x,y, and z direction matrices such that [x y z] is a unit
%vector
xdir=sqrt(1-thinit.^2).*cosd(phiinit);
ydir=sqrt(1-thinit.^2).*sind(phiinit);
zdir=thinit;

%choosing photon starting positions on the fiber as a uniform distribution
%covering the area of the circular aperture of the fiber
u=rand(1,num)+rand(1,num);
rinit=zeros(1,num);
rinit(u>1)=2-u(u>1);
rinit(u<=1)=u(u<=1);
rinit=rinit*r;
thinitpos=-180+2*rand(1,num)*180;

%initializing x, y, and z coordinates on the surface of the fiber
xcor=rinit.*cosd(thinitpos);
ycor=rinit.*sind(thinitpos);
zcor=ones(1,num)*zcenter;
cors=[xcor;ycor;zcor];
dirs=[xdir;ydir;zdir];
end