%Example usage of Light and Heat functions
%Generate 3 plots for 532 nm light:
%Figure 1 - Intensity relative to maximum intensity (log scale) from a 31
%nm, .22 NA fiber
%Figure 2 - Steady state temperature change for the same fiber
%Figure 3 - Time course of temperature change for the same fiber during and
%after 60 seconds of illumination
%
%Written by Joseph M. Stujenske
%Stujenske et al., Cell Reports
%
%Anticipate this to take >30 minutes to compute.
NA=.22;
wavelength=473; % nm
Radius=.105/2; % mm
Power=2; % mW
LightOn =3;%[.000 .125 .250 .375 .500 .625 .750 .875 1.000 1.125 1.250 1.375 1.500 1.625 1.750 1.875 2.000 2.125 2.250 2.375 2.500 2.625 2.750 2.875]; % s
LightOff=6;%[.003 .128 .253 .378 .503 .628 .753 .878 1.003 1.128 1.253 1.378 1.503 1.628 1.753 1.878 2.003 2.128 2.253 2.378 2.503 2.628 2.753 2.878]; % s
PlotDuration=9; % s


[frac_abs,frac_trans,params,r,depth]=MonteCarloLight(1,Radius,NA,wavelength,100);
figure; LightHeatPlotter(r,depth,log10(frac_trans./max(frac_trans(:))),wavelength);caxis([-4 0])

[u_time,u_space,t,r,depth]=HeatDiffusionLight(frac_abs,params,.03,LightOn,LightOff,PlotDuration,Power,6,.25);
figure; LightHeatPlotter(r,depth,u_space,'hot');caxis([0 4.2])
figure; imagesc(t,depth,u_time); colormap('hot');caxis([0 2.2]); xlabel('t (s)'); ylabel('Depth (mm)')
