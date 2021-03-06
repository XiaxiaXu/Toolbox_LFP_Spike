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

close all
clear all
savepath='D:\Projects\Project-OptogeneticHP-IUE\Final results\Fiber light\50um\';
savepath='D:\Projects\Project-OptogeneticHP-IUE\Final results\Fiber light\105um\';

%wavelength=594; % nm yellow laser
wavelength=473; % nm
NA=.22; % numerical aperture

Radius=.105/2; %[.050/2 .105/2] in mm REMEMBER TO CHANGE SAVEPATH

%light
[frac_abs,frac_trans,params,r,depth]=MonteCarloLight(1,Radius,NA,wavelength,100,[-2 0 4],3,.01);

PowerList=[7]; % mW
for m2=1:length(PowerList)
Power=PowerList(m2)
frac_trans_in_mWpermm2=Power*[frac_trans(:,end:-1:2) frac_trans]; %mW/mm^2 % Fluence rate
            
end          
            

h=figure; 
LightHeatPlotter_Modified(r,depth,log10(frac_trans./max(frac_trans(:))),wavelength);caxis([-4 0]);c=colorbar;ylabel(c,'LightPower');title(strcat('WL',num2str(wavelength),' NA',num2str(NA),' r',num2str(Radius),' P',num2str(Power),' F',num2str(Freq)))
hold on;contour([-r(end:-1:2) r],depth,frac_trans_in_mWpermm2,[5 10]);hold off
title(strcat(num2str(wavelength),'nm,', num2str(Power),'mW'))



h=figure; imagesc([-r(end:-1:2) r],depth,[frac_trans(:,end:-1:2) frac_trans]); colormap('hot');caxis([0 4.2]); xlabel('t (s)'); ylabel('Depth (mm)');c=colorbar;ylabel(c,'TempChange°C');title(strcat('WL',num2str(wavelength),' NA',num2str(NA),' r',num2str(Radius),' P',num2str(Power),' F',num2str(Freq)))
hold on;contour([-r(end:-1:2) r],depth,frac_trans_in_mWpermm2,[10 1]);hold off
saveas(h,strcat(savepath,'frac_trans','WL',num2str(wavelength),'-NA',num2str(NA*100),'-r',num2str(Radius*2*1000),'-P',num2str(Power),'.fig'))
saveas(h,strcat(savepath,'frac_trans','WL',num2str(wavelength),'-NA',num2str(NA*100),'-r',num2str(Radius*2*1000),'-P',num2str(Power),'-F',num2str(Freq),'.epsc'))
close all

        
