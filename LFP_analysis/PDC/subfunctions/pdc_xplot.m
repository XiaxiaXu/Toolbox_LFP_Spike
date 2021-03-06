% Connectivity plot in matrix layout with power spectra in the main 
% diagonal.
%
%[hxlabel,hylabel] = pdc_xplot(c, ...
%                          flgPrinting,fs,w_max,chLabels,flgColors)
%
% input: c.{SS,Coh,L,Lpatnaik,LTra,L2vinf,L2vsup,metric}- data structure
%        flgPrinting= [1 1 1 1 1 0 1]; %An example: except coherence.
%         blue-line    | | | | | | 7 Spectra (0: w/o; 1: Linear; 2: Log)
%         gray         | | | | | 6 Coherence
%         dashed-blue  | | | | 5 Plot lower confidence limit
%         dashed-blue  | | | 4 Plot upper confidence limit
%         red          | | 3 Significant PDC in red line
%         dashed-black | 2 Patnaik threshold level in black dashed-line
%         green        1 PDC in green line or black w/o statistics
%
%        fs - sampling frequency
%        w_max - frequency scale upper-limit
%        chLabels - channel identification labels
%        flgColors - 0 or 1 for PDC 0.1 and 0.01 y-axis rescaling
%
% output: graphics
%         hxlabel, hylabel = graphic label's handles
%
% Examples: Calculate c using alg_pdc.
%
%           pdc_xplot(c); % Defaults flgPrinting, fs, w_max and flgColor
%
%           pdc_xplot(c,[1 1 1 0 0 1 1],200,100,[],0);
%                        % PDC, threshold, power spectra, coherence
%                        % plotted. fs=200 Hz; default channel label;
%                        % flgColor=0 => no color or rescalling is used.



function [hxlabel,hylabel] = pdc_xplot(c,...
                                flgPrinting,fs,w_max,chLabels,chLabelsY,flgColor)
knargin = 6;     % Number of input arguments
SS = c.SS;       % Spectra
Coh = c.coh;     % Coh^2
L = c.pdc;       % PDC^2 (N x N x freq)
Lpatnaik = c.th; % Patnaik threshold values for alpha
L2vinf = c.ic1; 
L2vsup = c.ic2;
LTra   = c.pdc_th; % Significant PDC^2 on freq range, otherwise = NaN.
metric = c.metric; % Valid optons: "euc", "diag" or "info"
[N,q,nFreqs]=size(L); % N=q = Number of channels/time series;
%                     % nFreqs = number of points on frequency scale
nodesett=1:N;

if nargin <  (knargin-4),  flgPrinting = [1 1 1 0 0 0 1]; end;
if nargin <= (knargin-4),  fs=1; end
if nargin < (knargin-2),   w_max = fs/2; end
if nargin < knargin,       flgColor = 0; end;

if w_max > fs/2 + eps,
   error(['The parameter w_max should be =< Nyquist frequency,' ...
      'i.e, w_max <= fs/2.'])
end;   

if c.alpha == 0,
   flgPrinting = flgPrinting .* [1 0 0 0 0 1 1];
end;

if exist('shadedplot.m','file') ~= 2, flgColor = 0; end;

w = 0:fs/(2*nFreqs):w_max-fs/(2*nFreqs);
nPlotPoints = length(w);
w_min = w(1);

if nargin < (knargin-1) || isempty(chLabels),
   if isfield(c,'chLabels'), 
      chLabels=c.chLabels; 
   else
      chLabels=[];
   end;
   if ~isempty(chLabels) && max(size(chLabels)) < N,
      error('1 NOT ENOUGH CHANNEL LABELS.');
   end;
elseif max(size(chLabels)) < N,
   if isfield(c,'chLabels'), 
      chLabels=c.chLabels; 
   else
      chLabels=[];
   end;
   if ~isempty(chLabels) && max(size(chLabels)) < N,
      error('2 NOT ENOUGH CHANNEL LABELS 2.');
   else
      disp('3 NOT ENOUGH CHANNEL LABELS. Default labels assumed.');
   end;
end;

hxlabel=0; % x-axis labels' handles
hylabel=0; % y-axis labels' handles
for j = 1:N
   s=nodesett(j);
   for i=1:N,
      r=nodesett(i);
      if j ~= i || ( j == i && flgPrinting(7) ~= 0)
         h=subplot2(N,N,(i-1)*N+j);
      end;
%==========================================================================
%                       Power spectrum plottinga
%==========================================================================
      if j == i, %Power spectrum
         y = abs(getCij(SS,r,s,nPlotPoints)); 
         pdc_tmp = abs(getCij(L,r,s,nPlotPoints));
         
         if flgPrinting(7), % Power spectrum linear scale
            if flgPrinting(7)==4,
               %Normalization of linear power scale
               y=y/max(y);
               [ax,h1,h2]=plotyy(h,w,y,w,pdc_tmp);
               if j == 2,
                  hh=ylabel('Spectra [a.u.]');
                  pos = get(hh,'Position');
                  pos(1) = 0.005;
                  set(hh,'Position',pos,'FontWeight','bold', ...
                     'FontSize',[8]);
               end;
            elseif flgPrinting(7)==5, % log Spectra
               y=log(y);
               y=(y-min(y))/(max(y)-min(y));
               [ax,h1,h2]=plotyy(h,w,y,w,pdc_tmp);
               if j == 2,
                  hh=ylabel('log Spectra [a.u.]');
                  pos = get(hh,'Position');
                  pos(1) = 0.005;
                  set(hh,'Position',pos,'FontWeight','bold', ...
                     'FontSize',[8]);
               end;
            elseif flgPrinting(7)==3,
               atrib = 'g-';
               plot(h,w,abs(getCij(L,r,s,nPlotPoints)),atrib,'LineWidth',2);
               ax(1)=gca;
               ax(2)=ax(1);
            elseif flgPrinting(7)==1,
               %Normalization of linear power scale
               y=y/max(y);
               h12 = plot(h,w,y);
               ax(1)=gca;
               ax(2)=ax(1);
               if j == 2,
                  hh=ylabel('Spectra [a.u.]');
                  pos = get(hh,'Position');
                  pos(1) = 0.005;
                  set(hh,'Position',pos,'FontWeight','bold', ...
                     'FontSize',[8]);
               end;
            elseif flgPrinting(7)==2,
               y=log(y);
               y=(y-min(y))/(max(y)-min(y));
               h12 = plot(h,w,y);
               ax(1)=gca;
               ax(2)=ax(1);
               if j == 2,
                  hh=ylabel('log Spectra [a.u.]');
                  pos = get(hh,'Position');
                  pos(1) = 0.005;
                  set(hh,'Position',pos,'FontWeight','bold', ...
                     'FontSize',[8]);
               end;
            end;
            hold on
            switch flgPrinting(7)
               case {1,2}
                  set(h12,'LineWidth',2,'Color',[0 0 1]);
               case 3
                  if flgPrinting(2)==1
                     atrib='k--'; % Patnaik signif. level in black line
                     plot(h,w,abs(getCij(Lpatnaik,r,s,nPlotPoints)), ...
                        atrib,'LineWidth',2);
                  end;
                  if flgPrinting(3)
                     atrib='r-'; % Significant PDCn in red line
                     plot(h,w,abs(getCij(LTra,r,s,nPlotPoints)), ... 
                        atrib,'LineWidth',2);
                     %hold on
                  end;
               case {4,5}
                  set(h1,'LineWidth',2,'Color',[0 0 1]);
                  set(h2,'LineWidth',2,'Color',[0 1 0]);
                  hold(ax(1), 'all'); hold(ax(2), 'all');
                  if flgPrinting(2)==1
                     atrib='k--'; % Patnaik signif. level in black line
                     plot(h,w,abs(getCij(Lpatnaik,r,s,nPlotPoints)), ...
                        atrib,'LineWidth',2);
                  end;
                  if flgPrinting(3)
                     atrib='r-'; % Significant PDCn in red line
                     plot(h,w,abs(getCij(LTra,r,s,nPlotPoints)), ...
                        atrib,'LineWidth',2.5);
                     %hold on
                  end;

               otherwise
                  disp('Unknown flgPrinting(7) option.')
            end;
            
            set(ax(1),'XLim',[w_min w_max],'XTick',[],'XTickLabel',[' '],...
               'YLim',[-0.02 1.05],'YTick',[]);
            if j == N && (flgPrinting(7) == 4 || flgPrinting(7) == 5)
               set(ax(2),'XLim',[w_min w_max],'XTick',[], ...
                  'XTickLabel',[' '],'YLim',[-0.02 1.05], ...
                  'YTick',[0 .5 1], 'YTickLabel',[' 0'; '.5'; ' 1']);
            else
               set(ax(2),'XLim',[w_min w_max],'XTick',[], ...
                  'XTickLabel',[' '],'YLim',[-0.02 1.05], ...
                  'YTick',[0 .5 1], 'YTickLabel',[' ']);
            end;

            set(h, 'Color', 0.8*[1 1 1]); % Background color

            if j == N,
               hxlabel(j)=labelitx(j,chLabels);
               set(h,'XTickLabel',[' ';' ';' ';' ';' ';' '])
            else
               set(h,'XTickLabel',[' ';' ';' ';' ';' ';' '])
            end;
         end;
         if j == 1,
            hylabel(i)=labelity(i,chLabelsY);
         end;
         if j == N,
            hxlabel(j)=labelitx(j,chLabels);
         end;
%==========================================================================
%                      PDC and coherence plotting
%==========================================================================
      else % PDC and coherence
         if flgPrinting(1),
            atrib='g-'; % PDCn in green lines
            Ltemp = abs(getCij(L,r,s,nPlotPoints));
            if flgColor
               if max(Ltemp) < 0.01,
                  shadedplot([w_min w_max],[0.01 0.01],[0.105 0.105], ...
                                                             [0.7 0.7 1]);
                  hold on
                  shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                                                           0.8*[1 0.7 1]);
               else
                  shadedplot([w_min w_max],[0.0 0.0],[0.1 0.1], ...
                                                             [0.7 0.7 1]);
                  hold on
                  shadedplot([w_min w_max],[0.0 0.0],[0.01 0.01], ...
                                                           0.8*[1 0.7 1]);
               end;
               hold on
            end;
            if c.alpha == 0, atrib='k-'; end; % Without statistics, PDC is
                                              % plotted in black lines.
            plot(h,w,abs(getCij(L,r,s,nPlotPoints)),atrib,'LineWidth',2);
            grid off

            if flgColor
               if max(Ltemp) < 0.01,
                  ylim=[-0.0002 0.0105];
               elseif max(Ltemp) < 0.1 && max(Ltemp) >= 0.01,
                  ylim=[-0.002 0.105];
               else
                  ylim=[-0.02 1.05];
               end;
            else
               ylim=[-0.02 1.05];
            end;
         end;
%======================================================================
%                         Labeling axis
%======================================================================
         if i == N,
            if j == 1,
               hylabel(i)=labelity(i,chLabelsY);
               hxlabel(j)=labelitx(j,chLabels);
               set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                  'XTick',[0 .1 .2 .3 .4 .5], ...
                 'XTickLabel',[' 0';'  ';'  ';'  ';'  ';'.5'], ...
                  'FontSize',[10],'FontWeight','bold', ...
                  'YTick',[0 .5 1],'YTickLabel',[' 0';'.5';' 1']);
%    'XTickLabel',[' 0';'  ';'  ';'  ';'  ';sprintf('%4.1f',w_max)], ...
            else
               hxlabel(j)=labelitx(j,chLabels);
               set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
                  'XTick',[0 .1 .2 .3 .4 .5], ...
                  'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
                  'YTick',[0 .5 1],'YTickLabel',['  ';'  ';'  ']);
            end;
         elseif i==1 && j == 2 && flgPrinting(7)==0,
            hylabel(i)=labelity(i,chLabelsY);
            set(h,'XLim', [w_min w_max],'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'YTick',[0 .5 1], 'YTickLabel',['  '; '  '; '  ']);
         elseif j == 1,
            hylabel(i)=labelity(i,chLabelsY);
            set(h,'XLim', [w_min w_max],'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'YTick',[0 .5 1], 'YTickLabel',['  '; '  '; '  '], ...
               'FontSize',[10],'FontWeight','bold');
            if i == N,
               set(h,'FontSize',[10],'FontWeight','bold', ...
                  'YTick',[0 .5 1], 'YTickLabel',[' 0'; '.5'; ' 1']);
            end;
         else
            set(h,'XLim', [w_min w_max], 'YLim',ylim, ...
               'XTick',[0 .1 .2 .3 .4 .5], ...
               'XTickLabel',[' ';' ';' ';' ';' ';' '], ...
               'FontSize',[10],'FontWeight','bold', ...
               'YTick',[0 .5 1], 'YTickLabel',['  '; '  '; '  ']);
         end;

         hold on

%==========================================================================
%==========================================================================

         if flgPrinting(2),
            atrib='k--'; % Patnaik significance level in black line
            plot(h,w,abs(getCij(Lpatnaik,r,s,nPlotPoints)),atrib, ...
                                                        'LineWidth',2);
         end;

% Plot lower and upper bound of significance level
      if flgPrinting(4) || flgPrinting(5),
         kSignif=abs(getCij(LTra,r,s,nPlotPoints));
         flgSignif = sum(~isnan(kSignif));
         indexSignif=~isnan(kSignif);
         indexNotSignif=isnan(kSignif);
         if flgPrinting(4) && flgSignif
            atrib='b--'; % Lower-bound significant interval in blue -- line
            temp= getCij(L2vinf,r,s,nPlotPoints); 
            % lower-bound can be negative.
            temp(indexNotSignif)=NaN;
            plot(h,w,temp, atrib,'LineWidth',2);
         end;
         if flgPrinting(5) && flgSignif
            atrib='b--'; % Upper-bound significant interval in blue -- line
            temp = abs(getCij(L2vsup,r,s,nPlotPoints));
            temp(indexNotSignif)=NaN;
            plot(h,w,temp, atrib,'LineWidth',2);
         end;
      end;
      if flgPrinting(6), %Coh plot in gray-line
         plot(w,abs(getCij(Coh,r,s,nPlotPoints)),'-','LineWidth',2, ...
                                               'Color',[.7 .7 .7]);
      end;
      if flgPrinting(3),
         atrib='r-'; % Significant PDCn in red line
         plot(h,w,abs(getCij(LTra,r,s,nPlotPoints)),atrib,'LineWidth',2);
         end;
      end; % PDC and coherence
   end;
end;

supAxes=[.08 .08 .84 .84];
[ax,h1]=suplabel('Frequency','x',supAxes);
set(h1, 'FontSize',[14], 'FontWeight', 'bold')

switch metric
   case 'euc'
      [ax,h1]=suplabel('{\mid\pi_{\it{{i}{j}}}{(\lambda)\mid}^{2}}', ...
                       'y',supAxes);
   case 'diag'
      [ax,h1]=suplabel('{\mid{_g}\pi_{\it{{i}{j}}}{(\lambda)\mid}^{2}}',...
         'y',supAxes);
   case 'info'
      [ax,h1]=suplabel('{\mid{_i}\pi_{\it{{i}{j}}}{(\lambda)\mid}^{2}}',...
         'y',supAxes);
   otherwise
      error('Unknown metric.')
end;
set(h1, 'FontSize',[14], 'FontWeight', 'bold')
pos=get(h1,'Position');
pos(1)=pos(1)+0.05; %0.0545
set(h1,'Position',pos) % Adjust ylabel position

% Adjusting axis labels positions.
for k=1:N,
   set(hxlabel(k),'Units','normalized');
   set(hylabel(k),'Units','normalized');
   pos=get(hylabel(k),'Position');
   pos(1)=-0.135*N/3;
   set(hylabel(k),'Position',pos);
   pos=get(hxlabel(k),'Position');
   pos(2)=-0.145*N/3;
   set(hxlabel(k),'Position',pos);
end;

%==========================================================================
function [hxlabel]=labelitx(j,chLabels) % Labels x-axis plottings
if isempty(chLabels)
   hxlabel=xlabel(['\bf{\it{j}} \rm{\bf{ = ' int2str(j) '}}']);
   %    hxlabel=xlabel(['$\bf{\it{j}} \rm{\bf{ = ' int2str(j) '}}$'],...
   %       'interpreter','LaTeX');
   set(hxlabel,'FontSize',[14],'FontWeight','bold', ...
      'FontName','Times')
else
   % hxlabel=xlabel(chLabels{j});

   hxlabel=xlabel(['$' chLabels{j} '$'],'interpreter','LaTeX');
   set(hxlabel,'FontSize',[12],'FontWeight','bold')
end;

%==========================================================================
function [hylabel]=labelity(i,chLabels) % Labels y-axis plottings
if isempty(chLabels)
   hylabel=ylabel(['\bf{\it{i}} \rm{\bf{ = ' int2str(i) '}}'],...
      'Rotation',90);
   set(hylabel,'FontSize',[14],'FontWeight','bold', ...
      'FontName','Times')
else
   % hylabel=ylabel(chLabels{i});

   hylabel=ylabel(['$' chLabels{i} '$'],'interpreter','LaTeX');
   set(hylabel,'FontSize',[12],'FontWeight','bold','Color',[0 0 0])
end;

%==========================================================================
