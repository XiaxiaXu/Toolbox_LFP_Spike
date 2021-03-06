function c=PDC_computation(u,nFreqs,metric,maxIP,alg,alpha,criterion,flag_mvarresidue)
% input:
%        u     - data in columns, for example [HP PFC]
%               note: to gain a good VAR model, I apply data
%               length=1s,fs=400Hz (Xiaxia)
%        nFreqs  - point numbers of the nyqist frequency(fs/2), usually 512 or 1024
%                  for exampl=512, every frquency point is (fs/2)/512:(fs/2)/512:(fs/2)
%        metric   euc  - Euclidean ==> original PDC
%                 diag - diagonal ==> gPDC (generalized )
%                 info - informational ==> iPDC
%        maxIP - externally defined maximum IP
%        alg   - for algorithm (1: Nutall-Strand),(2: mlsm)
%                              (3: Vieira Morf),  (4: ARfit)
%        criterion - for AR order selection =>
%                                   1: AIC; 2: Hanna-Quinn; 3: Schwartz;
%                                   4: FPE, 5: fixed order in MaxIP
%        alpha  -  PDC test significance level
%
% output:
%         c.pdc - original/generalized/informational PDC
%         c.th -  threshold level by Patnaik approximation
%         c.pdc_th - above threshold pdc values otherwise equal NaN, I used
%                    this value as PDC results (Xiaxia)
%         c.ic1,c.ic2 - superior and inferior confidence interval
%         c.p - VAR model order

% important output
%         c.c12 - the strength of channel 2 driving channel 1
%         c.c21 - the strength of channel 1 driving channel 2

[nSegLength,nChannels]=size(u);
if nSegLength < nChannels, error('The data might be transposed.'); end;

flgDetrend=1; % detrend the signal
flgStandardize=1;
%<***> Usually it's recommended to detrend the time series
gct_signif = alpha;  % Granger causality test significance. Choose other
% value if you have good reason for using different
% one from PDC statistical testing.
igct_signif = alpha; % Instantaneous Granger causality test significance level.
VARadequacy_signif = 0.05; % VAR model adequacy significance level
%==========================================================================

if flgDetrend,
    for i=1:nChannels, u(:,i)=detrend(u(:,i)); end;
end;

[nChannels,nSegLength]=size(u);
if nChannels > nSegLength, u=u.';
    [nChannels,nSegLength]=size(u);
end;

if flgStandardize,
    for i=1:nChannels, u(i,:)=u(i,:)/std(u(i,:)); end;
end;

%                            VAR model estimation
%==========================================================================
[IP,pf,A,~,~,ef,~,~,~] = mvar(u,maxIP,alg,criterion);

%    Testing for adequacy of MAR model fitting through Portmanteau test
%==========================================================================
h = 20; % testing lag
aValueVAR = 1 - VARadequacy_signif;
flgPrintResults = 1;

[~,Portmanteau,~,~]=mvarresidue(ef,nSegLength,IP,aValueVAR,h,...
    flgPrintResults);

flgPrintResults = 0;
[~, ~, Tr_igct, ~] = gct_alg(u,A,pf,igct_signif, ...
    flgPrintResults);
nPairsIGC = (sum(sum(Tr_igct==1)))/2;
if nPairsIGC == 0, % no instantaneous GCT test,next perform gPDC calculation
    %==========================================================================
    %            PDC, threshold and confidence interval calculation.
    %==========================================================================
    c=asymp_pdc(u,A,pf,nFreqs,metric,gct_signif);
    % Power spectra and coherence calculation
    %c.SS = ss_alg(A, pf, nFreqs);
    %c.coh = coh_alg(c.SS);
    %Statistically significant PDC on frequency scale
    if alpha ~= 0,
        pdc_temp = ((abs(c.pdc)-c.th) > 0).*c.pdc + ((abs(c.pdc)-c.th) <= 0)*(-1);
        pdc_temp(ind2sub(size(pdc_temp),find(pdc_temp == -1))) = NaN;
        c.pdc_th = pdc_temp;
        c.c12 = getCij(c.pdc_th,1,2,nFreqs);
        c.c21 = getCij(c.pdc_th,2,1,nFreqs);
    else
        
        c.c12 = getCij(c.pdc,1,2,nFreqs);
        c.c21 = getCij(c.pdc,2,1,nFreqs);
    end
    
else
    c=[];
end

if flag_mvarresidue
    if ~Portmanteau
        c=[];
    end
end

