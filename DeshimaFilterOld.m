close all;
clear all;
clc;
  

%% CONSTANTS
mu_0   = pi*4e-7;             % []    Permeability of FS
eps_0  = 8.854187817620e-12;  % []    Permittivity of FS
c_0    = 1/sqrt(eps_0*mu_0);  % [m/s] Speed of light in FS


%% Plotting parameters
mindB = -25;

plotSanityChecks        = 0;
plotAllFilters          = 0;
plotEnvelope            = 1;
plotLthru               = 1;
plotQfactors            = 1;
plotOverlappingPeaks    = 0;


%% OBSERVATION BW
f_min = 200e9; % lowest frequency
f_max = 460e9; % highest frequency
nF    = 4e3; %number of frequency samples
f     = linspace(f_min,f_max,nF); % frequency
w     = 2*pi*f; % angular freq


%% FILTER-BANK PROPERTIES
Qi=Inf;%3300;Inf; % It can be a vector
sep=1/4;%0.133;%0.1625;%1/4;%1/2;%3/4;
Ql_target = 100;
% Qc=10000;
% Ql_expected=1./(1./Qi+2./Qc);
lossesOnlyInRes = 0;
S = 1; % Oversampling

f0_lowest  = 220e9;
f0_highest = 2*f0_lowest;
nFilters = floor( 1+log10(f0_highest/f0_lowest)/log10(1+S/Ql_target) );
f0(1)    = f0_lowest;
for filter_i=2:nFilters
   f0(filter_i) = f0(filter_i-1) + f0(filter_i-1)*S/Ql_target;
end

f0 = sort(f0,'descend');


%% THROUGH-LINE
Z0_thru      = 79.04;
eps_eff_thru = 49.93;
eps_eff_thru_old = 49.93;
v_eff_thru   = c_0/sqrt(eps_eff_thru);
l_eff_thru   = v_eff_thru./f;
l0_eff_thru  = v_eff_thru./f0;
b_eff_thru   = 2*pi./l_eff_thru;
a_eff_thru   = b_eff_thru/2/Qi;
if lossesOnlyInRes
    g_eff_thru   = 1i*b_eff_thru;
else
    g_eff_thru   = a_eff_thru+1i*b_eff_thru;
end
l_thru = sep*l0_eff_thru;


%% TRX. LINE CONNECTING THRU AND RES - SHOULD BE ZERO
Z0_thru2res      = Z0_thru;
eps_eff_thru2res = eps_eff_thru;
v_eff_thru2res   = c_0/sqrt(eps_eff_thru2res);
l_eff_thru2res   = v_eff_thru2res./f;
l0_eff_thru2res  = v_eff_thru2res./f0;
b_eff_thru2res   = 2*pi./l_eff_thru2res;
a_eff_thru2res   = b_eff_thru2res/2/Qi;
if lossesOnlyInRes
    g_eff_thru2res   = 1i*b_eff_thru2res;
else
    g_eff_thru2res   = a_eff_thru2res+1i*b_eff_thru2res;
end
l_thru2res       = 0e-6*ones(1,nFilters);
    

%% KID
Z0_KID = 54.81;
 

%% FILTER-BANK TERMINATION
Z0_termination = Z0_thru;


%% RESONATOR PARAMETERS
Z0_res      = 54.81;
eps_eff_res = 47.28;
v_eff_res   = c_0/sqrt(eps_eff_res);
l_eff_res   = v_eff_res./f;
l0_eff_res  = v_eff_res./f0;
b_eff_res   = 2*pi./l_eff_res;
a_eff_res   = b_eff_res/2/Qi;
g_eff_res   = a_eff_res+1i*b_eff_res;
b0_eff_res  = 2*pi./l0_eff_res;
a0_eff_res  = b0_eff_res/2/Qi;
g0_eff_res  = a0_eff_res+1i*b0_eff_res;

% COUPLER C FOR A GIVEN QL
C_coup_thru = findCcoupler_Ql('SERIES','L2',Ql_target,Qi,f0,Z0_thru/2,Z0_res);
C_coup_KID  = findCcoupler_Ql('SERIES','L2',Ql_target,Qi,f0,Z0_res,Z0_KID);
% C_coup_thru = findCcoupler_Qc('SERIES','L2',Qc,f0,Z0_thru/2,Z0_res);
% C_coup_KID  = findCcoupler_Qc('SERIES','L2',Qc,f0,Z0_res,Z0_KID);

% RESONATOR LENGTH
% Faster option - Forces Zin_filter=Z0_thru/2 as input impedance of the resonators and thus forces the optimum of the filters
l_res = findL2res(Z0_thru,Z0_KID,Z0_res,eps_eff_res,f0,C_coup_KID,C_coup_thru);
% Slower option - Im(Zin)=0
% for fi=1:length(f0)
%     w0 = 2*pi*f0(fi);
%     Z_coup1 = -1i./(w0.*C_coup_KID(fi));
%     Z_coup2 = -1i./(w0.*C_coup_thru(fi));
%     A = Z_coup1 + Z0_KID;
%     syms l;
%     eqn = imag(Z_coup2+ Z0_res*(A+1i*Z0_res*tan(b0_eff_res(fi)*l))/(Z0_res+1i*A*tan(b0_eff_res(fi)*l)))==0;
%     l_res(fi) = double(vpasolve(eqn,l,[0,l0_eff_res(fi)/2]));
% end

%% Color for plots
colors = nipy_spectral(nFilters);


%% FILTER-BANK CALCULATIONS

% Initialization
ABCD_eachFilterBranch          = zeros(2,2,nFilters,nF);
ABCD_eachFilterUnit            = zeros(2,2,nFilters,nF);
ABCD_wholeFB                   = repmat(eye(2),1,1,nF);
ABCD_downToEachFilter          = zeros(2,2,nFilters+1,nF);
ABCD_downToEachFilter(:,:,1,:) = repmat(eye(2),1,1,nF);
ABCD_upToEachFilterBranch      = zeros(2,2,nFilters,nF);

% From antenna to termination
for filter_i=1:1:nFilters
    
    %% Parameters of each filter
    Z_coup_KID            = -1i./(w*C_coup_KID(filter_i));
    Z_coup_thru           = -1i./(w*C_coup_thru(filter_i));
    L_res                 = l_res(filter_i);
    L_thru2res            = l_thru2res(filter_i);
    L_toNextFilter        = l_thru(filter_i);
    
    %% ABCD matrix of each filter branch
    ABCD_eachFilterBranch(:,:,filter_i,:)  = filterBranch_ABCD(Z_coup_KID,Z0_res,g_eff_res,L_res,Z_coup_thru);
    
    %% ABCD matrix of each filter branch as a shunted load
    Zin_filter = Zin_fromABCD(squeeze(ABCD_eachFilterBranch(:,:,filter_i,:)),Z0_KID);
    ABCD_eachFilterBranchShuntLoad = shuntLoad_ABCD(Zin_filter);
    
%     figure(10); hold on;
%     niceplot(f/1e9,real(Zin_filter),'Color','b');
%     niceplot(f/1e9,imag(Zin_filter),'Color','r');
    
    %% ABCD matrix of each thru trx. line
    ABCD_eachThru = trxLine_ABCD(Z0_thru,g_eff_thru,L_toNextFilter);
    
    %% ABCD matrices of all filter units: filter+thru
    ABCD_upToEachFilterBranch(:,:,filter_i,:) = ABCD_eachThru;
    ABCD_eachFilterUnit(:,:,filter_i,:) = mmat(ABCD_eachFilterBranchShuntLoad,squeeze(ABCD_upToEachFilterBranch(:,:,filter_i,:)),[1,2]);
    
    %% ABCD matrices of the filters leading the FB down to each filter
    ABCD_downToEachFilter(:,:,filter_i+1,:) = mmat(ABCD_downToEachFilter(:,:,filter_i,:),ABCD_eachFilterUnit(:,:,filter_i,:),[1,2]);
    
    %% ABCD matrix of the whole FB
    ABCD_wholeFB = mmat(ABCD_wholeFB,squeeze(ABCD_eachFilterUnit(:,:,filter_i,:)),[1,2]);
    
end


% Initialization
S_eachFilter = zeros(2,2,nFilters,nF);
ABCD_upToEachFilter = zeros(2,2,nFilters+1,nF);
ABCD_upToEachFilter(:,:,nFilters+1,:) = repmat(eye(2),1,1,nF);
% From termination to antenna
for filter_i=nFilters:-1:1
        
    %% Calculate for the next loop iteration
    ABCD_upToEachFilter(:,:,filter_i,:) = mmat(ABCD_eachFilterUnit(:,:,filter_i,:),ABCD_upToEachFilter(:,:,filter_i+1,:),[1,2]);
    
    %%
    ABCD_upToEachFilterBranchFromFBEnd = mmat(ABCD_upToEachFilterBranch(:,:,filter_i,:),ABCD_upToEachFilter(:,:,filter_i+1,:),[1,2]);
    Zin = Zin_fromABCD(ABCD_upToEachFilterBranchFromFBEnd,Z0_termination);
    ABCD_upToEachFilterBranchShuntLoad = shuntLoad_ABCD(Zin);

    %% ABCD of each filter
    ABCD_temp       = mmat(ABCD_upToEachFilterBranchShuntLoad,squeeze(ABCD_eachFilterBranch(:,:,filter_i,:)),[1,2]);
    ABCD_eachFilter = mmat(squeeze(ABCD_downToEachFilter(:,:,filter_i,:)),ABCD_temp,[1,2]);
    
    S_eachFilter = a2s(ABCD_eachFilter,[Z0_thru,Z0_KID]);
    S31_eachFilter(filter_i,:) = squeeze(S_eachFilter(2,1,:));
    
    [Si1_f0_absSq(filter_i),f0_actual(filter_i),FWHM(filter_i),~] = findpeaks(abs(S31_eachFilter(filter_i,:)).^2,f,'WidthReference','HalfHeight','MinPeakProminence',0.5*max(abs(S31_eachFilter(filter_i,:)).^2),'NPeaks',1,'Annotate','extents');%NPeaks=1 takes the 1st one (low freq)
    Q_actual(filter_i)= f0_actual(filter_i)/FWHM(filter_i);
    f_aroundPeak = linspace(f0_actual(filter_i)-FWHM(filter_i)/2,f0_actual(filter_i)+FWHM(filter_i)/2,500);
    Si1_aroundPeak_absSq = interp1(f,abs(S31_eachFilter(filter_i,:)).^2,f_aroundPeak);
    df_aroundPeak = f_aroundPeak(end)-f_aroundPeak(end-1);
    df=f(end)-f(end-1);
    
    avgRelPower_inBand_eachFilter(filter_i)  = nansum(Si1_aroundPeak_absSq*df_aroundPeak) / FWHM(filter_i);
    avgRelPower_offBand_eachFilter(filter_i) = 1-avgRelPower_inBand_eachFilter(filter_i);
    
    relPower_inBandVSoffBand(filter_i) = nansum(Si1_aroundPeak_absSq*df_aroundPeak)./nansum(abs(S31_eachFilter(filter_i,:)).^2*df);
    
    relPower_inBandVStotal(filter_i)   = nansum(Si1_aroundPeak_absSq*df_aroundPeak)./nansum(abs(S31_eachFilter(filter_i,:)).^2*df);
    
    relPower_inBandVSinput(filter_i)   = nansum(Si1_aroundPeak_absSq*df_aroundPeak)./FWHM(filter_i);
    
end

% Calculate S-parameters
S_FB = a2s(ABCD_wholeFB,[Z0_thru,Z0_termination]);
S11_FB = squeeze(S_FB(1,1,:));
S21_FB = squeeze(S_FB(2,1,:));

% Calculate FB envelope
[envelopeAllChannels,indexMainChannelPerFreq] = max((abs(S31_eachFilter).^2),[],1);

% Total power extracted power by filters
totalExtractedPowerFilters=nansum(abs(S31_eachFilter).^2,1);


%% SANITY CHECKS
if plotSanityChecks
    test1 = abs(S11_FB).^2 + abs(S21_FB).^2 + nansum(abs(permute(S31_eachFilter,[2 1])).^2,2);
    test2 = nansum(abs(permute(S31_eachFilter,[2 1])).^2,2);
    test3 = 1-(abs(S11_FB).^2 + abs(S21_FB).^2);
    test4 = test2-test3;
    test5 = test4+1;
    figure(); hold on;
    niceplot(f/1e9,100*test1);
    niceplot(f/1e9,100*test2);
    niceplot(f/1e9,100*test3);
    niceplot(f/1e9,100*test4);
    niceplot(f/1e9,100*test5);
    ylabel('%');

    figure(); hold on;
    niceplot(f0/1e9,f0_actual/1e9);
    xlabel('Designed f  [GHz]','Interpreter','LaTeX');
    ylabel('Measured f  [\%]','Interpreter','LaTeX');
    axis square;
end


%% PLOTTING S-PARAMETERS
fig=figure();
fig.WindowState = 'maximized';
hold on;

if plotAllFilters
    for filter_i=1:1:nFilters
        nicearea(f/1e9,10*log10(abs(S31_eachFilter(filter_i,:)).^2),'FaceColor',colors(filter_i,:),'FaceAlpha',1,'Basevalue',-Inf,'EdgeColor',[0,0,0],'Linewidth',0.1);
    end
end
hLegend(1) = niceplot(f/1e9,10*log10(abs(S11_FB).^2),'Color',[0 1 1],'LineWidth',1);
hLegend(2) = niceplot(f/1e9,10*log10(abs(S21_FB).^2),'Color',[1 0 1],'LineWidth',1);
hLegend(3) = niceplot(f/1e9,10*log10(totalExtractedPowerFilters),'Color',[0.5 0.5 0.5],'LineWidth',1);
lgnd = {'$|S_{11}|^2\quad$','$|S_{21}|^2\quad$','$\Sigma_{i=3}^{N+2}|S_{i1}|^2\quad$'};

if plotEnvelope
    hLegend(4) = niceplot(f/1e9,10*log10(envelopeAllChannels),'Color',[0 0 0],'LineWidth',1);
    lgnd = {'$|S_{11}|^2\quad$','$|S_{21}|^2\quad$','$1 - (|S_{11}|^2+|S_{21}|^2)\quad$','$|S_{i1}|^2\,\mathrm{envelope}$'};
end
lgd=legend(hLegend,lgnd,'Interpreter','LaTeX','Location','NorthOutside','Orientation','Horizontal');
xlim([min(f),max(f)]/1e9);
ylim([mindB,0]);
xlabel('Frequency [GHz]','Interpreter','LaTeX');
ylabel('S-parameters [dB]','Interpreter','LaTeX');


%% PLOTTING Q-FACTORS
if plotQfactors
    figure();
    hold on;
    niceplot(f0_actual/1e9,Q_actual,'LineStyle','none','LineWidth',2,'Marker','o','MarkerSize',3,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
    xlim([min(f),max(f)]/1e9);
    ylim([0,2*Ql_target]);
    xlabel('$f_i$ [GHz]','Interpreter','LaTeX');
    ylabel('$Q_l$','Interpreter','LaTeX');
    
end


%% PLOTTING INTER-FILTER SPACING
if plotLthru
    figure();
    hold on;
    niceplot(f0/1e9,l_thru/1e-6,'LineStyle','none','LineWidth',2,'Marker','o','MarkerSize',3,'MarkerFaceColor',[0,0,0]);
    xlim([min(f),max(f)]/1e9);
    xlabel('Frequency  [GHz]','Interpreter','LaTeX');
    ylabel('$L_\mathrm{thru}$ [$\mu{}m$]','Interpreter','LaTeX');   
end


%% PLOTTING OVERLAPPING PEAKS
if plotOverlappingPeaks
    figure;
    hold on;
    for filter_i=1:1:nFilters
        niceplot(f/f0(filter_i),10*log10(abs(S31_eachFilter(filter_i,:)).^2),'Color',colors(filter_i,:));
    end
    xlim(1+10*[-1/2,+1/2]*1/Ql_target);
    ylim([mindB,0]);
    xlabel('$f/f_i$','Interpreter','LaTeX');
    ylabel('S-parameters  [dB]','Interpreter','LaTeX');
end


%% POWER IN-BAND VS TOTAL INPUT POWER AND QL
colors=lines(2);
figure();
yyaxis left;
hold on;
niceplot(f0/1e9,100*relPower_inBandVStotal,'LineStyle','none','LineWidth',2,'Marker','o','MarkerSize',3,'MarkerFaceColor',colors(1,:));
xlim([min(f),max(f)]/1e9);
ylim([0,100]);
xlabel('$f_i$ [GHz]','Interpreter','LaTeX');
ylabel('$P_\mathrm{in}/P_\mathrm{tot}$ [\%]','Interpreter','LaTeX');
yyaxis right;
niceplot(f0/1e9,Q_actual,'LineStyle','none','LineWidth',2,'Marker','o','MarkerSize',3,'MarkerFaceColor',colors(2,:));
ylim([0,2*Ql_target]);
ylabel('$Q_l$','Interpreter','LaTeX');


%% AVERAGE POWER IN-BAND VS AVERAGE POWER OFF-BAND
figure(); hold on;
nicebar(f0/1e9,100*[avgRelPower_inBand_eachFilter;avgRelPower_offBand_eachFilter],'BarLayout','stacked','barwidth',1);
niceplot(f0/1e9,100*avgRelPower_inBand_eachFilter,'Color',[0,0,0]);
ylim([0,100]);
xlabel('$f_i$ [GHz]','Interpreter','LaTeX');
ylabel('Avg power [\%]','Interpreter','LaTeX');

%% POWER IN-BAND VS TOTAL INPUT POWER
figure();
niceplot(f0/1e9,10*log10(relPower_inBandVStotal),'LineStyle','none','LineWidth',2,'Marker','o','MarkerSize',3,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
xlim([min(f),max(f)]/1e9);
ylim([-15,0]);
xlabel('$f_i$ [GHz]','Interpreter','LaTeX');
ylabel('$P_\mathrm{in}/P_\mathrm{tot}$ [dB]','Interpreter','LaTeX');



%% POWER IN-BAND VS AVAILABLE POWER IN-BAND. CF ISOLATED FILTER
S31_absSq_max_isolated=(Qi-Ql_target)^2/(2*Qi^2);
if ~isfinite(Qi) && ~isnan(Qi)
    S31_absSq_max_isolated=0.5;
end
relPower_inBandVSinput_isolated = S31_absSq_max_isolated*pi/4;

figure();
hold on;
niceplot(f0/1e9,10*log10(relPower_inBandVSinput),'LineStyle','none','LineWidth',2,'Marker','o','MarkerSize',3,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none');
niceplot([min(f),max(f)]/1e9,10*log10(relPower_inBandVSinput_isolated*[1,1]),'LineStyle',':','Color',0.7*[1,1,1]);
xlim([min(f),max(f)]/1e9);
ylim([-15,0]);
xlabel('$f_i$ [GHz]','Interpreter','LaTeX');
ylabel('$\eta_{\mathrm{filter},i}^\mathrm{FWHM}$ [dB]','Interpreter','LaTeX');


%% MEAN OF FIGURES OF MERIT
S31_mean_dB = 10*log10(nanmean(Si1_f0_absSq(filter_i)))

[~,f_1_i] = min(abs(f-f0(1)));
[~,f_N_i] = min(abs(f-f0(end)));
totalExtractedPowerFilters_mean_dB = 10*log10(nanmean(totalExtractedPowerFilters(f_N_i:f_1_i)))

Q_actual_mean = nanmean(Q_actual)

relPower_inBandVStotal_dB = 10*log10(nanmean(relPower_inBandVStotal))

relPower_inBandVSinput_mean_dB     = 10*log10(nanmean(relPower_inBandVSinput))
relPower_inBandVSinput_isolated_dB = 10*log10(relPower_inBandVSinput_isolated)

