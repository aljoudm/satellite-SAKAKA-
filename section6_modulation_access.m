%% SECTION 6 - MODULATION AND MULTIPLE ACCESS
%  GEO Ku-band link (satellite at 45 deg E) for rural broadband to Al-Jouf,
%  Saudi Arabia (reference site: Sakaka).
%
%  Forward link : DVB-S2, 8PSK code rate 3/4, Rs = 30 Msym/s, RRC alpha = 0.20
%  Return link  : MF-TDMA, 14 carriers x 2 Msym/s, QPSK 3/4
%
%  Requires MATLAB + Communications Toolbox + Signal Processing Toolbox
%  (no Satellite Communications Toolbox needed).
%
%  Outputs: figures/Fig6_1 ... Fig6_6 (PNG, 300 dpi),
%           results/summary.txt and results/results.mat
%
%  Usage: run this file from any folder; paths are resolved relative to it.
clear; close all; clc;
rng(2026);                                   % reproducible random numbers

%% 0. Settings, design values and paths
runLDPC = false;                             % optional LDPC check (slow) - set true to run

projDir = fileparts(mfilename('fullpath'));
figDir  = fullfile(projDir, 'figures');
resDir  = fullfile(projDir, 'results');
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(resDir, 'dir'), mkdir(resDir); end

% ---- Fixed design values (do not change) ----
Rs        = 30e6;                            % symbol rate [sym/s]
alpha     = 0.20;                            % RRC roll-off
Bxp       = 36e6;                            % transponder bandwidth [Hz]
EsN0_clr  = 17.35;                           % overall Es/N0, clear sky [dB]
EsN0_999  = 16.44;                           % overall Es/N0, 99.9 %  availability [dB]
EsN0_9999 = 14.72;                           % overall Es/N0, 99.99 % availability [dB]
implMargin = 1.5;                            % implementation margin [dB]

% DVB-S2 ideal Es/N0 (ETSI EN 302 307-1, Table 13) and spectral efficiency
modName = {'QPSK 1/2','QPSK 3/4','8PSK 2/3','8PSK 3/4','8PSK 5/6'};
modEsN0 = [1.00, 4.03, 6.62, 7.91, 9.35];            % [dB]
modEta  = [0.988858, 1.487473, 1.980636, 2.228124, 2.478562];  % [bit/symbol]
chosen  = 4;                                 % index of the chosen design: 8PSK 3/4

% Return link (MF-TDMA)
nCarriers = 14;  carrierRs = 2e6;  carrierRate = 2.97e6;  totalRate = 41.6e6;

% ---- Plot style (formal report) ----
navy = [0.00 0.20 0.50];                     % 8PSK / chosen design
gry  = [0.55 0.55 0.55];                     % alternatives
blk  = [0 0 0];                              % theory
set(groot, 'defaultAxesFontName','Arial', 'defaultAxesFontSize',11, ...
    'defaultTextFontName','Arial', 'defaultLegendFontName','Arial', ...
    'defaultLineLineWidth',1.5, 'defaultFigureColor','w', 'defaultAxesColor','w');

R = struct();                                % container for saved results

%% 1. Figure 6.1 - Received QPSK and 8PSK constellations (Es/N0 = 16.44 dB)
nSym = 5000;
mods = [4 8];  phs = [pi/4 0];  names1 = {'QPSK','8PSK'};
fig = newFig([100 100 1000 480]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
for i = 1:2
    M = mods(i);
    d  = randi([0 M-1], nSym, 1);
    tx = pskmod(d, M, phs(i), 'gray');                   % unit average energy
    N0 = 10^(-EsN0_999/10);                              % Es = 1 -> N0 = 1/(Es/N0)
    rx = tx + sqrt(N0/2)*(randn(nSym,1) + 1i*randn(nSym,1));   % AWGN
    ideal = pskmod(0:M-1, M, phs(i), 'gray');
    ax = nexttile(tl);  hold(ax,'on');
    col = gry; if M == 8, col = navy; end
    plot(ax, real(rx), imag(rx), '.', 'Color', col, 'MarkerSize', 5, 'DisplayName','Received symbols');
    plot(ax, real(ideal), imag(ideal), 'o', 'MarkerSize', 11, 'LineWidth',1.5, ...
        'MarkerFaceColor','w', 'MarkerEdgeColor', [0.85 0.1 0.1], 'DisplayName','Ideal points');
    axis(ax,'equal'); xlim(ax,[-1.5 1.5]); ylim(ax,[-1.5 1.5]);
    xlabel(ax,'In-phase'); ylabel(ax,'Quadrature');
    title(ax, sprintf('%s, E_s/N_0 = %.2f dB (%d symbols)', names1{i}, EsN0_999, nSym));
    legend(ax,'Location','northeast'); styleAx(ax);
end
exportgraphics(fig, fullfile(figDir,'Fig6_1_constellation.png'), 'Resolution',300);

%% 2. Figure 6.2 - Uncoded BER vs Eb/N0 (theory + Monte-Carlo)
EbN0_th = 0:0.05:16;                         % fine grid for theory curves [dB]
EbN0_mc = 0:1:16;                            % Monte-Carlo points [dB]
minBits = 1e6;  maxBits = 2e7;  errStop = 100;   % >=1e6 bits; stop after 100 errors
pc = [0 1 1 2 1 2 2 3];                      % bit-population count lookup (3-bit labels)
ber_th = zeros(2, numel(EbN0_th));  ber_mc = nan(2, numel(EbN0_mc));
cross  = zeros(1,2);
for i = 1:2
    M = mods(i);  k = log2(M);
    ber_th(i,:) = berawgn(EbN0_th, 'psk', M, 'nondiff');   % Gray-coded theory
    % 1e-6 crossing by log-domain interpolation
    [bu, iu] = unique(log10(ber_th(i,:)));
    cross(i) = interp1(bu, EbN0_th(iu), -6);
    for j = 1:numel(EbN0_mc)
        EsN0 = EbN0_mc(j) + 10*log10(k);              % Es/N0 = k * Eb/N0
        sig  = sqrt(0.5 * 10^(-EsN0/10));             % noise std per dimension (Es=1)
        errs = 0; bits = 0; chunk = ceil(5e5/k);
        while true
            d  = randi([0 M-1], chunk, 1);
            rx = pskmod(d, M, phs(i), 'gray') + sig*(randn(chunk,1) + 1i*randn(chunk,1));
            dh = pskdemod(rx, M, phs(i), 'gray');
            errs = errs + sum(pc(bitxor(d, dh) + 1));
            bits = bits + chunk*k;
            if (bits >= minBits && errs >= errStop) || bits >= maxBits, break; end
        end
        if errs > 0, ber_mc(i,j) = errs/bits; end    % zero-error points are not plotted
    end
end
fprintf('\n--- Uncoded BER = 1e-6 crossing (theory) ---\n');
fprintf('QPSK : Eb/N0 = %.2f dB\n8PSK : Eb/N0 = %.2f dB\n', cross(1), cross(2));
R.berCrossing_QPSK_dB = cross(1);  R.berCrossing_8PSK_dB = cross(2);

fig = newFig([100 100 760 560]);  ax = axes(fig);  hold(ax,'on');
semilogy(ax, EbN0_th, ber_th(1,:), '-',  'Color',blk, 'DisplayName','QPSK theory');
semilogy(ax, EbN0_th, ber_th(2,:), '--', 'Color',blk, 'DisplayName','8PSK theory');
semilogy(ax, EbN0_mc, ber_mc(1,:), 's', 'Color',gry,  'MarkerFaceColor','w', 'MarkerSize',7, 'DisplayName','QPSK Monte-Carlo');
semilogy(ax, EbN0_mc, ber_mc(2,:), 'o', 'Color',navy, 'MarkerFaceColor',navy,'MarkerSize',6, 'DisplayName','8PSK Monte-Carlo');
yline(ax, 1e-6, ':', 'Color',[0.75 0.1 0.1], 'LineWidth',1.5, 'DisplayName','BER = 10^{-6}');
set(ax,'YScale','log'); ylim(ax,[1e-7 1]); xlim(ax,[0 16]);
xlabel(ax,'E_b/N_0 (dB)'); ylabel(ax,'Bit error rate (BER)');
title(ax,'Uncoded BER of Gray-coded QPSK and 8PSK in AWGN');
legend(ax,'Location','southwest'); styleAx(ax);
exportgraphics(fig, fullfile(figDir,'Fig6_2_BER.png'), 'Resolution',300);

%% 3. Figure 6.3 - Spectral efficiency vs Eb/N0 (Shannon bound and MODCODs)
eta = linspace(0.3, 3, 300);
EbN0_sh_dB = 10*log10((2.^eta - 1)./eta);            % Shannon bound
EbN0_mod   = modEsN0 - 10*log10(modEta);             % Eb/N0 = Es/N0 - 10log10(eta)
EbN0_min   = 10*log10((2.^modEta - 1)./modEta);      % Shannon limit at each MODCOD's eta
gap        = EbN0_mod - EbN0_min;                    % gap to Shannon [dB]
fprintf('\n--- Gap to Shannon bound ---\n');
for m = 1:5
    fprintf('%-9s eta = %.4f  Eb/N0 = %5.2f dB  Shannon min = %5.2f dB  gap = %4.2f dB\n', ...
        modName{m}, modEta(m), EbN0_mod(m), EbN0_min(m), gap(m));
end
R.shannon_min_8PSK34_dB = EbN0_min(chosen);  R.gap_dB = gap;  R.EbN0_modcod_dB = EbN0_mod;

fig = newFig([100 100 760 560]);  ax = axes(fig);  hold(ax,'on');
plot(ax, EbN0_sh_dB, eta, '-', 'Color',blk, 'DisplayName','Shannon bound');
for m = 1:5
    if m == chosen
        plot(ax, EbN0_mod(m), modEta(m), 'o', 'MarkerSize',11, 'MarkerFaceColor',navy, ...
            'MarkerEdgeColor',navy, 'DisplayName','DVB-S2 8PSK 3/4 (design)');
    else
        plot(ax, EbN0_mod(m), modEta(m), 's', 'MarkerSize',8, 'MarkerFaceColor',gry, ...
            'MarkerEdgeColor',gry, 'HandleVisibility','off');
    end
    text(ax, EbN0_mod(m)+0.25, modEta(m)-0.06, modName{m}, 'FontSize',10, ...
        'FontWeight', ternary(m==chosen,'bold','normal'), 'Color', ternary(m==chosen,navy,[0.3 0.3 0.3]));
end
plot(ax, nan, nan, 's', 'MarkerFaceColor',gry, 'MarkerEdgeColor',gry, 'MarkerSize',8, ...
    'DisplayName','Other DVB-S2 MODCODs');
xlim(ax,[-2 7]); ylim(ax,[0 3.2]);
xlabel(ax,'E_b/N_0 (dB)'); ylabel(ax,'Spectral efficiency (bit/s/Hz)');
title(ax,'DVB-S2 MODCODs against the Shannon bound');
legend(ax,'Location','southeast'); styleAx(ax);
exportgraphics(fig, fullfile(figDir,'Fig6_3_shannon.png'), 'Resolution',300);

%% 4. Figure 6.4 - Link margin per MODCOD
reqEsN0   = modEsN0 + implMargin;                    % required Es/N0 = ideal + implementation margin
margin999 = EsN0_999  - reqEsN0;
margin9999= EsN0_9999 - reqEsN0;
fprintf('\n--- Margin per MODCOD (required = ideal + %.1f dB) ---\n', implMargin);
for m = 1:5
    fprintf('%-9s required Es/N0 = %5.2f dB  margin@99.9%% = %5.2f dB  margin@99.99%% = %5.2f dB\n', ...
        modName{m}, reqEsN0(m), margin999(m), margin9999(m));
end
R.reqEsN0_dB = reqEsN0;  R.margin999_dB = margin999;  R.margin9999_dB = margin9999;

fig = newFig([100 100 820 520]);  ax = axes(fig);  hold(ax,'on');
b = barh(ax, 1:5, reqEsN0, 0.6, 'FaceColor','flat', 'EdgeColor','none');
b.CData = repmat(gry, 5, 1);  b.CData(chosen,:) = navy;
xline(ax, EsN0_999,  '--', 'Color',[0.1 0.1 0.1], 'LineWidth',1.5, ...
    'Label',sprintf('99.9%%: %.2f dB',EsN0_999),  'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','right');
xline(ax, EsN0_9999, '--', 'Color',[0.75 0.1 0.1], 'LineWidth',1.5, ...
    'Label',sprintf('99.99%%: %.2f dB',EsN0_9999), 'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left');
for m = 1:5
    text(ax, reqEsN0(m)+0.2, m, sprintf('margin %.1f dB', margin999(m)), ...
        'VerticalAlignment','middle', 'FontSize',10, ...
        'FontWeight', ternary(m==chosen,'bold','normal'), 'Color', ternary(m==chosen,navy,[0.2 0.2 0.2]));
end
set(ax,'YTick',1:5,'YTickLabel',modName,'YDir','reverse'); xlim(ax,[0 20]);
xlabel(ax,'Required E_s/N_0 = ideal + 1.5 dB (dB)'); ylabel(ax,'DVB-S2 MODCOD');
title(ax,'Required E_s/N_0 and margin at 99.9% availability');
styleAx(ax);
exportgraphics(fig, fullfile(figDir,'Fig6_4_margin.png'), 'Resolution',300);

%% 5. Figure 6.5 - Transmitted 8PSK spectrum (RRC, alpha = 0.20, 30 Msym/s)
sps  = 8;  fs = Rs*sps;                              % 240 MHz sampling rate
nS   = 2e5;
rrc  = rcosdesign(alpha, 10, sps, 'sqrt');           % root-raised-cosine, span 10 symbols
sym  = pskmod(randi([0 7], nS, 1), 8, 0, 'gray');
x    = upfirdn(sym, rrc, sps);                       % pulse-shaped baseband signal
[pxx, f] = pwelch(x, hann(4096), 2048, 4096, fs, 'centered');
pdb  = 10*log10(pxx/max(pxx));                       % peak normalised to 0 dB
f_MHz = f/1e6;
[obwHz, ~, ~] = obw(x, fs, [], 99);                  % 99 % occupied bandwidth
fprintf('\n--- Spectrum ---\n99%% occupied bandwidth = %.2f MHz (theory Rs(1+alpha) = %.1f MHz)\n', ...
    obwHz/1e6, Rs*(1+alpha)/1e6);
R.obw99_MHz = obwHz/1e6;
inBand = 100*sum(pxx(abs(f) <= Bxp/2))/sum(pxx);     % share of power inside +/-18 MHz
fprintf('Power inside +/-%.0f MHz (transponder) = %.2f %%\n', Bxp/2e6, inBand);
R.powerInTransponder_pct = inBand;

fig = newFig([100 100 820 500]);  ax = axes(fig);  hold(ax,'on');
plot(ax, f_MHz, pdb, '-', 'Color',navy, 'DisplayName','8PSK, RRC \alpha = 0.20');
xline(ax,  Bxp/2e6, '--', 'Color',[0.75 0.1 0.1], 'LineWidth',1.5, 'DisplayName','Transponder edges (\pm18 MHz)');
xline(ax, -Bxp/2e6, '--', 'Color',[0.75 0.1 0.1], 'LineWidth',1.5, 'HandleVisibility','off');
xlim(ax,[-60 60]); ylim(ax,[-80 5]);
xlabel(ax,'Frequency offset from carrier centre (MHz)'); ylabel(ax,'Normalised PSD (dB)');
title(ax, sprintf('Transmitted 8PSK spectrum at 30 Msym/s (99%% occupied BW = %.1f MHz)', obwHz/1e6));
legend(ax,'Location','northeast'); styleAx(ax);
exportgraphics(fig, fullfile(figDir,'Fig6_5_spectrum.png'), 'Resolution',300);

%% 6. Figure 6.6 - MF-TDMA demand-assignment simulation (return link)
nHomes = 300;  pAct = 0.05;  reqMbps = 2;  capMbps = totalRate/1e6;   % 41.6 Mbps
nFrames = 10000;
capBursts = floor(capMbps/reqMbps);                  % 20 bursts of ~2 Mbps per frame
nSlots = 8;  slotsPerBurst = floor(nCarriers*nSlots/capBursts);   % 112 cells / 20 = 5 cells per burst

active   = rand(nHomes, nFrames) < pAct;             % activity of every home in every frame
nActive  = sum(active, 1);                           % active homes per frame
demand   = nActive * reqMbps;                        % offered demand [Mbps]
carried  = min(demand, capMbps);                     % traffic actually served
meanDemand = mean(demand);
utilOffer  = 100*meanDemand/capMbps;                 % offered load / capacity
utilCarry  = 100*mean(carried)/capMbps;              % carried traffic / capacity
pExceed    = 100*mean(demand > capMbps);
p99Demand  = prctile(demand, 99);   p99Homes = prctile(nActive, 99);
fdmaMbps   = nHomes*reqMbps;                         % fixed FDMA: one channel per home
fprintf('\n--- MF-TDMA demand assignment (%d homes, p = %.2f, %d frames) ---\n', nHomes, pAct, nFrames);
fprintf('Mean active homes        : %.2f\n', mean(nActive));
fprintf('Mean demand              : %.2f Mbps\n', meanDemand);
fprintf('Mean utilisation (offered): %.1f %% of %.1f Mbps (carried: %.1f %%)\n', utilOffer, capMbps, utilCarry);
fprintf('Frames with demand > cap : %.2f %%\n', pExceed);
fprintf('99th-percentile demand   : %.1f Mbps (%.0f homes)\n', p99Demand, p99Homes);
fprintf('Fixed FDMA requirement   : %d x %d = %d Mbps (%.1fx the MF-TDMA capacity)\n', ...
    nHomes, reqMbps, fdmaMbps, fdmaMbps/capMbps);
R.tdma = struct('meanActive',mean(nActive),'meanDemand_Mbps',meanDemand, ...
    'utilOffered_pct',utilOffer,'utilCarried_pct',utilCarry,'pExceed_pct',pExceed, ...
    'p99Demand_Mbps',p99Demand,'p99Homes',p99Homes,'fdma_Mbps',fdmaMbps,'capacity_Mbps',capMbps);

% Example frame: first frame whose active-home count equals the mean (rounded)
fEx = find(nActive == round(mean(nActive)), 1);
homesEx = find(active(:, fEx));
grid2 = zeros(nCarriers, nSlots);                    % 0 = empty cell
cellIdx = 1;                                            % illustrative fill: bursts placed consecutively
for h = 1:min(numel(homesEx), capBursts)
    for s = 1:slotsPerBurst
        [si, ci] = ind2sub([nSlots nCarriers], cellIdx);   % fill along time on each carrier
        grid2(ci, si) = homesEx(h);  cellIdx = cellIdx + 1;
    end
end
fig = newFig([100 100 1150 500]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
ax = nexttile(tl);
idx = zeros(size(grid2));  [~,~,ic] = unique(grid2(grid2>0));  idx(grid2>0) = ic;
cm = [1 1 1; interp1([0 0.5 1], [0.75 0.85 1; 0.35 0.55 0.85; 0.0 0.2 0.5], linspace(0,1,max(2,max(idx(:)))))];
shades = cm(2:end,:);  cm(2:end,:) = shades(randperm(size(shades,1)),:);   % shuffle so neighbours differ
image(ax, 1:nSlots, 1:nCarriers, idx+1); colormap(ax, cm);
set(ax,'YDir','normal'); axis(ax,'tight'); hold(ax,'on');
for c = 0.5:1:nCarriers+0.5, yline(ax,c,'-','Color',[0.8 0.8 0.8],'LineWidth',0.5); end
for s = 0.5:1:nSlots+0.5,    xline(ax,s,'-','Color',[0.8 0.8 0.8],'LineWidth',0.5); end
for ci = 1:nCarriers
    for si = 1:nSlots
        if grid2(ci,si) > 0
            text(ax, si, ci, num2str(grid2(ci,si)), 'HorizontalAlignment','center', 'FontSize',8, 'Color','k');
        end
    end
end
set(ax,'XTick',1:nSlots,'YTick',1:nCarriers);
xlabel(ax,'Time slot within frame'); ylabel(ax,'Carrier number (2 Msym/s each)');
title(ax, sprintf('(a) Example frame: %d active homes, numbers = home ID', numel(homesEx)));
set(ax,'Box','on','Layer','top');

ax = nexttile(tl);  hold(ax,'on');
histogram(ax, nActive, -0.5:1:max(nActive)+0.5, 'Normalization','probability', ...
    'FaceColor',navy, 'EdgeColor','w', 'DisplayName','Active homes per frame');
xline(ax, capMbps/reqMbps, '--', 'Color',[0.75 0.1 0.1], 'LineWidth',1.5, ...
    'DisplayName', sprintf('Capacity = %.1f Mbps (%d bursts)', capMbps, capBursts));
xlabel(ax,'Active homes per frame'); ylabel(ax,'Fraction of frames');
title(ax, sprintf('(b) Demand per frame: exceeds capacity in %.1f%% of frames', pExceed));
legend(ax,'Location','northeast'); styleAx(ax);
exportgraphics(fig, fullfile(figDir,'Fig6_6_mftdma.png'), 'Resolution',300);

%% 7. Optional - DVB-S2 LDPC check for 8PSK 3/4 (frame error rate)
if runLDPC
    H    = dvbsLDPCPCM(3/4);                         % normal frame, 64800 bits
    cfgE = ldpcEncoderConfig(H);  cfgD = ldpcDecoderConfig(H);
    EsN0_ldpc = 7.0:0.25:9.0;  nFrm = 30;  fer = zeros(size(EsN0_ldpc));
    for j = 1:numel(EsN0_ldpc)
        N0 = 10^(-EsN0_ldpc(j)/10);  nfe = 0;
        for fr = 1:nFrm
            msg = randi([0 1], cfgE.NumInformationBits, 1, 'int8');
            cw  = ldpcEncode(msg, cfgE);
            tx  = pskmod(double(cw), 8, 0, 'gray', 'InputType','bit');
            rx  = tx + sqrt(N0/2)*(randn(size(tx)) + 1i*randn(size(tx)));
            llr = pskdemod(rx, 8, 0, 'gray', 'OutputType','approxllr', 'NoiseVariance',N0);
            dec = ldpcDecode(llr, cfgD, 50);
            nfe = nfe + any(dec ~= msg);
        end
        fer(j) = nfe/nFrm;
        fprintf('LDPC Es/N0 = %.2f dB  FER = %.3f\n', EsN0_ldpc(j), fer(j));
    end
    fig = newFig([100 100 760 520]);  ax = axes(fig);  hold(ax,'on');
    plot(ax, EsN0_ldpc, max(fer,1/(2*nFrm)), 'o-', 'Color',navy, 'MarkerFaceColor',navy, 'DisplayName','Simulated FER (LDPC only)');
    xline(ax, 7.91, '--', 'Color',blk, 'LineWidth',1.5, 'DisplayName','ETSI ideal 7.91 dB');
    set(ax,'YScale','log'); ylim(ax,[1/(4*nFrm) 1.5]);
    xlabel(ax,'E_s/N_0 (dB)'); ylabel(ax,'Frame error rate');
    title(ax,'DVB-S2 8PSK 3/4 LDPC check (64800-bit frames, no BCH)');
    legend(ax,'Location','southwest'); styleAx(ax);
    exportgraphics(fig, fullfile(figDir,'Fig6_7_ldpc_check.png'), 'Resolution',300);
    R.ldpc_EsN0_dB = EsN0_ldpc;  R.ldpc_FER = fer;
end

%% 8. Summary table (printed and saved)
L = {};
L{end+1} = '================ SECTION 6 SUMMARY: MODULATION AND MULTIPLE ACCESS ================';
L{end+1} = sprintf('Forward link: DVB-S2 8PSK 3/4, eta = %.6f bit/sym, Rs = %.0f Msym/s, alpha = %.2f', modEta(chosen), Rs/1e6, alpha);
L{end+1} = sprintf('Link Es/N0: clear %.2f dB | 99.9%% %.2f dB | 99.99%% %.2f dB | implementation margin %.1f dB', EsN0_clr, EsN0_999, EsN0_9999, implMargin);
L{end+1} = '';
L{end+1} = sprintf('%-10s %6s %9s %9s %9s %10s %10s %9s', 'MODCOD','eta','Ideal','Required','Eb/N0','Gap to','Margin','Margin');
L{end+1} = sprintf('%-10s %6s %9s %9s %9s %10s %10s %9s', '','[b/s]','Es/N0 dB','Es/N0 dB','dB','Shannon dB','99.9% dB','99.99% dB');
for m = 1:5
    L{end+1} = sprintf('%-10s %6.3f %9.2f %9.2f %9.2f %10.2f %10.2f %9.2f', modName{m}, modEta(m), ...
        modEsN0(m), reqEsN0(m), EbN0_mod(m), gap(m), margin999(m), margin9999(m)); %#ok<*AGROW>
end
L{end+1} = '';
L{end+1} = sprintf('Uncoded BER = 1e-6 crossing (theory): QPSK %.2f dB, 8PSK %.2f dB (Eb/N0)', cross(1), cross(2));
L{end+1} = sprintf('Shannon minimum Eb/N0 at eta = %.3f: %.2f dB; 8PSK 3/4 gap = %.2f dB', modEta(chosen), EbN0_min(chosen), gap(chosen));
L{end+1} = sprintf('99%% occupied bandwidth: %.2f MHz; nominal Rs(1+alpha) = %.1f MHz; %.2f %% of power inside the %.0f MHz transponder', obwHz/1e6, Rs*(1+alpha)/1e6, inBand, Bxp/1e6);
L{end+1} = '';
L{end+1} = sprintf('MF-TDMA return link: %d carriers x %.0f Msym/s, QPSK 3/4, %.2f Mbps/carrier, %.1f Mbps total', nCarriers, carrierRs/1e6, carrierRate/1e6, totalRate/1e6);
L{end+1} = sprintf('  %d homes, p(active) = %.2f, %d Mbps each, %d frames', nHomes, pAct, reqMbps, nFrames);
L{end+1} = sprintf('  Mean demand            : %.2f Mbps (%.2f active homes)', meanDemand, mean(nActive));
L{end+1} = sprintf('  Mean utilisation       : %.1f %% offered, %.1f %% carried', utilOffer, utilCarry);
L{end+1} = sprintf('  Frames over capacity   : %.2f %%', pExceed);
L{end+1} = sprintf('  99th-percentile demand : %.1f Mbps (%.0f homes)', p99Demand, p99Homes);
L{end+1} = sprintf('  Fixed FDMA requirement : %d Mbps vs %.1f Mbps shared (%.1fx saving)', fdmaMbps, capMbps, fdmaMbps/capMbps);
L{end+1} = '===================================================================================';
fprintf('\n%s\n', strjoin(L, newline));
fid = fopen(fullfile(resDir,'summary.txt'), 'w');
fprintf(fid, '%s\n', strjoin(L, newline));  fclose(fid);
R.modName = modName;  R.modEsN0 = modEsN0;  R.modEta = modEta;
save(fullfile(resDir,'results.mat'), '-struct', 'R');
fprintf('\nFigures saved in %s\nResults saved in %s\n', figDir, resDir);

%% Local functions
function fig = newFig(pos)
% Create a white figure of given position/size (invisible friendly for batch runs)
fig = figure('Color','w', 'Position',pos);
end

function styleAx(ax)
% Report style: Arial 11, light major grid only, boxed axes
set(ax, 'FontName','Arial', 'FontSize',11, 'Box','on', 'Color','w', ...
    'XGrid','on', 'YGrid','on', 'XMinorGrid','off', 'YMinorGrid','off', ...
    'GridColor',[0.8 0.8 0.8], 'GridAlpha',1, 'LineWidth',0.75);
end

function out = ternary(c, a, b)
% Inline conditional helper
if c, out = a; else, out = b; end
end

