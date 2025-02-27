aesthetics
Jfac = 1e9; % Conversion from mol/s to nmol/s
rp = logspace(-9, -4.5, 400)'; % Fracture sizes, m
Jgc = 400e-9; % kg/m^2/s, around 400 nmol/m^2/s

r = 1.5e-5; % Cell radius
vbar = 1.805e-5; % Molar volume of water
RT = 298*8.314; % RT thermal energy scale
psi_g = -1.5; % Gland water potential
cg = 100; % Gland concentration, mM
vfac = vbar*8.64e6/(pi*r^2); % mol/s to mL/cm^2/d
cstar = 1400; % Gradient threshold, mM

% Run simulation
sa = saltgland(rp, Jgc, 'n', 1, 'r', r);

cb = sa.cb; % Brine concentration, mmol/kg
c = sa.c; % M
Pc = sa.Pc; % MPa
psi_c = sa.psi_c; % MPa
Jtot = sa.Jtot*vfac; % mL/cm^2/day
Deff = sa.Deff;
Jsalt = sa.Jsalt;
JD = sa.JD; % umol/day
JP = sa.JP; % umol/day
psi_b = -cb*RT*1e-6;


Jgc = (Jgc*1e-23)*(1e4)^2*Jfac; % umol/m^2/day

Rp = (rp*1e9); % nm
xlab = 'Fracture size (nm)';

% Transporter activity range
rt = Rp((c-(100/cb)-0.7) <= 0);
rt = rt([1 end]);

% Velocities
vel = @(J) J./(pi*rp.^2*1e12)*1e9;

% Find point of max adv velocity
rv = Rp(vel(JP) == max(vel(JP)));

%% water loss
try
    figure(rpJtotfig);
    clf
catch
    rpJtotfig = figure();
end

rpJtotfig.Position(3:4)=figsize;

hold on
j1 = loglog(Rp, Jtot(:,1), 'LineWidth', linewidth, 'MarkerSize', 2, 'Color', 'black');

hold off
set(gca, 'XScale', 'log', 'Layer', 'top', 'FontSize', tickFS, 'XTick', [1 1e2 1e4], 'XTickLabels', {'1', '100', '10000'})%, 'YScale', 'log')
xlabel(xlab, 'FontSize', axesFS, 'Color', 'none')
ylabel('Water loss (mL cm$^{-2}$ d$^{-1}$)', 'FontSize', axesFS)

axpos = get(gca, 'Position');
hold on
wax = gca;
plot([1 1]*rt(1), wax.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
plot([1 1]*rt(2), wax.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
plot([1 1]*rv, wax.YLim, ':', 'Color', cmap.gray, 'LineWidth', 2)
xlim([1 2*r*1e9])


box on
wax.Children = wax.Children([end 1:(end-1)]);

%% pressure
try
    figure(rpcfig)
    clf
catch
    rpcfig = figure();
end

rpcfig.Position(3:4) = figsize;

hold on
semilogx(Rp, Pc, 'LineWidth', linewidth, 'Color', 'black');

ylim([0 4.2])

ax3 = gca;
set(ax3, 'XScale', 'log', 'YTick', 0:4, 'FontSize', tickFS, ...
    'XTick', [1 1e2 1e4], 'XTickLabels', {'1', '100', '10000'})
plot([1 1]*rt(1), ax3.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
plot([1 1]*rt(2), ax3.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)


xlabel(xlab, 'FontSize', axesFS, 'Color', 'none')
ylabel('Chamber pressure (MPa)', 'FontSize', axesFS)
xlim([1 3e4])

ax3.Children = ax3.Children([end 1:(end-1)]);
box on

%% concentration/osmotic potential
try
    figure(rcdfig)
    clf
catch
    rcdfig = figure();
end

rcdfig.Position(3:4) = figsize;

hold on
semilogx(Rp, 1e-3*(c*cb - cg), 'LineWidth', linewidth, 'Color', 'black');

plot([1 1e5], [1 1]*(cstar)*1e-3, '--', 'Color', cmap.c, 'LineWidth', 2)

ax3 = gca;
set(ax3, 'XScale', 'log', 'FontSize', tickFS,...
    'XTick', [1 1e2 1e4], 'XTickLabels', {'1', '100', '10000'},...
    'YTick', 1:.2:2)
ax3.Position = [0.1300 0.1222 0.7750 0.8028];
colormap(ax3, cmapJ)
xlabel(xlab, 'FontSize', axesFS, 'Color', 'none')
ylabel('$c_c-c_g$ (M)', 'FontSize', axesFS)
ylim([1 2])
xlim([1 3e4])
box off



axcp = copyobj(ax3, gcf);
delete(axcp.Children)
delete(axcp.XLabel)
axcp.XTick = [];
axcp.YAxisLocation = 'right';
axcp.Position = ax3.Position;

set(axcp, 'Color', 'none')

axcp.YLim = (ax3.YLim*1e3 + 100)*RT*(1e-6);
set(axcp, 'YTick', round(axcp.YLim(1)):round(axcp.YLim(2)))
axcp.XLim = ax3.XLim;
% axcp.YTick = 2:5;

for i = 1:length(axcp.YTickLabels)
    if length(axcp.YTickLabels{i}) == 1
        axcp.YTickLabels{i} = strjoin({'$-', axcp.YTickLabels{i},'$'},'');
    else
        axcp.YTickLabels{i} = '';
    end
end
ylabel(axcp,'Chamber osmotic potential (MPa)', 'FontSize', axesFS, 'Rotation',270,'VerticalAlignment','bottom');
hold on
plot([1 1]*rt(1), ax3.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
plot([1 1]*rt(2), ax3.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
xlabel(xlab, 'FontSize', axesFS, 'Color', 'none')
ax3.Children = ax3.Children([end 1:(end-1)]);
boxax = axes;
set(boxax, 'Color', 'none', 'box', 'on', 'Position', ax3.Position,...
    'XScale', 'log', 'XTick', [], 'YTick', [])
xlim([1 1e5])
box on

%% Just total salt flux
try
    figure(saltfluxfig)
    clf
catch
    saltfluxfig = figure();
end
saltfluxfig.Position(3:4) = figsize;
semilogx(Rp, Jsalt*1e6/(pi*r^2)*8.64, 'LineWidth', linewidth', 'Color', 'black')
ylim([0 3.75])
set(gca, 'FontSize', tickFS, 'XTick', [1 1e2 1e4], 'XTickLabels', {'1', '100', '10000'},...
    'YTick', 0:4)
sax = gca;
hold on

plot([1 1]*rt(1), sax.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
plot([1 1]*rt(2), sax.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
xlabel(xlab, 'FontSize', axesFS)
ylabel('Net salt flux ($\mu$mol cm$^{-2}$ d$^{-1}$)', 'FontSize', axesFS)
sax.Children = sax.Children([end 1:(end-1)]);
xlim([1 3e4])

%% flux densities
try
    figure(velocityfig);
    clf;
catch
    velocityfig = figure();
end
velocityfig.Position(3:4) = figsize;


dif = semilogx(Rp, vel(JD), '--', 'LineWidth', linewidth, 'Color', cmap.diff);
hold on
adv = semilogx(Rp, vel(JP), 'LineWidth', linewidth, 'Color', 'black');

set(gca, 'FontSize', tickFS, 'XTick', [1 1e2 1e4], 'XTickLabels', {'1', '100', '10000'}, ...
    'YTick', -.1:.1:.4)
vax = gca;
ylim([-.1 .18])
xlim([1 3e4])

plot([1 1]*rt(1), vax.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
plot([1 1]*rt(2), vax.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)

plot([1 1]*rv, vax.YLim, ':', 'Color', cmap.gray, 'LineWidth', 2)
xlabel(xlab, 'FontSize', axesFS)
ylabel('Flux density (nmol $\mu$m$^{-2}$ s$^{-1}$)', 'FontSize', axesFS)

legend([adv dif], {'Advective', 'Diffusive'})
vax.Children = vax.Children([(end-1):end 1:(end-2)]);

%% Export Efficiency

try
    figure(efficiencyfig)
    clf
catch
    efficiencyfig = figure();
end

efficiencyfig.Position(3:4) = figsize;

exef = Jsalt./sa.Jtot;
semilogx(Rp, exef/vbar, 'Linewidth', linewidth, 'Color', 'black') % umol/mL

ylim([0 4.5])
xlim([1 3e4])


set(gca, 'FontSize', tickFS, 'XTick', [1 1e2 1e4], 'XTickLabels', {'1', '100', '10000'}, ...
    'YTick', 0:4)
eax = gca;
hold on
plot([1 1]*rt(1), eax.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
plot([1 1]*rt(2), eax.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
eax.Children = eax.Children([end 1:(end-1)]);
xlabel(xlab, 'FontSize', axesFS)
ylabel('Export efficiency ($\mu$mol mL$^{-1}$)', 'FontSize', axesFS)
box on


%% advection vs diffusion fig
try
    figure(advfig)
    clf
catch
    advfig = figure();
end

advfig.Position(3:4) = [700 480];

dif = semilogx(Rp, JD(:,1)*1e9, '--', 'LineWidth', linewidth, 'Color', cmap.diff);
hold on
adv = semilogx(Rp, JP(:,1)*1e9, 'LineWidth', linewidth, 'Color', 'black');

ax3 = gca;
set(gca, 'FontSize', tickFS, 'XTick', [1 1e2 1e4], 'XTickLabels', {'1', '100', '10000'})
plot([1 1]*rt(1), ax3.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)
plot([1 1]*rt(2), ax3.YLim, '--', 'Color', cmap.c, 'LineWidth', 2)

xlabel(xlab, 'FontSize', axesFS)
ylabel('Molar salt flux (nmol s$^{-1}$)', 'FontSize', axesFS)
xlim([1 2*r*1e9])

legend([adv dif], {'Advective', 'Diffusive'}, 'Location', 'southwest', 'FontSize', 0.8*tickFS)
box on