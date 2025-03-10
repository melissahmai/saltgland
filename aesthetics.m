set(groot, 'defaultFigureColor', 'white', ...
'defaultTextInterpreter', 'latex', ...
'defaultLegendInterpreter', 'latex', ...
'defaultAxesTickLabelInterpreter', 'latex', ...
'defaultTextFontSize', 24, ...
'defaultAxesLabelFontSizeMultiplier', 1.75, ...
'defaultAxesLineWidth', 1.5, ...
'defaultAxesFontSize', 12)

axesFS = 21;
tickFS = 15;
legFS = 12;
textFS = 18;
cmapJ = [.9 .6 0; 0 .3 .9];
figsize = [700 540];
linewidth = 3;
cmap = struct('psi', [.2 .5 .1], ...
    'p', [.6 0 .5], 'c', [.9 .6 0], ...
    'j', [0 .2 .8], 'gray', [1 1 1]*.5,...
    'diff', [.8 .2 0]);
