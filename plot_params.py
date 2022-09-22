import matplotlib as mpl
import matplotlib.font_manager as font_manager

# Add every font at the specified location
# arial ttf dir containing location
font_dir = ['']
for font in font_manager.findSystemFonts(font_dir):
    font_manager.fontManager.addfont(font)

# plot style
mpl.rcParams['backend'] = 'Agg'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['font.size'] = 8

# Set font family globally
mpl.rcParams['font.family'] = 'Arial'

COLOR_CYCLER = [(0.00, 0.00, 0.00), (0.00, 0.00, 0.40), (0.84, 1.00, 0.00),
                (1.00, 0.00, 0.34), (0.62, 0.00, 0.56), (0.05, 0.30, 0.63),
                (1.00, 0.90, 0.01), (0.00, 0.37, 0.22), (0.00, 1.00, 0.00),
                (0.58, 0.00, 0.23), (1.00, 0.58, 0.49), (0.64, 0.14, 0.00),
                (0.00, 0.08, 0.27), (0.57, 0.82, 0.80), (0.38, 0.05, 0.00),
                (0.42, 0.41, 0.51), (0.00, 0.00, 1.00), (0.00, 0.49, 0.71),
                (0.42, 0.51, 0.42), (0.00, 0.68, 0.49), (0.76, 0.55, 0.62),
                (0.75, 0.60, 0.44), (0.00, 0.56, 0.61), (0.37, 0.68, 0.31),
                (1.00, 0.00, 0.00), (1.00, 0.00, 0.96), (1.00, 0.01, 0.62),
                (0.41, 0.24, 0.23), (1.00, 0.45, 0.64), (0.59, 0.54, 0.91),
                (0.60, 1.00, 0.32), (0.65, 0.34, 0.25), (0.00, 1.00, 1.00),
                (1.00, 0.93, 0.91), (1.00, 0.54, 0.00), (0.74, 0.78, 1.00),
                (0.00, 0.82, 1.00), (0.73, 0.53, 0.00), (0.46, 0.27, 0.69),
                (0.65, 1.00, 0.82), (1.00, 0.65, 1.00), (0.47, 0.30, 0.00),
                (0.48, 0.28, 0.51), (0.15, 0.20, 0.00), (0.00, 0.28, 0.33),
                (0.26, 0.00, 0.17), (0.71, 0.00, 1.00), (1.00, 0.69, 0.40),
                (1.00, 0.86, 0.40), (0.56, 0.98, 0.57), (0.49, 0.18, 0.82),
                (0.74, 0.83, 0.58), (0.90, 0.44, 1.00), (0.87, 1.00, 0.45),
                (0.00, 1.00, 0.47), (0.00, 0.61, 1.00), (0.00, 0.39, 0.00),
                (0.00, 0.46, 1.00), (0.52, 0.66, 0.00), (0.00, 0.73, 0.09),
                (0.47, 0.51, 0.19), (0.00, 1.00, 0.78), (1.00, 0.43, 0.25),
                (0.91, 0.37, 0.75)]
