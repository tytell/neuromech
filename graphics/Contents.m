% Functions for graphics.
%
% addscalebars - Add a scale bar (or bars) to an axis.  Useful for plots in
%   which the scale of the axes is meaningful, but not the absolute
%   numbers.
% buildAVI - Builds an AVI movie out of a 3D image matrix.
% camcopy - Copies the camera settings from one figure to another.
% circle - Draws a circle.
% figblackbg - Sets a black background for a figure, changing black lines into white.
% figwhitebg - Sets a white background for a figure, changing white lines
%   into black.
% gplotmatreg - Plots multivariate data in a matrix of small axes.  Similar
%   to gplotmatrix in the Statistics toolbox, but more powerful.
% horizplot - Plots a horizontal line or lines.
% imshow6 - Shows an image using the syntax of Matlab version 6.
% latextable - Builds a table in LaTeX, including sparklines as suggested
%   by Edward Tufte.
% makePresentation - Formats a figure for a presentation by increasing font
%   sizes and line widths.
% makePublication - Formats a figure for publication.  Attempts to size the
%   figure correctly for different journal specifications.
% movieplot - Shows 3D data as a movie.  Can also export the movie as an
%   AVI.
% mplot - Slight update to plot, allowing different linespecs for multiple
%   lines to be specified easily.  Also includes an option for filled
%   symbols.
% multiplot - Overlays many different types of plots and allows the user to
%   step through multiple frames of data in an easy way.  Useful for
%   complex time series data.
% orthoplot - Plots volumetric data with orthogonal slices.  Similar to
%   functions provided on confocal microscopes.
% overlayplot - Provides a (relatively) simple mechanism for overlaying
%   plots with different y (or x) axes.  Related to plotyy, but much more
%   powerful.
% plotgroups - Fancy function for plotting data that are grouped.  Provides
%   options for multiple grouping axes, calculating means and error bars,
%   and many different coloring and labeling options.
% plotpolybounds - Plots error bounds on linear regressions.
% plottimeseries - Plots long timeseries data as multiple rows.  Attempts
%   to format the figure so that it will print nicely.
% printallfigs - Prints all open figures (or a selection of figures)
% saveallfigs - Saves all open figures (or a selection of them)
% showtable - Shows a table in a figure window, which may include sparklines 
%   (small plots).  Has lots of formatting options.
% splot - Useful for dealing with 3D datasets.  Runs squeeze on each of its
%   arguments bofer plotting.
% superboxplot - Like boxplot, but takes an x coordinate, and has some
%   other convenient options
% vertplot - Plots a vertical line or lines
% whitejet - Another colormap, like jet, but with white in the middle
%   rather than green.
% xtick - Adjusts the tick values and/or labels on the x axis.
% ytick - Adjusts the tick values and/or labels on the y axis.