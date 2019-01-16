User manual for hrd_plot and hrd_plot_old routines
********************************

Made by Oskari Tervo, summer 2017, University of Oulu
Questions about the code can be sent to tervo.oskari@gmail.com. I'll try to answer as best as I can.
I've written (perhaps an overabundance of) comments in the code which should help to better understand the routine, if one wishes to.





This routine plots particle data based on the detections made by Cassini's HRD-instrument during Enceladus flybys.


========================================

Files needed to run the routine are:



bindat.pro
flybys.pro
gentab.pro
modulus.pro
HRD-counter datafile, e.g. hrd_2010_350_365_prc.tab
"Kernels" folder, which includes the necessary CK, SPK etc. kernels for the timeframe the user is interested in.



HRD-data can be downloaded from https://sbn.psi.edu/archive/cohrd/
and kernels from https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/

Also included is hrd_plot_old routine, which is required for flybys before year 2007. For this routine all of the above is required, but also 
a separate pointing datafile must exist, e.g. cdapointing2005v2.tab, because pointing data before year 2007 is unobtainable through SPICE. 
These files can be downloaded with HRD-datafiles through the link above. Please also note that hrd_plot_old routine is considerably slower
to process than hrd_plot, so do not be alarmed if nothing happens immediately after running it.

I've also included the necessary datafiles from flybys done at 2005-07-14 and 2010-12-21 for convenience and testing purposes.

========================================

HOW TO USE:
************

At the bottom of the hrd_plot.pro file arguments in the draw_graph function can be changed to best suit user's needs:

1st argument is the HRD-datafile of the flyby date one is interested in.
2nd argument is the width of the plot timeframe in minutes. For example, a value of "10" will plot data from 10 minutes before AND after the closest approach.
3rd argument is the UTC date of the flyby in format 'YEAR-MONTH-DAY'.
4th argument defines how many plots of different sized particles are shown after the routine has ended:
	Value of "1" plots data for M1 and m1 particles
	Value of "2" plots data for M1, M2, m1, and m2 particles
	Value of "3" plots M1-M3 and m1-m3
	Value of "4" plots M1-M4 and m1-m4
5th argument defines what data the routine plots:
	Value "f" plots particles as a function of time
	Value "d" plots differentiated time-particle-function
	Value "n" plots differentiated data where HRD's pointing is taken to account, i.e. number density

For hrd_plot_old routine the 5th argument is optional. If given, it must be the above mentioned pointing data file, in which case the data plotted is
the number density. Otherwise differentiated function is plotted similar to option "d" mention above.


KNOWN BUGS:

-When pointing data is applied for number density plots, the graphs are flipped upside down. Can be solved with abs() function, although far from ideal.
-Plot histograms can act out weirdly. Most likely something to do with the previous bug.







