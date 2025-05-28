#!/usr/bin/env python

import optparse
import drawer
import time
import pylab

desc="""drawMS automatically converts the fringes seen in the visibilities to locuses in the sky.
As the fringe "produced" by each individual baseline are rotating on the sky, each source modulate the visibility, depending 
on its distance to the phase center. This code does an FFT of the visibility of a each given baseline in a particular timeslot, along the time axis,
and derives a line on the sky. It reflects all the possible places where the source producing the given detected modulation could be.
The lines are then gridded onto an image. This is preliminary version, that do not deal with data chunks yet, so that it first reads the whole MS
and puts the visibilities into memory. Therefore this software should be used on averaged datasets containing a few channels. This is an experimental code,
so there no guarantee there are not bugs!!! Questions and suggestions: cyril.tasse@obspm.fr"""

opt = optparse.OptionParser(usage='Usage: %prog --ms=somename.MS <options>',version='%prog version 1.0',description=desc)

group = optparse.OptionGroup(opt, "* Necessary options", "Won't work if not specified.")
group.add_option('--ms',help='Input MS to draw [no default]',default='')
opt.add_option_group(group)

group = optparse.OptionGroup(opt, "* Data selection options", "ColName is set to DATA column by default, and other parameters select all the data.")
group.add_option('--ColName',help='Name of the column to work on. Default is %default. For example: --ColName=CORRECTED_DATA',default="DATA")
group.add_option('--uvrange',help='UV range (in meters, not in lambda!). Default is %default. For example: --uvrange=100,1000',default="0,10000000")
group.add_option('--wmax',help='Maximum W distance. Default is %default. ',default=10000000)
group.add_option('--timerange',help='Time selection range, in fraction of total observing time. For example, --timerange=0.1,0.2 will select the second 10% of the observation. Default is %default. ',default="0,1")
group.add_option('--AntList',help='List of antennas to compute the lines for. Default is all. For example: --AntList=0,1,2 will plot 0-n, 1-n, 2-n',default="")
group.add_option('--FillFactor',help='The probability of a baseline/timeslot to be processed. Default is %default. Useful when large dataset are to be drawn. For example --FillFactor=0.1 will result in a random selection of 10% of the data',default=1.)
opt.add_option_group(group)

group = optparse.OptionGroup(opt, "* Algorithm options", "Default values should give reasonable results, but all of them have noticeable influence on the results")
group.add_option('--timestep',help='Time step between the different time chunks of which the drawer does the fft. Default is %default. ',default=500)
group.add_option('--timewindow',help='Time interval width centered on the time bin controlled by --timestep. If not defined then it is set to --timestep.',default=0)
group.add_option('--snrcut',help='Cut above which the fringe is drawn. Default is %default. ',default=5.)
group.add_option('--maskfreq',help='When a fringe is found, it will set the fft to zero in that 1D pixel range. Default is %default. ',default=2.)
group.add_option('--MaxNPeaks',help='Maximum number of fringes it will find per baseline and timeslot. Default is %default. ',default=7)
group.add_option('--NTheta',help='Number of angles in the l-m plane the algorithm will solve for. Default is %default. ',default=20)
opt.add_option_group(group)

group = optparse.OptionGroup(opt, "* Fancy options", "Plot NVSS sources, or make a movies.")
group.add_option('--RadNVSS',help='Over-plot NVSS sources within this radius. Default is %default (in beam diameter). ',default=0)
group.add_option('--SlimNVSS',help='If --RadNVSS>0, plot the sources above this flux density. Default is %default Jy. ',default=0.5)
group.add_option('--MovieName',help='Name of the directory that contains the movie (.mpg), the individual timeslots (.png), and the stack (.stack.png). Each page correspond to the data selected by --timewindow, separated by --timestep. For example --MovieName=test will create a directory "dMSprods.test". Default is None. ',default="")
opt.add_option_group(group)
options, arguments = opt.parse_args()

# if options.ms=="":
#     print("Give an MS name!")
#     exit()
    
drawer.print_logo()
drawer.init(options.ms,True,ColName=options.ColName)
drawer.set_options(options)
drawer.display_time("",init=True)
drawer.find_source(True,noplot=True)
drawer.display_time("Line finding")
if options.MovieName=="":
    drawer.plot_gridded()
    drawer.display_time("Line drawing")
    pylab.show()

if options.MovieName!="":
    drawer.movie_maker(options.MovieName)
    drawer.display_time("Movie making")
