import os
from anuga.utilities import plot_utils as util
from matplotlib import pyplot as pyplot
import numpy
import numpy as np
from anuga import get_flow_through_cross_section
from anuga.utilities import plot_utils
from anuga.geospatial_data.geospatial_data import Geospatial_data

verbose= True

#os.system('python /home/wollongong/anuga_core/anuga/utilities/sww_merge.py -f 2007 -np 4')

swwfile = '2007.sww'
## Hydrograph plot

xll = 382250.0
yll = 6354265.0
cross_section = [[130.+xll,80.+yll],[103+xll, 100.+yll]]

time, Q = get_flow_through_cross_section(swwfile,
                                         cross_section,
                                         verbose=True)
    
f = open('ANUGA_Direct_Rain_Hydrograph.csv', 'w')
for t in range(len(time)):
    f.write('%s %s\n' % (str(time[t]), str(Q[t])))
f.close()

for f in ('Observations/wbnm_hydrograph_burst_WRL.csv', 'ANUGA_Direct_Rain_Hydrograph.csv', 'Observations/wbnm_hydrograph_burst_dpm.csv'):
    data=np.loadtxt(f)
    X=data[:,0]/3600.
    Y=data[:,1]
    pyplot.plot(X,Y)
    pyplot.title("Flow through transect 1\nWBNM WRL (blue) -vs- ANUGA (orange) -vs- WBNM DPM (green)")
    pyplot.xlabel('time (hours)')
    pyplot.ylabel('flow (m3/s)')
    pyplot.savefig('flowrate.png',dpi=100,bbox_inches='tight')

p=util.get_output(swwfile)
p2=util.get_centroids(p)
# Time index at last time
tindex = len(p2.time)-1

if verbose: print 'calculating experimental transect'

x_data =    [ 0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0]
#vel =   [ 0.0, 0.0, 1.1, 3.2,  3.4, 2.4,  3.2,  3.2,  3.7,  3.1,  0.4,  0.0]
vel_data =   [ 0.0, 0.4, 3.1, 3.7,  3.2, 3.2,  2.4,  3.4,  3.2,  1.1,  0.0,  0.0]
#depth = [ 0.0, 0.0, 0.1, 0.5,  0.45, 0.4, 0.55, 0.1, 0.1,  0.05,  0.04, 0.0]
depth_data = [ 0.0, 0.04, 0.05, 0.1,  0.1, 0.55, 0.4, 0.45, 0.5,  0.1,  0.0, 0.0]

from scipy import interpolate

fvel = interpolate.interp1d(x_data, vel_data)
fdepth = interpolate.interp1d(x_data, depth_data)

if verbose: print 'calculating model heights at observation points'
# Get nearest wet points to 'point observations'
point_observations = numpy.genfromtxt(
    'Observations/ObservationPoints.csv',
    delimiter=",",skip_header=1)

nearest_points = (point_observations[:,0]*0. - 1).astype(int)

for i in range(len(nearest_points)):
    # Compute distance of ANUGA points to observation, and
    # if the ANUGA point is dry then add a large value
    # Then find index of minimum
    n = ( (p2.x+p2.xllcorner-point_observations[i,0])**2 + \
          (p2.y+p2.yllcorner-point_observations[i,1])**2 + \
          (p2.stage[tindex,:] <= p2.elev)*1.0e+06).argmin()
    nearest_points[i] = n

#import pdb; pdb.set_trace()

if verbose: print 'Plot transect'
## Plot transect 1 [need to guess appropriate end points as these are not so
## clear from the report]

x_adjust = 1240
y_adjust = 760
xx=util.near_transect(p2,[x_adjust + 106., y_adjust + 96.], [x_adjust + 132.,y_adjust + 78.],tol=0.5)
xx0 = numpy.array(xx[0])
xx1 = numpy.array(xx[1])

for tindex, time in enumerate(p2.time):
    # time at which event peak occurs for the burst rainfall run
    # if running the entire storm, time at peak is 139200 instead of 17700
    if time == 17700: 
        
        #f = open('Stage_point_comparison_%s.csv'% time,'w')
        f = open('Stage_point_comparison.csv','w')
        f.writelines( 'Field, ANUGA, TUFLOW, ANUGA minus Field, ANUGA minus TUFLOW \n' )
        for i in range(len(nearest_points)):
            po = point_observations[i,-2]
            tu = point_observations[i,-1]
            anuga_data = p2.stage[tindex, nearest_points[i]]
            newline = str(round(po,2)) + ', ' + str(round(anuga_data,2)) + ', ' + str(tu) + ', ' + \
                  str(round(anuga_data - po,2)) + ', ' + str(round(anuga_data - tu,2)) + '\n'
            f.writelines(newline)
        
        f.flush()
        f.close()
            
        if verbose: print 'Transect Time %s sec.'% time
        pyplot.clf()
        pyplot.figure(figsize=(16,10.5))
        pyplot.subplot(121)
        pyplot.scatter(p2.x, p2.y, c=p2.elev,edgecolors='none')
        # Add nice elevation data
        colVals = numpy.maximum(numpy.minimum(p2.elev, 25.), 19.)
        util.plot_triangles(p, values = colVals, edgecolors='none')
    
        pyplot.gca().set_aspect('equal')
        pyplot.scatter(p2.x[xx0],p2.y[xx0],color='green')
        pyplot.xlim( (x_adjust + 40., x_adjust + 160.))
        pyplot.ylim( (y_adjust + 0.,y_adjust + 140.))
        pyplot.title('Transect points in green')
    
        #import pdb; pdb.set_trace()
        pyplot.subplot(222)
        pyplot.scatter(xx1,p2.vel[tindex,xx0],color='green',label='model')
        pyplot.scatter(xx1,fvel(xx1),color='blue',label='data')
        pyplot.legend(loc='upper left')
        #pyplot.xlim(0,25)
        pyplot.title('Final flow speed along the transect at time %g'% time)
    
        pyplot.subplot(224)
        pyplot.scatter(xx1,p2.stage[tindex,xx0]-p2.elev[xx0],color='green',label='model')
        pyplot.scatter(xx1,fdepth(xx1),color='blue',label='data')
        pyplot.legend(loc='upper left')
        #pyplot.xlim(0,25)
        pyplot.title('Depth along the transect at time %g'% time)
    
        pyplot.savefig('Transect1.png', bbox_inches='tight')
    
        pyplot.close()
    
        if verbose: print 'Plot velocity field Time %g sec.'% time
        pyplot.clf()
        
        # Velocity vector plot
        pyplot.figure(figsize=(16,22))
        pyplot.scatter(p2.x,p2.y,c=(p2.elev>24.),edgecolors='none', s=0.2)
        pyplot.gca().set_aspect('equal')
        pyplot.xlim((1300.,1550.))
        pyplot.ylim(( 750.,1150.))
        #k=range(0,len(p2.x),2) # Thin out the vectors for easier viewing
        colVals = numpy.maximum(numpy.minimum(p2.elev, 25.), 19.)
        util.plot_triangles(p, values = colVals, edgecolors='white')
        #k = range(len(p2.x))
        # Thin out the triangles
        k = (((10.*(p2.x - p2.x.round())).round()%2 == 0.0)*((10.*(p2.y - p2.y.round())).round()%2 == 0.0)).nonzero()[0]
        pyplot.quiver(p2.x[k],p2.y[k],p2.xvel[tindex,k], p2.yvel[tindex,k],
                      scale_units='xy',units='xy',width=0.1,
                      color='black',scale=1.0)
        pyplot.savefig('velocity_stationary.png', dpi=100, bbox_inches='tight')
        
        
        
        ## Froude number plot
        if verbose: print 'Plot Froude number plot'
        pyplot.clf()
        pyplot.figure(figsize=(6,8))
        froude_number = p2.vel[tindex]/(numpy.maximum(p2.height[tindex], 1.0e-03)*9.8)**0.5
        froude_category = (froude_number>1.).astype(float) + (froude_number > 0.).astype(float)
        pyplot.scatter(p2.x,p2.y,edgecolors='none', s=0.2)
        
        ## Fake additions to plot to hack matplotlib legend
        pyplot.scatter(0.,0., color='FireBrick',label='>1', marker='s')
        pyplot.scatter(0.,0., color='PaleGreen',label='0-1', marker='s')
        pyplot.scatter(0.,0., color='blue',label='0',marker='s')
        
        pyplot.gca().set_aspect('equal')
        util.plot_triangles(p, values = froude_category, edgecolors='none')
        pyplot.xlim((1300.,1550.))
        pyplot.ylim(( 750.,1150.))
        pyplot.title("Froude Number zones: 0, (0,1], or >1")
        pyplot.legend(loc='upper left')
        pyplot.savefig('froudeNumber.png', dpi=100, bbox_inches='tight')
"""
#extract maximums from sww as GEOTIFFS
output_quantities=['depthIntegratedVelocity','depth', 'velocity', 'stage']

for i in output_quantities:
    print 'extracting ', [i]
    plot_utils.Make_Geotif(swwFile='2007.sww', output_quantities=[i], myTimeStep='max', CellSize=1.0, velocity_extrapolation=True, min_allowed_height=0.01, EPSG_CODE=28356, verbose=False, k_nearest_neighbours=3)
"""
