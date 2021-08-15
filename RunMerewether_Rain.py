"""
Prepared by Petar Milevski and Stephen Roberts

Pasha Bulka Storm Event 2007 - storm burst model

with option to run the entire 2 day event if you are not time limited
"""

import anuga, numpy, time, os, glob, shutil
from anuga import Rate_operator
from anuga import file_function, Polygon_function, read_polygon, create_mesh_from_regions, Domain
from anuga import distribute, myid, numprocs, finalize, barrier

basename = 'Model/DEM/M'
outname = '2007'
meshname = 'Model/DEM/M.msh'

if anuga.myid == 0:
    
    # First untar the large DEM file
    shutil.unpack_archive(os.path.join('Model','DEM','M.csv.tgz'), os.path.join('Model','DEM'))

    CatchmentDictionary = {os.path.join('Model','Bdy','catchment.csv'):25., os.path.join('Model','Bdy','fine.csv'):1.}
    
    bounding_polygon = anuga.read_polygon(os.path.join('Model','Bdy','catchment.csv'))
    interior_regions = anuga.read_polygon_dir(CatchmentDictionary, os.path.join('Model','Bdy'))
    
    BuildingDictionary = 'Model/Houses/houses.csv'
    interior_holes = anuga.read_multi_poly_file(BuildingDictionary)
    
    create_mesh_from_regions(bounding_polygon,
        boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
        maximum_triangle_area=100.,
        interior_regions=interior_regions,
        interior_holes=interior_holes,
        filename=meshname,
        verbose=True)
 
    domain = Domain(meshname, use_cache=True, verbose=True)
    print ('reading terrain grid')
    domain.set_quantity('elevation', filename=basename+'.csv', use_cache=True, verbose=True, alpha=0.99)
    domain.set_name(outname)
    print (domain.statistics())

    print ('Applying mannings roughness')
    Roughness_Polygons_directory = 'Model/Mannings'
    friction_list = anuga.get_polygon_value_list(Roughness_Polygons_directory)
    domain.set_quantity('friction', anuga.Polygon_function(friction_list, default=0.02, geo_reference=domain.geo_reference))
    domain.set_quantity('stage', 0.)

else:
	domain = None

verbose = True
if anuga.myid == 0 and verbose: print ('DISTRIBUTING DOMAIN')
domain = anuga.parallel.distribute(domain)
domain.set_minimum_storable_height(0.05)
if anuga.myid == 0 and verbose: print ('CREATING INLETS')  

#change from "rain_burst.tms" to "rain_entire_storm.tms" to run the 2 day event
polygon = anuga.read_polygon('Model/Rain/Gauge/catchment.csv')
rainfall = anuga.file_function('Model/Rain/Rainfall/rain_burst.tms', quantities='rate')
op1 = Rate_operator(domain, rate=rainfall, factor=1.0e-3, polygon=polygon, default_rate = 0.0) 

print ('Available boundary tags', domain.get_boundary_tags())
Br = anuga.Reflective_boundary(domain)
Bd = anuga.Dirichlet_boundary([0,0,0])
domain.set_boundary({'interior': Br, 'exterior': Bd, 'west': Bd, 'south': Bd, 'north': Bd, 'east': Bd})

if anuga.parallel.myid == 0 and verbose: print ('EVOLVE')
t0 = time.time()

#change final time from 46500. to 177600. to run the 2 day event
for t in domain.evolve(yieldstep = 300., finaltime = 46500.):  
    if anuga.parallel.myid == 0:
        domain.write_time()

if anuga.parallel.myid == 0:
    print ('Number of processors %g ' %anuga.parallel.numprocs)
    print ('That took %.2f seconds' %(time.time()-t0))
    print ('Communication time %.2f seconds'%domain.communication_time)
    print ('Reduction Communication time %.2f seconds'%domain.communication_reduce_time)
    print ('Broadcast time %.2f seconds'%domain.communication_broadcast_time)
    
domain.sww_merge(delete_old=True)
anuga.finalize()
