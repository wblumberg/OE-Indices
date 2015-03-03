import sys
sys.path.append('/home/greg.blumberg/python_pkgs/')
from pylab import *
from netCDF4 import Dataset, MFDataset
import indices_helper as sti
import sys
from datetime import datetime
import sharppy
print sys.argv[1]
#d = MFDataset(sys.argv[1])


l = 2
d = Dataset(sys.argv[1])
height = d.variables['height'][:]
pres = d.variables['pressure'][:l]
bt = d.variables['base_time'][:]
to = d.variables['time_offset'][:l]
Xop = d.variables['Xop'][:l]
Sop = d.variables['Sop'][:l]
converged_flag = d.variables['converged_flag'][:l]
height = height * 1000.

num_perts = 700 # The ideal number of profiles to compute the index percentiles
cush = 600 # The number of profiles where indices MUST be computed from (allows for crashing)
ideb, details = sti.makeIndicesErrors(Xop, Sop, height, pres, num_perts, cush, converged_flag)

fn = sys.argv[1]

#spl2 = '.'.join(fn.split('/')[-1].split('.'))
name = fn.replace('aerioe1turn', 'aerioeidx1blum')
name = name.replace('c1', 'c2')
#sp =  '/'.join(fn.split('/')[:-1])
#std = sp + '/' + name  

out = Dataset(name, 'w', format='NETCDF3_CLASSIC')

out.Date_created = datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S")
out.description = "This file contains convective indices generated from thermodynamic profiles produced by an OE (i.e. AERIoe/MWRoe) algorithm.  The convective indicies within this file also contain estimates of uncertainity.  These were produced by Monte Carlo sampling the retrieved profile via its covariance matrix."
out.comment = "Uses the SHARPpy libraries to calculate the convective indices."
out.author = "Greg Blumberg; OU/CIMMS"
out.input_file = fn
out.num_perts = str(num_perts)
out.cushon = str(cush)
out.SHARPpy_version = sharppy.__version__

out.createDimension('time', len(ideb))
out.createDimension('percentiles', 3)
out.createDimension('single', 1)

var = out.createVariable('base_time', 'f4', ('single',))
var[:] = bt
var.long_name = d.variables['base_time'].long_name
var.units = str(d.variables['base_time'].units)

var = out.createVariable('time_offset', 'f4', ('time', ))
var[:] = to
var.long_name = d.variables['time_offset'].long_name
var.units = str(d.variables['time_offset'].units)

var = out.createVariable('percentiles', 'f4', ('percentiles', ))
var[:] = [25,50,75]
var.long_name = "Percentiles for the errorbars"
var.units = "%"

# loop through the indices and add them to the netCDF file
for var_name in np.sort(details.keys()):
	longname = details[var_name][0]
	unit = details[var_name][1]
	var = out.createVariable(var_name, 'f4', ('time', 'percentiles',))
	var.long_name = longname
	var.units = unit
	data = np.empty((len(ideb),3))
	for i in xrange(len(ideb)):
		data[i,:] = ideb[i][var_name]
	var[:] = data

var = out.createVariable('hatchOpen', 'i4', ('time', ))
var[:] = d.variables['hatchOpen'][:l]
var.long_name = d.variables['hatchOpen'].long_name
var.units = str(d.variables['hatchOpen'].units)

var = out.createVariable('converged_flag', 'i4', ('time', ))
var[:] = d.variables['converged_flag'][:l]
var.long_name = d.variables['converged_flag'].long_name
var.units = str(d.variables['converged_flag'].units)

var = out.createVariable('rms', 'f4', ('time', ))
var[:] = d.variables['rms'][:l]
var.long_name = d.variables['rms'].long_name
var.units = str(d.variables['rms'].units)

var = out.createVariable('lwp', 'f4', ('time', ))
var[:] = d.variables['lwp'][:l]
var.long_name = d.variables['lwp'].long_name
var.units = str(d.variables['lwp'].units)

var = out.createVariable('qc_flag', 'i4', ('time', ))
var[:] = d.variables['qc_flag'][:l]
var.long_name = d.variables['qc_flag'].long_name
var.units = str(d.variables['qc_flag'].units)


out.close()

