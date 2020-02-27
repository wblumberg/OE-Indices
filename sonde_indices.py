print 'python sonde_indices.py [-a or -w] <input_file>'

import sys
sys.path.append('/home/greg.blumberg/python_pkgs/')
from pylab import *
from netCDF4 import Dataset, MFDataset
import indices_helper as sti
import sys
from datetime import datetime
import sharppy
import warnings
warnings.filterwarnings("ignore")

'''
    sonde_indices.py

    AUTHOR: Greg Blumberg (OU/CIMMS)

    This script takes in ARM formatted radiosonde files and performs 
        a.) Monte Carlo sampling of the retrievals
        to generate approximately 500 retrievals
        b.) uses the SHARPpy libraries to perform convection 
        index calculations on all 500 retrievals.
    
    After this the convection indices are sorted, the 25th, 50th, and 75th percentiles
    are pulled out for each index and then saved in a netCDF file.

    Running this script requires the netCDF4-python libraries, the SHARPpy libraries,
    and IPython to perform the calculations in an efficient manner.

    To run this script, use this command:
        python sonde_indices.py [-a or -w] <input_file>

    where <input_file> is the path to the OE netCDF file and the flag -a means append to an existing file
    and the -w means overwrite any existing file.

    Before creating a new OE indices netCDF file, this program first checks for an existing file
    because it may need to be appended to (for example if this script is running in real-time.

'''

temp_sigma = 0.1 # C
rh_sigma = 3 # %
 
def rh2q(temp, pres, rh):
    """
    Inputs
    ------
    temp [K]
    pres [Pa]
    rh [fraction]
    """
    Rv = 461.
    L = 2.453 * 10**6
    es = 6.11 * np.exp((L/Rv)*((1./(273.15)) - (1./temp)))
    e = rh * es
    q = (0.622*e) / ((pres/100.) - e)
    return q

def perturb_radiosonde(tdry, rh, pres, num_perts):
    new_t = np.empty((num_perts,len(tdry)))
    new_q = np.empty((num_perts,len(tdry)))  
    for i in range(num_perts):                                                                                                 
        new_t[i] = temp_sigma * np.random.normal(0, 1, len(tdry)) + tdry  
        new_rh = rh_sigma * np.random.normal(0, 1, len(tdry)) + rh
        idx = np.where(new_rh < 0)[0]
        new_rh[idx] = 0
        idx = np.where(new_rh > 100)[0]
        new_rh[idx] = 100 
        new_q[i] = rh2q(tdry+273.15, pres*100., new_rh/100.)
    return new_t, 1000.*new_q

def run(fn, out):
    #flag = sys.argv[1]
    fn = sys.argv[1]
    out = sys.argv[2]
    # Replace this for the radiosonde data formats
    name = fn.split('/')[-1].replace('sondered', 'sondemcidx1blum')
    name = name.replace('c1', 'c2')
    name = out + '/' + name
    stride = 1
    d = Dataset(fn)
    bt = d.variables['base_time'][:]
    to = d.variables['time_offset'][::stride]
    height = d.variables['alt'][::stride]
    pres = np.asarray([d.variables['pres'][::stride]])
    rh = d.variables['rh'][::stride]
    tdry = d.variables['tdry'][::stride]
    print tdry, pres, rh
    q = rh2q(tdry+273.15, 100.*pres[0], rh/100.)#d.variables['q'][:]
    print d.variables['q'][:]
    print q
    stop
    print "QC check of profile."
    print np.max(height), np.min(height)
    print np.max(pres[0]), np.min(pres[0])
    print np.max(tdry), np.min(tdry)
    print np.max(rh), np.min(rh)
    print np.max(q), np.min(q)
    stop
    print height.shape, pres.shape, tdry.shape, q.shape
    t_dist, q_dist = perturb_radiosonde(tdry, rh, pres, 1000)
    Sop = np.asarray([np.cov(np.hstack((t_dist, q_dist)).T)])
    Xop = np.asarray([np.concatenate((tdry,q))])
    print Xop.shape, Sop.shape
    print np.max(Sop), np.min(Sop)
    print np.sqrt(np.diag(Sop[0]))
    #stop 
    #converged_flag = d.variables['converged_flag'][beg_idx:end_idx]
    converged_flag = [1]
    height = height# * 1000.

    num_perts = 2 # The ideal number of profiles to compute the index percentiles
    cush = 1 # The number of profiles where indices MUST be computed from (allows for crashing)
    ideb, details = sti.makeIndicesErrors(Xop, Sop, height, pres, num_perts, cush, converged_flag)

    #if flag == '-w' and len(glob.glob(name)) != 0:
    #    os.system('rm ' + name)
    #else:
    #    print "Appending to the file."

    out = Dataset(name, 'w', format='NETCDF3_CLASSIC')

    out.Date_created = datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S")
    out.description = "This file contains convective indices generated from thermodynamic profiles produced by an OE (i.e. AERIoe/MWRoe) algorithm.  The convective indicies within this file also contain estimates of uncertainity.  These were produced by Monte Carlo sampling the retrieved profile via its covariance matrix."
    out.comment = "Uses the SHARPpy libraries to calculate the convective indices."
    out.author = "Greg Blumberg; OU/CIMMS"
    out.email = 'wblumberg@ou.edu'
    out.input_file = fn
    out.num_perts = str(num_perts)
    out.cushon = str(cush)
    out.SHARPpy_version = sharppy.__version__

    out.createDimension('time', len(ideb))
    out.createDimension('samples', num_perts)
    out.createDimension('single', 1)

    var = out.createVariable('base_time', 'f4', ('single',))
    var[:] = bt
    var.long_name = 'beginning time'#d.variables['base_time'].long_name
    var.units = str(d.variables['base_time'].units)

    var = out.createVariable('time_offset', 'f4', ('time', ))
    var[:] = [0]
    var.long_name = 'time offset since base_time'
    var.units = str(d.variables['time_offset'].units)

    #var = out.createVariable('percentiles', 'f4', ('percentiles', ))
    #var[:] = [25,50,75]
    #var.long_name = "Percentiles for the errorbars"
    #var.units = "%"

    # loop through the indices and add them to the netCDF file
    for var_name in np.sort(details.keys()):
        longname = details[var_name][0]
        unit = details[var_name][1]
        var = out.createVariable(var_name, 'f4', ('time', 'samples',))
        var.long_name = longname
        var.units = unit
        data = np.empty((len(ideb),num_perts))
        print longname
        for i in xrange(len(ideb)):
            data_subset = ideb[i][var_name]
            #print len(data_subset)
            print num_perts, len(data_subset)
            if len(data_subset) > num_perts:
                data_subset = data_subset[:num_perts]
            data_subset = np.concatenate((data_subset, np.ones(num_perts - len(data_subset))*-9999))
            #print data_subset
            data[i,:] = data_subset
        var[:] = data
    """
    var = out.createVariable('hatchOpen', 'i4', ('time', ))
    var[:] = d.variables['hatchOpen'][beg_idx:end_idx]
    var.long_name = d.variables['hatchOpen'].long_name
    var.units = str(d.variables['hatchOpen'].units)

    var = out.createVariable('converged_flag', 'i4', ('time', ))
    var[:] = d.variables['converged_flag'][beg_idx:end_idx]
    var.long_name = d.variables['converged_flag'].long_name
    var.units = str(d.variables['converged_flag'].units)

    var = out.createVariable('rmsr', 'f4', ('time', ))
    var[:] = d.variables['rmsr'][beg_idx:end_idx]
    var.long_name = d.variables['rmsr'].long_name
    var.units = str(d.variables['rmsr'].units)

    var = out.createVariable('rmsa', 'f4', ('time', ))
    var[:] = d.variables['rmsa'][beg_idx:end_idx]
    var.long_name = d.variables['rmsa'].long_name
    var.units = str(d.variables['rmsa'].units)

    var = out.createVariable('cbh', 'f4', ('time', ))
    var[:] = d.variables['cbh'][beg_idx:end_idx]
    var.long_name = d.variables['cbh'].long_name
    var.units = str(d.variables['cbh'].units)

    var = out.createVariable('lwp', 'f4', ('time', ))
    var[:] = d.variables['lwp'][beg_idx:end_idx]
    var.long_name = d.variables['lwp'].long_name
    var.units = str(d.variables['lwp'].units)

    var = out.createVariable('qc_flag', 'i4', ('time', ))
    var[:] = d.variables['qc_flag'][beg_idx:end_idx]
    var.long_name = d.variables['qc_flag'].long_name
    var.units = str(d.variables['qc_flag'].units)
    """
    out.close()

if __name__ == '__main__':
    fn = sys.argv[1]
    out = sys.argv[2]
    run(fn, out)
