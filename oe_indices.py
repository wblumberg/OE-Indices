print('python oe_indices [-a or -w] <input_file>')

import sys
sys.path.append('/home/wblumberg/.conda/envs/local/site-packages/')
#sys.path.append('/home/greg.blumberg/python_pkgs/')
#from pylab import *
from netCDF4 import Dataset, MFDataset, num2date
import indices_helper as sti
import sys
from datetime import datetime
import sharppy
import warnings
import numpy as np
warnings.filterwarnings("ignore")

'''
    oe_indices.py

    AUTHOR: Greg Blumberg (OU/CIMMS)

    This script takes in OE-retrievals of thermodynamic profiles from
    either AERIoe or MWRoe and performs a.) Monte Carlo sampling of the retrievals
    to generate approximately 500 retrievals, and then b.) uses the SHARPpy
    libraries to perform convection index calculations on all 500 retrievals.
    After this the convection indices are sorted, the 25th, 50th, and 75th percentiles
    are pulled out for each index and then saved in a netCDF file.

    Running this script requires the netCDF4-python libraries, the SHARPpy libraries,
    and IPython to perform the calculations in an efficient manner.

    To run this script, use this command:
        python oe_indices.py [-a or -w] <input_file>

    where <input_file> is the path to the OE netCDF file and the flag -a means append to an existing file
    and the -w means overwrite any existing file.

    Before creating a new OE indices netCDF file, this program first checks for an existing file
    because it may need to be appended to (for example if this script is running in real-time.

'''
 
def run(fn, out):
    #flag = sys.argv[1]
    fn = sys.argv[1]
    out = sys.argv[2]
    name = fn.split('/')[-1].replace('aerioe1turn', 'aerioeidx1blum')
    name = name.replace('c1', 'c2')
    name = out + '/' + name
    d = Dataset(fn)
    bt = d.variables['base_time'][:]
    to = d.variables['time_offset'][:]

    retrieval_times = num2date(bt+to, 'seconds since 1970-01-01 00:00') 

    #existing_file = glob.glob(name)
    #if len(existing_file) == 1 and flag == '-a':
    #    temporary = Dataset(existing_file[0])
    #    temp_bt = temporary.variables['base_time'][:] + temporary.variables['time_offset'][:]
    #    beg_idx = len(to) - len(temp_bt)
    #    end_idx = len(to)
    #elif flag == '-w':
    beg_idx = 0
    end_idx = len(to)

    to = to[beg_idx:end_idx]
    height = d.variables['height'][:]
    pres = d.variables['pressure'][beg_idx:end_idx]
    Xop = d.variables['Xop'][beg_idx:end_idx]
    Sop = d.variables['Sop'][beg_idx:end_idx]
    converged_flag = d.variables['converged_flag'][beg_idx:end_idx]
    height = height * 1000.

    print("Generating output file:", name)
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
    var.long_name = d.variables['base_time'].long_name
    var.units = str(d.variables['base_time'].units)

    var = out.createVariable('time_offset', 'f4', ('time', ))
    var[:] = to
    var.long_name = d.variables['time_offset'].long_name
    var.units = str(d.variables['time_offset'].units)

    print("Copying over AERIoe QC flags to the output file...")
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

    num_perts = 500 # The ideal number of profiles to compute the index percentiles
    cush = 100 # The number of profiles where indices MUST be computed from (allows for crashing)
 
    print("-----------------------------------------------------------")
    print("Starting to generate the errors for the convective indices.")
    print("-----------------------------------------------------------")
    for i in range(len(Xop)): # Loop over each retrieval
        print("Starting with " + retrieval_times[i].strftime("%Y-%m-%d %H:%M:%S UTC") + " Retrieval No. " + str(i) + "/" + str(len(Xop)))
        ideb, details = sti.makeIndicesErrors(Xop[i], Sop[i], height, pres[i], num_perts, cush, converged_flag)
        print("Completed Monte Carlo sampling of the convective indices.  Saving to file...")
        for var_name in np.sort(list(details.keys())): # Loop over all of the index names
            if i == 0: # if we're working on the first retrieval, make the netCDF variables
                # Obtain the index description and units
                longname = details[var_name][0]
                unit = details[var_name][1]
                var = out.createVariable(var_name, 'f4', ('time', 'samples',))
                var.long_name = longname
                var.units = unit
            else:
                var = out.variables[var_name]
            data = np.empty((len(ideb),num_perts))
            for i in range(len(ideb)):
                data_subset = ideb[i][var_name]
                if len(data_subset) > num_perts:
                    data_subset = data_subset[:num_perts]
                data_subset = np.concatenate((data_subset, np.ones(num_perts - len(data_subset))*-9999))
            data[i,:] = data_subset
        var[i,:] = data
        

    #if flag == '-w' and len(glob.glob(name)) != 0:
    #    os.system('rm ' + name)
    #else:
    #    print "Appending to the file."

    #var = out.createVariable('percentiles', 'f4', ('percentiles', ))
    #var[:] = [25,50,75]
    #var.long_name = "Percentiles for the errorbars"
    #var.units = "%"

    # loop through the indices and add them to the netCDF file
   out.close()

if __name__ == '__main__':
    fn = sys.argv[1]
    out = sys.argv[2]
    run(fn, out)
