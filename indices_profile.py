import sharppy.sharptab.params as params
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import numpy as np
from netCDF4 import Dataset
from io import StringIO
from datetime import datetime
 
def add_qc(ncdf, dictionary):
    for qc_key in ['lwp', 'rmsa', 'rmsr', 'cbh', 'converged_flag']:
        var = ncdf.createVariable(qc_key, 'f4', ('time',))
        var[:] = dictionary[qc_key]
    
    ncdf.comment = dictionary['comment']
    ncdf.input_file = dictionary['input_file'] 
    return ncdf

def writeIndexProfiles(filename, dictionary):
    ncdf = Dataset(filename, 'w', format='NETCDF3_CLASSIC')
    ncdf.createDimension('height', len(dictionary['hght']))
    ncdf.createDimension('time', len(dictionary['cape']))
    ncdf.createDimension('single', 1)
    ncdf.date_created = datetime.strftime(datetime.utcnow(), "%Y-%m-%d %H:%M:%S UTC")
    
    var = ncdf.createVariable('base_time', 'f4', ('single',))
    var[:] = dictionary['bt']
     
    var = ncdf.createVariable('time_offset', 'f4', ('time',))
    var[:] = dictionary['to']
    
    var = ncdf.createVariable('lat', 'f4', ('single',))
    var[:] = dictionary['lat']
     
    var = ncdf.createVariable('lon', 'f4', ('time',))
    var[:] = dictionary['lon']

    try:
        # add the QC data
        ncdf = add_qc(ncdf, dictionary)
    except:
        print("NO QC DATA ADDED TO THIS FILE")

    var = ncdf.createVariable('height', 'f4', ('height',))
    var[:] = dictionary['hght']
    var.units = 'm MSL'

    var = ncdf.createVariable('pressure', 'f4', ('time','height',))
    var[:] = dictionary['pres']
    var.units = 'mb'

    var = ncdf.createVariable('cape', 'f4', ('time','height',))
    var[:] = dictionary['cape']
    var.units = 'J/kg'


    var = ncdf.createVariable('cin', 'f4', ('time','height',))
    var[:] = dictionary['cin']
    var.units = 'J/kg'


    var = ncdf.createVariable('lfchght', 'f4', ('time','height',))
    var[:] = dictionary['lfchght']
    var.units = 'm AGL'


    var = ncdf.createVariable('lclhght', 'f4', ('time','height',))
    var[:] = dictionary['lclhght']
    var.units = 'm AGL'

    var = ncdf.createVariable('li5', 'f4', ('time','height',))
    var[:] = dictionary['li5']
    var.units = 'C'

    var = ncdf.createVariable('bmin', 'f4', ('time','height',))
    var[:] = dictionary['bmin']
    var.units = 'C'

    var = ncdf.createVariable('bminhght', 'f4', ('time','height',))
    var[:] = dictionary['bminhght']
    var.units = 'm AGL'

    ncdf.close()


def writeBuoyancyFile(filename, dictionary):
    ncdf = Dataset(filename, 'w', format='NETCDF3_CLASSIC')
    ncdf.createDimension('height', len(dictionary['hght']))
    ncdf.pcl_type = dictionary['pcl_type']    
    ncdf.date_created = datetime.strftime(datetime.utcnow(), "%Y-%m-%d %H:%M:%S UTC")
    # add the QC data
    ncdf = add_qc(ncdf, dictionary)

    var = ncdf.createVariable('height', 'f4', ('height',))
    var[:] = dictionary['hght']
    var.units = 'm MSL'
   
    var = ncdf.createVariable('buoyancy_profile', 'f4', ('height',))
    var[:] = dictionary['buoyancy_profile']
    var.units = 'C'
    var.comment = "virtual temperature correction applied"

    var = ncdf.createVariable('pressure', 'f4', ('height',))
    var[:] = dictionary['pres']
    var.units = 'mb'

    ncdf.close()

def index_profiles(prof):
    # performs the SHARPpy lifting routines on each level of a sounding and then
    # generates profiles of the indices calculated from each parcel
    # returns a dictionary with those profiles along with information about the height
    # and pressure grids they are on
    # prof is a SHARPpy profile
    top_idx = np.where(interp.to_agl(prof, prof.hght) >= 4000)[0][0]
    top_idx = prof.top + 1
    bot_idx = prof.sfc
    all_parcels = []
    limit = np.where(prof.hght > 4000)[0][0]
    for i in range(bot_idx, top_idx, 1):
        if i > limit:
            continue
        pcl = params.parcelx(prof, tmpc=prof.tmpc[i], dwpc=prof.dwpc[i], pres=prof.pres[i])
        all_parcels.append(pcl)

    parcels = list(range(bot_idx, top_idx, 1))
    lclhght = np.ones(len(parcels)) * -9999
    cape = np.ones(len(parcels)) * -9999
    cin = np.ones(len(parcels)) * -9999
    lfchght = np.ones(len(parcels)) * -9999
    bmin = np.ones(len(parcels)) * -9999
    bminhght = np.ones(len(parcels)) * -9999
    limax = np.ones(len(parcels)) * -9999
    li5 = np.ones(len(parcels)) * -9999
    limaxhght = np.ones(len(parcels)) * -9999
    for i in range(len(all_parcels)):
        pcl = all_parcels[i]
        lclhght[i] = pcl.lclhght
        lfchght[i] = pcl.lfchght
        cape[i] = pcl.bplus
        cin[i] = pcl.bminus
        bmin[i] = pcl.bmin
        bminhght[i] = interp.to_agl(prof, interp.hght(prof, pcl.bminpres))
        limax[i] = pcl.limax
        li5[i] = pcl.li5
        limaxhght[i] = interp.to_agl(prof, interp.hght(prof, pcl.limaxpres))
    limax = np.ma.masked_invalid(limax)
    print(cape.shape)
    data = {}
    data['pcl'] = all_parcels
    data['hght'] = prof.hght[list(range(bot_idx, top_idx, 1))]
    data['pres'] = prof.pres[list(range(bot_idx, top_idx, 1))]
    data['cape'] = cape
    data['cin'] = cin
    data['lfchght'] = lfchght
    data['lclhght'] = lclhght
    data['bmin'] = bmin
    data['bminhght'] = bminhght
    data['limax'] = limax
    data['li5'] = li5
    data['limaxhght'] = limaxhght

    return data

def buoyancy_profile(prof, pcl_flag=1):
    # returns the buoyancy profile for a specific parcel type
    # flag = 1 (sfc pcl)
    # flag = 3 (most-unstable)
    # flag = 4 (100 mb mean layer)
    pcl = params.parcelx(prof, flag=pcl_flag)

    ttrace = pcl.ttrace
    ptrace = pcl.ptrace
    height = interp.hght(prof, ptrace)
   
    #interpolate to the prof grid (new grid)
    ttrace_ng = np.interp(prof.hght, height, ttrace)
    ttrace_ng = np.ma.masked_where(prof.hght < interp.hght(prof, pcl.lplvals.pres), ttrace_ng)
    data = {}
    data['hght'] = prof.hght
    data['buoyancy_profile'] = ttrace_ng - prof.vtmp
    data['pres'] = prof.pres
        
    # Tell the user which parcel
    if pcl_flag == 1:
        pcl_type = 'sb'
    elif pcl_flag == 3:
        pcl_type = 'mu'
    elif pcl_flag == 4:
        pcl_type = "ml"
    data['pcl_type'] = pcl_type
    
    return data

def parseSPC(spc_file):
    data = np.array([l.strip() for l in spc_file.split('\n')])
   
    title_idx = np.where( data == '%TITLE%' )[0][0]
    start_idx = np.where( data == '%RAW%')[0] + 1
    finish_idx = np.where( data == '%END%')[0]

    data_header = data[title_idx + 1].split()
    location = data_header[0]
    time = data_header[1][:11]
    print(title_idx, start_idx, finish_idx)
    print(data[start_idx : finish_idx])
    full_data = '\n'.join(data[start_idx : finish_idx][:])
    sound_data = StringIO(full_data)

    p, h, T, Td, wdir, wspd = np.genfromtxt( sound_data, delimiter=',', comments='%', unpack=True)
    prof = profile.create_profile(profile='default', pres=p, hght=h, tmpc=T, dwpc=Td, wdir=wdir, wspd=wspd,\
                                  missing=-9999, strictQC=True)
    return prof

#def writeIndices(buoyancy_data, index_profile_data):
     

def main():
    # Test code here to make sure that the two functions works
    spc_file = open('sample_profiles/14061619.OAX', 'r').read()
    prof = parseSPC(spc_file)
    print(buoyancy_profile(prof, 1))
    print(index_profiles(prof))

if __name__ == "__main__":
    main()
