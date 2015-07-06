from netCDF4 import Dataset
import sys
#sys.path.append('/home/greg.blumberg/python_pkgs/')
from pylab import *
from sharppy.sharptab import interp, profile, params, thermo
import numpy as np
import scipy
import scipy.stats
from IPython.parallel import Client
import subprocess
import time
from datetime import datetime
import pandas as pd

from sharppy.sharptab.profile import BasicProfile
from sharppy.sharptab import params, interp

no_process = 8

class AERIProfile(BasicProfile):
    '''
    The AERI custom data class using the SHARPPy Profile framework.
    
    This class inherits from the Profile object.

    '''
    def __init__(self, **kwargs):
        '''
        Create the sounding data object
        
        Parameters
        ----------
        Mandatory Keywords
        pres : array_like
        The pressure values (Hectopaschals)
        hght : array_like
        The corresponding height values (Meters)
        tmpc : array_like
        The corresponding temperature values (Celsius)
        dwpc : array_like
        The corresponding dewpoint temperature values (Celsius)
            
        Optional Keyword Pairs (must use one or the other)
        wdir : array_like
        The direction from which the wind is blowing in
        meteorological degrees
        wspd : array_like
        The speed of the wind
        
        OR
            
        u : array_like
        The U-component of the direction from which the wind
        is blowing
            
        v : array_like
        The V-component of the direction from which the wind
        is blowing.
            
        Optional Keywords
        missing : number (default: sharppy.sharptab.constants.MISSING)
        The value of the missing flag

        Returns
        -------
        A profile object
        '''
        ## call the constructor for Profile
        super(AERIProfile, self).__init__(**kwargs)
        
        ## generate various parcels
        self.get_parcels()

        ## calculate thermodynamic window indices
        self.get_thermo()


    def get_parcels(self):
        '''
        Function to generate various parcels and parcel
        traces.
        Returns nothing, but sets the following
        variables:

        self.mupcl : Most Unstable Parcel
        self.sfcpcl : Surface Based Parcel
        self.mlpcl : Mixed Layer Parcel
        self.fcstpcl : Forecast Surface Parcel
        self.ebottom : The bottom pressure level of
            the effective inflow layer
        self.etop : the top pressure level of
            the effective inflow layer
        self.ebotm : The bottom, meters (agl), of the
            effective inflow layer
        self.etopm : The top, meters (agl), of the
            effective inflow layer
    
        Parameters
        ----------
        None

        Returns
        -------
        None
        '''

        self.mupcl = params.parcelx( self, flag=3 )
        if self.mupcl.lplvals.pres == self.pres[self.sfc]:
            self.sfcpcl = self.mupcl
        else:
            self.sfcpcl = params.parcelx( self, flag=1 )
        #self.fcstpcl = params.parcelx( self, flag=2 )
        self.mlpcl = params.parcelx( self, flag=4 )

        """
        ## get the effective inflow layer data
        self.ebottom, self.etop = params.effective_inflow_layer( self, mupcl=self.mupcl )

        ## if there was no effective inflow layer, set the values to masked
        if self.etop is ma.masked or self.ebottom is ma.masked:
            self.ebotm = ma.masked; self.etopm = ma.masked

        ## otherwise, interpolate the heights given to above ground level
        else:
            self.ebotm = interp.to_agl(self, interp.hght(self, self.ebottom))
            self.etopm = interp.to_agl(self, interp.hght(self, self.etop))
        """

    def get_thermo(self):
        '''
        Function to generate thermodynamic indices.
        
        Function returns nothing, but sets the following
        variables:

        self.k_idx - K Index, a severe weather index
        self.pwat - Precipitable Water Vapor (inches)
        self.lapserate_3km - 0 to 3km AGL lapse rate (C/km)
        self.lapserate_2km - 0 to 1km AGL lapse rate (C/km)
        self.lapserate_1km - 0 to 1km AGL lapse rate (C/km)
        self.convT - The Convective Temperature (F)
        self.maxT - The Maximum Forecast Surface Temp (F)
        self.mean_mixr - Mean Mixing Ratio
        self.totals_totals - Totals Totals index, a severe weather index

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        '''
        ## either get or calculate the indices, round to the nearest int, and
        ## convert them to strings.
        ## K Index
        self.k_idx = params.k_index( self )
        ## precipitable water
        self.pwat = params.precip_water( self )

        ## 0-3km agl lapse rate
        self.lapserate_3km = params.lapse_rate( self, 0., 3000., pres=False )
        ## 0-1km agl lapse rate
        self.lapserate_1km = params.lapse_rate( self, 0., 1000., pres=False )
        ## 850-500mb lapse rate
        self.lapserate_2km = params.lapse_rate( self, 0., 2000., pres=False )

        ## convective temperature
        self.convT = params.convective_temp( self ) 
        ## sounding forecast surface temperature
        self.maxT = params.max_temp( self )

        ## 100mb mean mixing ratio
        self.mean_mixr = params.mean_mixratio( self )
        ## calculate the totals totals index
        self.totals_totals = params.t_totals( self )

        self.ppbl_top = params.pbl_top(self)
        self.pbl_h = interp.to_agl(self, interp.hght(self, self.ppbl_top))
        self.sfc_rh = thermo.relh(self.pres[self.sfc], self.tmpc[self.sfc], self.dwpc[self.sfc])
        pres_sfc = self.pres[self.sfc]
        pres_1km = interp.pres(self, interp.to_msl(self, 1000.))
        self.pblrh = params.mean_relh(self, pbot=pres_sfc, ptop=self.ppbl_top)


def Td(pres, mxr):
    '''
        Dewpoint Conversion

        Description:
        This function converts mixing ratio to dewpoint.

        Parameters
        ----------
        pres : pressure (mb)
        mxr : water vapor mixing ratio (g/kg)

        Returns
        -------
        Td : dewpoint (K)
    '''
    L = (2.5*10.**6)
    R_v = 461.5
    T_o = 273.15
    e_so = 6.11
    mxr = mxr / 1000.
    e = (mxr*pres)/(.622+mxr)
    a = (L/R_v)
    b = (1./T_o)
    temp = np.power( (-1.) * ((1./a)*np.log(e/e_so) - b), -1)

    return temp

def makeIndicies(temp, dwpt, pres, height, cushon, parallel=False):
    """
        Input:
        temperature (2D array)
        dewpoint (2D array)
        pressure (2D array)
        height (1D array)

        Output:
        A dictionary containing various indicies for each profile passed through the
        input 2D arrays.
    """
    missing = np.ones(len(temp[0])) # This is the missing profile to be passed to the wind arrays in the Profile object
    
    if parallel == True:
        child = subprocess.Popen('ipcluster start -n ' + str(no_process), shell=True)
        #time.sleep(10)
        height = np.tile([height],(len(temp),1))
        dt = datetime.now()
        cli = Client()
        dview = cli[:]
        lbview = cli.load_balanced_view()
        lbview.block = True
        with dview.sync_imports():
            from sharppy.sharptab import interp, profile, params, thermo
            import indices_helper as sti

        results = lbview.map(makeProf, temp, dwpt, pres, height) # Returns a list of AERIProfile objects
        print datetime.now() - dt
    else:       
        from sharppy.sharptab import interp, profile, params, thermo

        results = []
        for i in range(len(temp)):
            prof = makeProf(temp[i], dwpt[i], pres[i], height)
            #lts = lowerTropStab(prof, 850.)
            results.append(prof)
    results = np.asarray(results, dtype=object)
    good_idx = np.where(results != -9999)[0]
    results = results[good_idx]
    indices, details = extractFields(results, (25,50,75), cushon)

    return indices, details

def hypsometric(temp, alt, sfc_press_ts):
    R = 287. # J/kg*K
    #note that temp is in celsius
    temp = temp + 273.15
    g = 9.81 #m/s^2

    pres_arr = -999*np.ones((temp.shape))
    pres_arr[:,0] = sfc_press_ts

    for l in np.arange(1,len(pres_arr.T),1):
        avg_temp = (temp[:,l] + temp[:,l-1])/2.
        delta_z = (alt[l-1] - alt[l])*(1000.) #To m from km
        a = (g/(R*avg_temp))
        p_2 = pres_arr[:,l-1]
        pres_arr[:,l] = p_2*np.exp(a*delta_z)

    return pres_arr

def monteCarlo(X,S,i):
    '''
        monteCarlo

        Description:
        This function takes in the X_op matrix from the AERIoe/MWRoe retrievals, the 
        retrieval covariance matrix S_op, and the dimension of how many profiles
        you want returned.

        Parameters
        ----------
        X : X_op matrix from AERIoe/MWRoe
        S : S_op matrix from AERIoe/MWRoe
        i : number of profiles generated by the Monte Carlo sampling

        Returns
        -------
        temp : the temperature array (C) of dimension (i,55)
        mxr : the water vapor mixing ratio array (g/kg) of dimension (i,55)
    '''
    Z = np.random.normal(0,1, (i, S.shape[1]))
    u,l,v = np.linalg.svd(S)
    Ssqrt = np.dot(np.dot(u, np.diag(np.sqrt(l))), v)
    Z_hat = np.dot(Z, Ssqrt) + X
    temp = Z_hat[:,:55]
    mxr = Z_hat[:,55:55+55]
    return temp, mxr

def extractFields(profs, percentiles, cushon):
    '''
        extractFields

        Description:
        This function extracts all of the necessary variables from a list of AERIProfile objects
        created by passing all of the profiles created via Monte Carlo sampling into the SHARPpy
        package via the AERIProfile object.

        Parameters
        ----------
        profs : a Numpy list of AERIProfile objects

        Returns
        -------
        Still need to think of this part!
    '''
    good_idx = np.where(profs != -9999)[0]
    profs = profs[good_idx]

    indices_dictionary = {}

    print "Extracting indices."
    parcels = ['mupcl', 'sfcpcl', 'mlpcl']

    prefixes = ['mu', 'sb', 'ml']
    prefix_descriptions = ['Most Unstable', 'Surface Based', '100 mb Mean Layer']
    pcl_att = ['bplus', 'b3km', 'bminus', 'li5', 'lfchght', 'lclhght', 'elhght', 'bmin', 'bminpres', 'cap', 'cappres', 'tmpc', 'dwpc', 'pres']
    att_name = ['cape', 'cape03', 'cin', 'li5', 'lfchght', 'lclhght', 'elhght', 'bmin', 'bminpres', 'cap', 'cappres', 'lpltmpc', 'lpldwpc', 'lplpres']
    descriptions = ['Convective Available Potential Energy', "0-3 km Convective Available Potential Energy", "Convective Inhibition",\
                   '500 mb Lifted Index', 'Level of Free Convection Height', "Lifted Condensation Level Height", "Equilibrium Level", "Buoyancy Minimum Value", "Buoyancy Minimum Pressure", \
                   'Cap Strength', "Pressure at Maximum Cap Strength", "Lifted Parcel Level Temperature", "Lifted Parcel Level Dewpoint",\
                   "Lifted Parcel Level Pressure"]
    units = ['J/kg', 'J/kg', 'J/kg', 'C', 'meters AGL', 'meters AGL', 'meters AGL', 'C', 'mb', 'C', 'mb', 'C', 'C', 'mb']

    var_details = {}

    # Add all parcel information to the dictionaries.
    for i in xrange(len(profs)):
        prof = profs[i]
        for p, prefix, prefix_desc in zip(parcels, prefixes, prefix_descriptions):
            pcl = getattr(prof, p)
            for att, name, desc, unit in zip(pcl_att, att_name, descriptions, units):
                try:
                    indices_dictionary[prefix + name] = indices_dictionary[prefix + name] + [getattr(pcl, att)]
                except:
                    indices_dictionary[prefix + name] = [getattr(pcl, att)]
                var_details[prefix + name] = [prefix_desc + ' ' + desc, unit]


    var_names = ['k_idx', 'pwat', 'lapserate_1km', 'lapserate_2km', 'lapserate_3km', 'convT', 'maxT', 'totals_totals',\
                 'ppbl_top']
    var_descriptions = ['K-Index', "Precipitable Water Vapor", '0-1 km Lapse Rate', '0-2 km Lapse Rate', '0-3 km Lapse Rate',\
                        'Convective Temperature', 'Maximum Forecasted Temperature', 'Total Totals', 'Planetary Boundary Layer Top']
    var_units = ['unitless', 'inches', 'C/km', 'C/km', 'C/km', 'C', 'C', 'unitless', 'mb']
    # Add additional information (lapse rates, etc.) to the dictionaries
    for i in xrange(len(profs)):
        prof = profs[i]
        for att, desc, unit in zip(var_names, var_descriptions, var_units):
            try:
                indices_dictionary[att] = indices_dictionary[att] + [getattr(prof, att)]
            except:
                indices_dictionary[att] = [getattr(prof, att)]
            var_details[att] = [desc, unit]

    dt = datetime.now()
    for dic in indices_dictionary.keys():
        filtered_indices = np.ma.masked_invalid(indices_dictionary[dic])
        if len(filtered_indices) > cushon and len(filtered_indices[~filtered_indices.mask]) != 0:
            #indices_dictionary[dic] = np.nanpercentile(filtered_indices, percentiles)
            indices_dictionary[dic] = np.percentile(filtered_indices[~filtered_indices.mask], percentiles)
        else:
            indices_dictionary[dic] = [-9999,-9999,-9999]
    print "Time to sort indices:", datetime.now() - dt

    return indices_dictionary, var_details
    

def makeProf(temp, dwpt, pres, height):
    missing = np.ones((len(height)))
    #prof = AERIProfile(pres=pres, hght=height, tmpc=temp, dwpc=dwpt, wdir=missing, wspd=missing)
    try:
        prof = AERIProfile(pres=pres, hght=height, tmpc=temp, dwpc=dwpt, wdir=missing, wspd=missing)
    except Exception,e:
        print e
    #    prof = e
        prof = -9999

    return prof


def makeIndicesErrors(X, S, height, pressure, num_perts, cushon, flag):#, keys):
    temp = X

    sblcls = -999*np.ones((X.shape[0], 3))
    sblfcs = -999*np.ones((X.shape[0], 3))
    sbcapes = -999*np.ones((X.shape[0], 3))
    sbcins = -999*np.ones((X.shape[0], 3))
    sbcape03s = -999*np.ones((X.shape[0], 3))

    mllcls = -999*np.ones((X.shape[0], 3))
    mllfcs = -999*np.ones((X.shape[0], 3))
    mlcapes = -999*np.ones((X.shape[0], 3))
    mlcins = -999*np.ones((X.shape[0], 3))
    mlcape03s = -999*np.ones((X.shape[0], 3))

    mulcls = -999*np.ones((X.shape[0], 3))
    mulfcs = -999*np.ones((X.shape[0], 3))
    mucapes = -999*np.ones((X.shape[0], 3))
    mucins = -999*np.ones((X.shape[0], 3))
    mucape03s = -999*np.ones((X.shape[0], 3))

    pws = -999 * np.ones((X.shape[0], 3))
    ctss = -999 * np.ones((X.shape[0], 3))
    #lr03s = -999*np.ones(X.shape[0])
    ltss = -999.* np.ones((X.shape[0], 3))

    print "\nThis calculation has " + str(len(X)) + " samples."

    all_indices = np.empty((len(X)), dtype=dict)
    for i in range(len(X)): # Loop for all of the indices 
        print "Beginning Convective Index Calculation of Sample: " + str(i)

        try:
            temp, mxr = monteCarlo(X[i], S[i], num_perts)
        except Exception as e:
            print "\tMonte Carlo sampling of Sop failed."
            print e
            print "\t" + e
            continue
        print "\tSuccess performing the Monte Carlo sampling."

        sfc_pres = np.repeat(pressure[i][0], num_perts)
        pres = hypsometric(temp, height/1000., sfc_pres)
        dwpt = Td(pres, mxr)-273.15

        #Take care of super saturated profiles that show up
        sat_idx = np.where(temp-dwpt <= 0)
        dwpt[sat_idx] = temp[sat_idx] - .001
        
        indices, details = makeIndicies(temp, dwpt, pres, height, cushon, parallel=True)
        all_indices[i] = indices
   
    return all_indices, details
