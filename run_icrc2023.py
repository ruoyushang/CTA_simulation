#!/usr/bin/python3

# ctobssim draws random events using a model of the celestial gamma-ray intensity distribution 
# that was convolved with the instrument response function.
# The tool adds random events from the expected distribution of the residual background to the data.

import os
import sys
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import gammalib
from datetime import datetime


scriptdir = "/nevis/tehanu/home/capasso/Analysis/CTA/scripts/dribeiro/"

layout = sys.argv[1]
expo_hrs = float(sys.argv[2])
model_type = sys.argv[3]

outdir = f'/a/data/serret/shang/condor_output/CTAsim_icrc2023_{layout}_{int(expo_hrs)}hr_{model_type}'
print(outdir)
try:
    os.makedirs(outdir)
except FileExistsError:
    pass
os.chdir(outdir)

fitting_only = False
#fitting_only = True

################################### SIMULATION #######################################

center_string = '05h34m30.9s +22d00m44.5s'#from http://tevcat.uchicago.edu/?mode=1;id=84

source_name = f"crab_nebula"
inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.

if 'two' in model_type:
    source_name = f"two_extended_sources"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
elif 'dragonfly' in model_type:
    source_name = f"dragonfly_xray"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '20h19m25.0s +36d48m14s'#from http://tevcat.uchicago.edu/?mode=1;id=84
elif 'rx_j1713_0p25kpc' in model_type:
    source_name = f"rx_j1713.7-3946_gas_0p25kpc"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '17h13m33.6s -39d45m36s'#from http://tevcat.uchicago.edu/?mode=1;id=84
elif 'rx_j1713_0p5kpc' in model_type:
    source_name = f"rx_j1713.7-3946_gas_0p5kpc"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '17h13m33.6s -39d45m36s'#from http://tevcat.uchicago.edu/?mode=1;id=84
elif 'rx_j1713_1p0kpc' in model_type:
    source_name = f"rx_j1713.7-3946_gas_1p0kpc"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '17h13m33.6s -39d45m36s'#from http://tevcat.uchicago.edu/?mode=1;id=84
elif 'rx_j1713_2p0kpc' in model_type:
    source_name = f"rx_j1713.7-3946_gas_2p0kpc"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '17h13m33.6s -39d45m36s'#from http://tevcat.uchicago.edu/?mode=1;id=84
elif 'rx_j1713_3p0kpc' in model_type:
    source_name = f"rx_j1713.7-3946_gas_3p0kpc"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '17h13m33.6s -39d45m36s'#from http://tevcat.uchicago.edu/?mode=1;id=84
elif 'rx_j1713_4p0kpc' in model_type:
    source_name = f"rx_j1713.7-3946_gas_4p0kpc"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '17h13m33.6s -39d45m36s'#from http://tevcat.uchicago.edu/?mode=1;id=84
elif 'ic443' in model_type:
    source_name = f"ic443_optical"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '06h16m51.13886s +22d30m24.0996s'
elif 'lmc' in model_type:
    source_name = f"lmc_radio"
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/{source_name}.xml' # Input model XML file.
    center_string = '05h38m42.396s -69d06m03.36s'
elif 'gal_center' in model_type:
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/models/model_galactic_pwn.xml' # Input model XML file.
    center_string = '18h07m45.72774s -20d17m24.0976s' # Gal (lon,lat) = (10,0)
if 'bkg' in model_type:
    inmodel = f'/nevis/tehanu/home/ryshang/CTA_sim/Ruo/condor_sub/icrc2023/cta_background.xml' # Input model XML file.

random_seed = int(datetime.now().timestamp())

center = SkyCoord(center_string, frame='icrs') 
#side_skymap_x = 1.5
#side_skymap_y = 1.5
side_skymap_x = 6.0
side_skymap_y = 6.0
xref_skymap = center.ra.deg
yref_skymap = center.dec.deg

print ('inmodel = %s'%(inmodel))
obsdef_xml = "obsdef.xml"
logfile_obsdef = "csobsdef.log"

instrument = layout
#if layout=='f4':
#    instrument = "mult4_zd20_f4" # F4 layout
#if layout=='m4':
#    instrument = "mult3_zd20_m4interp" # M4 layout
#if layout=='c0':
#    instrument = "mult4_zd20_c0"
#if layout=='f5':
#    instrument = "mult4_zd20_f5"
caldb = instrument

#sys.exit("End simulation.")

irf = "South_50h" # Instrumental response function.
radius_obs = 4
emin_obs = 0.1 # min energy threshold TeV
emax_obs = 100.0 # max energyy threshold TeV
duration = expo_hrs*60*60 # duration of each pointing in seconds
#n_obs = int(expo_hrs*60.*60./duration) # number of runs
n_obs=4 # number of runs
ra_offset_obs = 0. 
dec_offset_obs = 0.

outevents_obs = "events.xml"
logfile_obs = "ctobssim.log"
nthreads = 0
#nthreads = 10
    
#sys.exit("End simulation.")

pointing_file = 'pointings.txt'
# open file
f = open(pointing_file, 'w')
# header
f.write('id,ra,dec\n')
# pointings
for i in range(n_obs):

    if(i%4==0):
        ra_obs = center.ra.deg - ra_offset_obs
        dec_obs = center.dec.deg

    if(i%4==1):
        ra_obs = center.ra.deg + ra_offset_obs
        dec_obs = center.dec.deg

    if(i%4==2):
        ra_obs = center.ra.deg
        dec_obs = center.dec.deg + dec_offset_obs

    if(i%4==3):
        ra_obs = center.ra.deg
        dec_obs = center.dec.deg - dec_offset_obs

    obs_num=i+1
    f.write(f'{obs_num},{ra_obs:.2f},{dec_obs:.2f}\n')
# close file
f.close()


cmd = f'python {scriptdir}/ObsDef.py --inpnt {pointing_file} --outobs {obsdef_xml} --caldb {caldb} --irf {irf} --emin {emin_obs} --emax {emax_obs} --duration {duration} --rad {radius_obs} --logfile {logfile_obsdef} --clobber'
print(cmd)
if not fitting_only:
    os.system(cmd)

#sys.exit("End simulation.")

#obsdef_xml = "/a/data/serret/shang/1dc/obs/obs_gps_baseline_ruo_short.xml"

cmd = f'python {scriptdir}/Simulation.py --inobs {obsdef_xml} --inmodel {inmodel} --logfile {logfile_obs} --outevents {outevents_obs} --clobber --nthreads {nthreads} --seed {random_seed}'
print(cmd)
if not fitting_only:
    os.system(cmd)

#sys.exit("End simulation.")



################################### EVENT SELECTION #######################################

# This created a new event list FITS file selected_events.fits that contains the selected data. 

##inobs_select = 'xml_obs' # Input event list or observation definition XML file
#radius_select = 4
#ra_select = center.ra.deg
#dec_select = center.dec.deg
#emin_select = 0.2
#emax_select = 20.0
## Start/stop time for event selection (UTC string, JD, MJD or MET in seconds). Start times given in MET seconds are counted with respect to the time reference of the input observation(s). If INDEF, NONE, UNDEF or UNDEFINED is passed as value, no time selection will be performed.
#tmin_select = 0 
#tmax_select = 30*60
#usepnt_select = "yes"
#
#logfile_select = "ctselect_emin={0:.2f}_emax={1:.2f}_radselect={2:.2f}.log".format(emin_select,emax_select,radius_select) 
#prefix_select = "selected_emin={0:.2f}_emax={1:.2f}_radselect={2:.2f}".format(emin_select,emax_select,radius_select)
#outobs_select = "{0}events.xml".format(prefix_select)
#
##cmd = 'python '+scriptdir+'/Select.py --inobs {0} --rad {1} --ra {2:.2f} --dec {3} --emin {4} --emax {5} --tmin {6} --tmax {7} \
##--logfile {8} --outobs {9} --clobber --prefix {10} --usepnt {11}'.format(inobs_select,radius_select,ra_select,dec_select,
##                                             emin_select,emax_select,tmin_select,tmax_select,
##                                            logfile_select,outobs_select,prefix_select,usepnt_select)
#cmd = 'python '+scriptdir+'/Select.py --rad {0} --ra {1} --dec {2} --emin {3} --emax {4} --tmin {5} --tmax {6} \
#--logfile {7} --outobs {8} --clobber --prefix {9} --usepnt {10}'.format(radius_select,ra_select,dec_select,
#                                             emin_select,emax_select,tmin_select,tmax_select,
#                                            logfile_select,outobs_select,prefix_select,usepnt_select)
#print(cmd)
#os.system(cmd)


################################### SKYMAP PRODUCTION #######################################

inobs_skymap = outevents_obs
coordsys = "CEL"
#coordsys = "GAL"
proj = "CAR"
emin_skymap = emin_obs
emax_skymap = emax_obs
usepnt_skymap = "no"
bkgsubtract = "IRF" # Options: NONE, IRF, RING
binsz_skymap = 0.02
#side_skymap = radius_obs*np.sqrt(2)
nxpix_skymap = int(np.floor(side_skymap_x/binsz_skymap)) #round always to the floor, to make sure 
nypix_skymap = int(np.floor(side_skymap_y/binsz_skymap)) #round always to the floor, to make sure 
    
outmap_skymap = f'skymap_selected_events_emin={emin_skymap:.2f}_emax={emax_skymap:.2f}_side={side_skymap_x:.2f}_bkgsub={bkgsubtract}.fits'
logfile_skymap = f'ctskymap_selected_events_emin={emin_skymap:.2f}_emax={emax_skymap:.2f}_side={side_skymap_x:.2f}_bkgsub={bkgsubtract}.log'
cmd = f'python {scriptdir}/Skymap.py --inobs {inobs_skymap} --caldb {caldb} --irf {irf} --outmap {outmap_skymap} --coordsys {coordsys} --proj {proj} --emin {emin_skymap} --emax {emax_skymap} --xref {xref_skymap} --yref {yref_skymap} --bkgsubtract {bkgsubtract} --logfile {logfile_skymap} --binsz {binsz_skymap} --nxpix {nxpix_skymap} --nypix {nypix_skymap} --usepnt {usepnt_skymap} --clobber'
print(cmd)
if not fitting_only:
    os.system(cmd)

#sys.exit("End simulation.")

################################### BINNING EVENT DATA #######################################
# A counts cube is a 3-dimensional data cube that is spanned by Right Ascension (or Galactic longitude), Declination (or Galactic latitude), and energy (by default logarithmically spaced, but this is under your control).

inobs_bins = outevents_obs
ebinalg_bins = "LOG"
xref_bins = xref_skymap
yref_bins = yref_skymap
emin_bins = emin_obs
emax_bins = emax_obs
enumbins_bins = int(np.round(np.abs(np.log10(emax_bins/emin_bins)))*25) 
#--> I got a warning that 10 bins where too few and I should take at least 25/decade, which for the energy range 0.2-10.0 means 42, 0.2-1.0 means 17
usepnt_bins = "no"
binsz_bins = 0.02
nxpix_bins = int(np.floor(side_skymap_x/binsz_bins)) #round always to the floor, to make sure 
nypix_bins = int(np.floor(side_skymap_y/binsz_bins)) #round always to the floor, to make sure 
    
outobs_bins = f'cntcube_emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}.fits'
logfile_bins = f'ctbin_emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}.log'

prefix_bins = f'cntcube_emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}'

cmd = f'python {scriptdir}/Binning.py --inobs {inobs_bins} --outobs {outobs_bins} --ebinalg {ebinalg_bins} --coordsys {coordsys} --proj {proj} --emin {emin_bins} --emax {emax_bins} --enumbins {enumbins_bins} --xref {xref_bins} --yref {yref_bins} --logfile {logfile_bins} --clobber --prefix {prefix_bins} --usepnt {usepnt_bins} --binsz {binsz_bins} --nxpix {nxpix_bins} --nypix {nypix_bins}'
    
print(cmd)
if not fitting_only:
    os.system(cmd)
    
sys.exit("End simulation.")

#################################### Pre-computing the binned response #######################################
# Pre-computing the instrument response will speed-up the model fitting later.
# To speed-up the model fitting it is recommended to pre-compute the instrument response for the counts cube. 
# Specifically, you should compute an exposure cube, a point spread function cube and a background cube.
# The exposure cube is computed using the ctexpcube tool and provides the effective area multiplied by the livetime of the observation. 

#inobs_exp = outevents_obs
inobs_exp = obsdef_xml
incube_exp = outobs_bins # exposure cube definition is copied from counts cube
outcube_exp = f'expcube_emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}.fits'
logfile_exp = f'ctexpcube_emin={emin_bins:.2f}_emax={emax_bins:.2f}_radselect={side_skymap_x:.2f}.log'

cmd = f'python {scriptdir}/Exposure.py --inobs {inobs_exp} --incube {incube_exp} --caldb {caldb} --irf {irf} --outcube {outcube_exp} --logfile {logfile_exp} --clobber --usepnt {usepnt_bins}'

print(cmd)
if not fitting_only:
    os.system(cmd)

# Next, we use the ctpsfcube tool to compute the point spread function (PSF) cube for the counts cube.
# This produces the FITS file psfcube.fits that contains the point spread function as function of sky position and energy. 

#inobs_psf = outevents_obs
inobs_psf = obsdef_xml
binsz_psf = 1.
nxpix_psf = 10
nypix_psf = 10
ebinalg_psf = "LOG"
enumbins_psf = enumbins_bins #what is the optimal? Keeping same as ctbin for the moment
outcube_psf = f'psfcube_emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}.fits'
logfile_psf = f'ctpsfcube_emin={emin_bins:.2f}_emax={emax_bins:.2f}_radselect={side_skymap_x:.2f}.log'

cmd = f'python {scriptdir}/PSF.py --inobs {inobs_psf} --caldb {caldb} --irf {irf} --coordsys {coordsys} --proj {proj} --xref {xref_bins} --yref {yref_bins} --binsz {binsz_psf} --nxpix {nxpix_psf} --nypix {nypix_psf} --ebinalg {ebinalg_psf} --emin {emin_bins} --emax {emax_bins} --enumbins {enumbins_psf} --outcube {outcube_psf} --logfile {logfile_psf} --clobber --usepnt {usepnt_bins}'
print(cmd)
if not fitting_only:
    os.system(cmd)

# Finally, we use the ctbkgcube tool to compute the background cube.
# This produces the FITS file bkgcube.fits that contains the predicted background rate as function of sky position and energy. 
# The tool also produces the model definition file models.xml on output that will serve as input for the maximum likelihood analysis that will follow. 

#inobs_bkg = outevents_obs
inobs_bkg = obsdef_xml
incube_bkg = outobs_bins
outcube_bkg = f'bkgcube_emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}.fits'
outmodel_bkg = f'emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}_models.xml'
logfile_bkg = f'ctbkgcube_emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}.log'

cmd = f'python {scriptdir}/Background.py --inobs {inobs_bkg} --incube {incube_bkg} --inmodel {inmodel} --caldb {caldb} --irf {irf} --outcube {outcube_bkg} --outmodel {outmodel_bkg} --logfile {logfile_bkg} --clobber --seed {random_seed}'
print(cmd)
if not fitting_only:
    os.system(cmd)


################################### Fitting binned data #######################################
# You will learn how to adjust a parametrised model to the counts cube by means of a maximum likelihood fit.
# The data are fitted in binned mode, which means that the events have been binned into a 3D counts cube and the fit computes the log-likelihood function by summing over all 200 x 200 x 20 bins of the counts cube. 
# The ctlike tool produced a model definition XML output file crab_results.xml that contains the fit results. 
# The file is a copy of the model definition XML input file models.xml where the parameter values were replaced by the fit results. 
# In addition, the statistical uncertainties were added for each fitted parameter using the attribute error.

# The ctlike tool has the ability to estimate the detection significance for sources in the XML model. 
# This is done by computing the Test Statistic value which is defined as twice the log-likelihood difference between fitting a source at a given position on top of a (background) model or fitting no source. 
# As a rule of thumb, the square root of the Test Statistic value gives the source detection significance in Gaussian sigmas, although the actual conversion depends somewhat on the formulation of the statistical problem and the number of degrees of freedom associated with the source.

# To instruct ctlike to compute the Test Statistic value for a given source you need to add the attribute tscalc="1" to the XML file


#model_name_lep = "RX J1713.7-3946"
#model_name_bkg = "BackgroundModel"
model_name_lep = "Dragonfly PWN"
model_name_bkg = "CTABackgroundModel"


spatial_type = "LepOnly"
data_model_container = gammalib.GModels(outmodel_bkg) #--> can try to have bkg cube as model
bkgd_model = data_model_container[model_name_bkg]
init_model_lep = data_model_container[model_name_lep]

var_name = 'Prefactor'
bkgd_model[var_name].fix()
var_name = 'Index'
bkgd_model[var_name].fix()
var_name = 'PivotEnergy'
bkgd_model[var_name].fix()

var_name = 'Prefactor'
init_model_lep[var_name].value(init_model_lep[var_name].value())
init_model_lep[var_name].free()
var_name = 'Index'
init_model_lep[var_name].value(init_model_lep[var_name].value())
init_model_lep[var_name].free()
var_name = 'CutoffEnergy'
init_model_lep[var_name].value(init_model_lep[var_name].value())
init_model_lep[var_name].fix()
var_name = 'PivotEnergy'
init_model_lep[var_name].value(init_model_lep[var_name].value())
init_model_lep[var_name].fix()
var_name = 'Normalization'
init_model_lep[var_name].value(init_model_lep[var_name].value())
init_model_lep[var_name].fix()

init_model_lep.tscalc(True)

inmodel_xml = f'{outmodel_bkg[:-4]}_{spatial_type}_tscalc.xml' # this is the initial values of the fit function
data_model_container.save(inmodel_xml)

inmodel_like = inmodel_xml
outmodel_like = inmodel_xml[:-4]+"_results.xml" # this is the fit result
logfile_like = f'ctlike_emin={emin_bins:.2f}_emax={emax_bins:.2f}_side={side_skymap_x:.2f}_{spatial_type}.log'

cmd = f'python {scriptdir}/Like.py --inobs {outobs_bins} --expcube {outcube_exp} --psfcube {outcube_psf} --bkgcube {outcube_bkg} --inmodel {inmodel_like} --outmodel {outmodel_like} --logfile {logfile_like} --clobber'

print(cmd)
os.system(cmd)


