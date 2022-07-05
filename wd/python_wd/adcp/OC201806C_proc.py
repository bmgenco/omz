cruiseid = "OC201806C"  # for titles
yearbase = 2018  # usually year of first data logged
uhdas_dir = "/home/brandon/adcpproc/adcp_OC1806C/no_edits"



## for processing
##----------------
## ship name: shipname = ''
## at-sea "proc_cfg.py" initialized date = ''
##
## This file starts as /home/adcp/config/proc_cfg.py and
## includes the following information.  Uncomment, left-justify
## and fill these in if you are attempting to generate proc_cfg.py
## from this template.  The file must be named {cruiseid}_proc.py
## or for this example, kk1105_proc.py
##
## example values: fill in for your cruise...
#
# yearbase = 2011                  # usually year of first data logged
# uhdas_dir = '/home/data/kk1105'  # path to uhdas data directory
# shipname = 'Ka`imikai O Kanaloa' # for documentation
# cruiseid = 'kk1105'              # for titles
#
#

#======== serial inputs =========

# choose position instrument (directory and rbin message)

pos_inst = 'gpsnav'
pos_msg = 'gps'

# choose attitude  instruments (directory and rbin message)

pitch_inst = 'adu5'     # pitch is recorded, but NOT used in transformation
pitch_msg = 'adu'      # disable with '' (not None)

roll_inst = 'adu5'      # roll is recorded, but NOT used in transformation
roll_msg = 'adu'       # disable with '' (not None)

hdg_inst = 'gyro'       # reliable heading, used for beam-earth transformation
hdg_msg = 'hdg'


## heading correction
## all heading+msg pairs, for hbin files
hdg_inst_msgs = [
    ('gyro', 'hdg'),
    ('adu5', 'adu'),
    ('adu800', 'adu'),
]

## instrument for heading corr to ADCP data (dir and msg)
hcorr_inst = 'adu5'       # disable with '' (not None)
hcorr_msg = 'adu'        # disable with '' (not None)
hcorr_gap_fill = 3.0   ## fallback correction for hcorr gaps
                     ## calculate hdg_inst - hcorr_inst, eg gyro - ashtech
                     ## SAME SIGN CONVENTION as cal/rotate/ens_hcorr.ang

## if there is a posmv
acc_heading_cutoff = 0.02

# =========== ADCP transformations========
# historically, values were substituted into cruise_proc.m
# now, in Python, they are used directly

# heading alignment:  nominal - (cal/watertrack)
h_align = dict(
     wh300 = 51.9,
      os75 = -43.9,
)

# transducer depth, meters
ducer_depth = dict(
     wh300 = 5,
      os75 = 5,
)

# velocity scalefactor
# see SoundspeedFixer in pycurrents/adcp/pingavg.py
scalefactor = dict(
     wh300 = 1.0,
    os75bb = 1.0,
    os75nb = 1.0,
)

# soundspeed
# Soundspeed is usually None, and should ALWAYS be left as None for Ocean Surveyor
# (it is remotely possible that soundspeed for a WH, BB, or NB might need to
#           bet set to a nnumber, but usually that just results in an erroneous
#           scale factor.
soundspeed = dict(
     wh300 = None,
    os75bb = None,
    os75nb = None,
)

# salinity
salinity = dict(
     wh300 = None,
    os75bb = None,
    os75nb = None,
)

#=================================================================
# =========           values for quick_adcp.py          ==========
# ========= These are set here for at-sea procesing,    ==========
# ========= but are REQUIRED in quick_adcp.py control   ==========
# =========  file for batch mode or reprocessing.       ==========

## choose whether or not to use topography for editing
## 0 = "always use amplitude to guess the bottom;
##          flag data below the bottom as bad
## -1 = "never search for the bottom"
## positive integer: Only look for the bottom in deep water
##      "deep water" defined as "topo database says greater than this"

max_search_depth = dict(
     wh300 = 500,
    os75bb = 2000,
    os75nb = 2000,
)

# special: weakprof_numbins
weakprof_numbins = dict(
     wh300 = None,
    os75bb = None,
    os75nb = None,
)

# set averaging intervals
enslength = dict(
     wh300 = 120,
    os75bb = 300,
    os75nb = 300,
)

# Estimate of offset between ADCP transducer and gps:
# - Specify integer values for 'xducer_dx' and 'xducer_dy' for each instrument
# - 'xducer_dx' = ADCP's location in meters, positive starboard with the GPS
#   location as origin
# - 'xducer_dy' = ADCP's location in meters, positive forward with the GPS
#   location as origin
#
# There should be one set of 'xducer_dx'/'xducer_dy' values per instrument
# Ex.:
#   xducer_dx = dict(
#   wh300 = -2,
#   os38 = 16,
#   )
#   xducer_dy = dict(
#   wh300 = 1,
#   os38 = 6,
#   )
# Note that estimates of 'xducer_dx'/'xducer_dy' can be found in
# cal/watertrk/guess_xducerxy


xducer_dx = dict(
      os75 = 1,
     wh300 = 1,
)
xducer_dy = dict(
      os75 = -30,
     wh300 = -30,
)


## if there is a bad beam, create a dictionary modeled after
## enslen (i.e. Sonar-based, not instrument based) and use the
## RDI number (1,2,3,4) to designate the beam to leave out.


