#!/usr/bin/env python 
### This script runs the necessary post-processing of the MIDAS output 
### across all species in serial. 

import os
import sys
import config
import state_utils

if len(sys.argv) > 1:
    argument=sys.argv[1]
else:
    argument = 'all'

# Call create_StrainFinderInput for each species

os.system('python loop_over_species_wrapper.py %s python create_StrainFinderInput.py -o ${SF_DIR} --species' % argument)
