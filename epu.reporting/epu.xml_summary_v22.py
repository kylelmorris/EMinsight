#!/usr/bin/env python3

# These programs make up the package EMinsight
#
# Copyright (c) 2022
# Author(s): Kyle Morris, Jaehoon Cha, Dan Hatton
#
# xml parsing method follows that of Dan Hatton (DLS)
# https://github.com/d-j-hatton/python-smartem/blob/main/src/smartem/parsing/epu.py
#
# DBSCAN eps optimisation function developed by Jaehoon Cha (SciML STFC)
#
# Cryolo cbox file parser by Stephen Riggs (DLS)
#

# BSD 3-Clause License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# General modules
import os
import fnmatch
from os.path import normpath, basename, isfile
import pwd
import shutil
import sys
import glob
from glob import iglob
import time
import math
import datetime
import dateutil.parser
from tqdm import tqdm
tqdm.pandas()
import secrets
import re
import random
from pathlib import Path
import argparse
parser = argparse.ArgumentParser()

from itertools import groupby

# Plotting
import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.cm import get_cmap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from scipy import stats

# Data handling and analysis
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.preprocessing import StandardScaler

from time import time
from functools import reduce

from matplotlib.patches import Circle, Rectangle

# XML parsing
import json
from rich.pretty import pprint
from typing import Any, Dict, Tuple
import xmltodict

import xml.etree.ElementTree as ET

# Star and MRC file handling
from starparser import fileparser
import mrcfile

from gemmi import cif
import starfile # Alistair Burt's starfile
import pyem

#from cctbx.array_family import flex
#import iotbx.cif

# PDF and image handling
from fpdf import FPDF
from PIL import Image

# Suppress python warnings....
import warnings
warnings.filterwarnings("ignore")

# Deposition
import hashlib

# The following needs to be depreciated and only the dataframe used!
# Create dictionary of Krios calibrated magnified pixel sizes
# This will need to be cutomised for your facilities instruments
# https://stackabuse.com/python-dictionary-tutorial/
magTable = {'GLACIOS-9953831' : {'#name' : 'GlcsII', 'number' : '2', 'mXX' : 'm12', 'camera' : 'F4i', 'cameraWpx' : '4096', 'cameraHpx' : '4096', \
            'optimalDose' : '4.5', 'ID' : '9953831', 'PC' : 'GLACIOS-9953831', '120000' : '1.192', '150000' : '0.930'}, \
            'ARCTICA52203594' : {'#name' : 'Talos', 'number' : '1', 'mXX' : 'm04', 'camera' : 'K3', 'cameraWpx' : '5760', 'cameraHpx' : '4092', \
            'optimalDose' : '15', 'ID' : '3594', 'PC' : 'ARCTICA52203594', '63000' : '1.315', '79000' : '1.049', '100000' : '0.819', '130000' : '0.640', '165000' : '0.497'}, \
            'TITAN52334190' : {'#name' : 'Krios1', 'number' : '1', 'mXX' : 'm02', 'camera' : 'K3', 'cameraWpx' : '5760', 'cameraHpx' : '4092', \
            'optimalDose' : '15', 'ID' : '3419', 'PC' : 'TITAN52334190', '64000' : '1.350', '81000' : '1.072', '105000' : '0.831', '130000' : '0.653', '165000' : '0.513'}, \
            'TITAN52334150' : {'#name' : 'Krios2', 'number' : '2', 'mXX' : 'm03', 'camera' : 'K3', 'cameraWpx' : '5760', 'cameraHpx' : '4092', \
            'optimalDose' : '15', 'ID' : '3415', 'PC' : 'TITAN52334150', '64000' : '1.350', '81000' : '1.072', '105000' : '0.831', '130000' : '0.653', '165000' : '0.513'}, \
            'TITAN52335130' : {'#name' : 'Krios3', 'number' : '3', 'mXX' : 'm06', 'camera' : 'F4i', 'cameraWpx' : '4096', 'cameraHpx' : '4096', \
            'optimalDose' : '7', 'ID' : '3513', 'PC' : 'TITAN52335130', '53000' : '2.325', '64000' : '1.9', '81000' : '1.501', '85000' : '1.501', '105000' : '1.171', '130000' : '0.921', '165000' : '0.723', '215000' : '0.576'}, \
            'TITAN52336320' : {'#name' : 'Krios4', 'number' : '4', 'mXX' : 'm07', 'camera' : 'K3', 'cameraWpx' : '5760', 'cameraHpx' : '4092', \
            'optimalDose' : '15', 'ID' : '3632', 'PC' : 'TITAN52336320', '53000' : '1.630', '64000' : '1.340', '81000' : '1.060', '105000' : '0.829', '130000' : '0.651', '165000' : '0.508'}}

# Create a dataframe of Krios calibrated magnified pixel sizes and equipment parameters
# Optmal dose for Falcon 4 is 3-6 e/px/s, for the Falcon 4i is 7 e/px/s
magTableData = [{'ID':'GLACIOS-9953831', '#name' : 'GlcsII', 'number' : '2', 'mXX' : 'm12', 'camera' : 'F4i', 'cameraWpx' : '4096', 'cameraHpx' : '4096', \
            'optimalDose' : '4.5', 'serial' : '9953831', 'PC' : 'GLACIOS-9953831', '120000' : '1.192', '150000' : '0.930'}, \
            {'ID':'TITAN52334190', '#name' : 'Krios1', 'number' : '1', 'mXX' : 'm02', 'camera' : 'K3', 'cameraWpx' : '5760', 'cameraHpx' : '4092', \
            'optimalDose' : '15', 'serial' : '3419', 'PC' : 'TITAN52334190', '64000' : '1.350', '81000' : '1.072', '105000' : '0.831', '130000' : '0.653', '165000' : '0.513'}, \
            {'ID':'TITAN52334150', '#name' : 'Krios2', 'number' : '2', 'mXX' : 'm03', 'camera' : 'K3', 'cameraWpx' : '5760', 'cameraHpx' : '4092', \
            'optimalDose' : '15', 'serial' : '3415', 'PC' : 'TITAN52334150', '64000' : '1.350', '81000' : '1.072', '105000' : '0.831', '130000' : '0.653', '165000' : '0.513'}, \
            {'ID':'ARCTICA52203594', '#name' : 'Talos', 'number' : '1', 'mXX' : 'm04', 'camera' : 'K3', 'cameraWpx' : '5760', 'cameraHpx' : '4092', \
            'optimalDose' : '15', 'serial' : '3594', 'PC' : 'ARCTICA52203594', '63000' : '1.315', '79000' : '1.049', '100000' : '0.819', '130000' : '0.640', '165000' : '0.497'}, \
            {'ID':'TITAN52335130', '#name' : 'Krios3', 'number' : '3', 'mXX' : 'm06', 'camera' : 'F4i', 'cameraWpx' : '4096', 'cameraHpx' : '4096', \
            'optimalDose' : '7', 'serial' : '3513', 'PC' : 'TITAN52335130', '53000' : '2.325', '64000' : '1.9', '81000' : '1.9', '85000' : '1.501', '105000' : '1.171', '130000' : '0.921', '165000' : '0.723', '215000' : '0.576'}, \
            {'ID':'TITAN52336320', '#name' : 'Krios4', 'number' : '4', 'mXX' : 'm07', 'camera' : 'K3', 'cameraWpx' : '5760', 'cameraHpx' : '4092', \
            'optimalDose' : '15', 'serial' : '3632', 'PC' : 'TITAN52336320', '53000' : '1.630', '64000' : '1.340', '81000' : '1.060', '105000' : '0.829', '130000' : '0.651', '165000' : '0.508'}]
magDfTable = pd.DataFrame(magTableData)

bagTableData = [{'ID':'bi23268', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'2'},
            {'ID':'bi31336', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'2'},
            {'ID':'bi22238', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'3'},
            {'ID':'bi31589', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'3'},
            {'ID':'bi28713', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'4'},
            {'ID':'bi26703', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'5'},
            {'ID':'bi28654', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'6'},
            {'ID':'bi25127', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'7'},
            {'ID':'bi21404', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'8'},
            {'ID':'bi27980', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'8'},
            {'ID':'bi22006', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'9'},
            {'ID':'bi30374', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'9'},
            {'ID':'bi24557', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'10'},
            {'ID':'bi31827', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'10'},
            {'ID':'bi25452', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'11'},
            {'ID':'bi26876', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'12'},
            {'ID':'bi23872', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'13'},
            {'ID':'bi31586', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'13'},
            {'ID':'bi21809', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'14'},
            {'ID':'bi28549', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'14'},
            {'ID':'bi22724', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'15'},
            {'ID':'bi29255', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'15'},
            {'ID':'bi28576', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'16'},
            {'ID':'bi25222', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'17'},
            {'ID':'bi24039', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'18'},
            {'ID':'bi25832', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'19'},
            {'ID':'bi23047', 'name':'anonymous', 'institute':'anonymous', 'anonymous':'20'}]
bagTable = pd.DataFrame(bagTableData)

# Argument parsing
# https://stackoverflow.com/questions/11604653/how-to-add-command-line-arguments-with-flags-in-python3
parser.add_argument("-i", "--input", help="Input directory, will search for EPU sessions, and Atlas unless specified by --atlas")
parser.add_argument("-e", "--epu", help="Depreciated: Custom: Manual entry of EPU session top level directory")
parser.add_argument("-a", "--atlas", help="Optional: Force use of an Atlas session, top level directory or ScreeningSession.dm")
parser.add_argument("-d", "--dose", help="Recommended: Dose rate over vacuum (e-/px/s) (default: estimated based on detector)")
parser.add_argument("-o", "--output", help="Optional: Directory for outputs (default: EPU session directory)")
parser.add_argument("-r", "--report", help="Optional: Y = Generate a PDF report (default: Y)")
parser.add_argument("-cc", "--corcount", help="Optional: Y = Count motion correction mics and timestamps - slow (default: N)")
parser.add_argument("-t", "--timemethod", help="Optional: X = XML interal timestamp (slow/accurate), F = file timestamp (fast/inaccurate), M = MRC header timestamp (default = X)")
parser.add_argument("-apix", "--apix", help="Optional: Microscope magnified A pixel size (default: automatically detected)")
parser.add_argument("-p", "--print", help="Optional: Y = Only print xml and exit")
parser.add_argument("--screening", help="Optional: Y = Create screening report PDF (square analysis) - (default: N)")
parser.add_argument("--silent", help="Optional: Y = Print less to the terminal")
parser.add_argument("--doppio", help="Optional: Y = Only run minimal metadata capture for CCPEM doppio deposition")
parser.add_argument("--sizecalc", help="Optional: Y = Calculate sizes of Supervisor/raw/processed directories")
args = parser.parse_args()
print( "inputs {}".format(
        args.input,
        args.epu,
        args.atlas,
        args.dose,
        args.output,
        args.apix,
        args.report,
        args.timemethod,
        args.print,
        args.screening,
        args.silent,
        args.doppio
        ))

# Main processes the input command but run calls the functions that do the work
def main():

    # Get absolute dirpath of this script as it is run
    # https://note.nkmk.me/en/python-script-file-path/
    #script = os.path.abspath(__file__)
    main.scriptDir = os.path.dirname(os.path.abspath(__file__))
    main.noAtlas = main.scriptDir+'/slot_empty.png'

    # Default behaviour for finding timestamps of exposures is accurate method
    if not args.timemethod:
        args.timemethod = 'X'

    # Default behaviour for generating report
    if not args.report:
        args.report = 'Y'

    # Default behaviour for counting motion corrected micrographs and checking timestamps
    if not args.corcount:
        args.corcount = 'N'

    # Default behaviour for gridsquare screening report
    if not args.screening:
        args.screening = 'N'

    # needs to be seperate function so can be used to reset directories when looping over run function
    defineOutputFolders()
    report_dir = Path(main.report_dir)
    plot_dir = Path(main.plot_dir)
    csv_dir = Path(main.csv_dir)
    pdf_dir = Path(main.pdf_dir)
    dep_dir = Path(main.dep_dir)

    ## Main inputs
    # Are we looking in a directory automatically or going for direct inputs
    # The following is all a bit hideous but works in tests
    # Principle development done assuming automatic Supervisor directory searching is preferred method
    # Have tested manual input of directories and this is working 30/08/22
    if not args.input:
        print('Top level directory input not provided: not performing search for EPU and Atlas directories.')
        if not args.epu:
            print("EPU Session directory also not supplied, don't know where to look for EPU sessions.")
            print("Check your script inputs, exiting...")
            exit(1)
        else:
            # Look for EpuSession.dm in defined top level EPU session directory
            path_to_file = str(args.epu+'/EpuSession.dm')
            path = Path(path_to_file)

            directory = os.path.realpath(args.epu)
            realPath = Path(directory).parents[0]
            p = str(realPath)
            realDLSName = basename(normpath(p))
            main.visitNameDLS = realDLSName

            if path.is_file():
                print(f'The file {path_to_file} exists')
                # This is the path to the epusession xml
                main.xml = path_to_file
                # This is the path to the visit directory
                main.visitPath = realPath
                # This is the top level EPU session directory
                main.epu = args.epu
                main.epuSessionPath = main.epu
                # Define atlas from specific input
                main.atlas = args.atlas
                main.atlasSessionPath = main.atlas
            else:
                print(f'The file {path_to_file} does not exist, exiting')
                exit(1)
            # This will parse the xml in the target directory, wow this is getting messy
            run()
    else:
        #################################################
        ## This is the main intended way to use EMinsight
        #################################################
        # In here you need to do directory searching

        # Input directory to search for data directories
        main.dirInput = args.input

        # Location of EPU session directory on which this script was ran
        realPath = os.path.realpath(args.input)

        # The user may provide input straight to an EpuSession.dm
        # We then don't know where to look for other directories
        # Loop up through directories to find parent visit directory
        # DEV DEV DEV THIS ASSUMES /Processing IS FOUND IN TOP LEVEL DIRECTORY, BIG ASSUMPTION FOR OTHER FACILITIES
        # CAREFUL: INFINITE LOOPS HAVE HAPPENED HERE
        i = 0
        while i == 0:
            if os.path.isdir(realPath+'/processing'):
                print('Visit parent directory found: '+str(realPath))
                i = 1
                continue
            else:
                # Move up a level
                realPath = os.path.dirname(realPath)
                # Deal with the case where you move up to top level and don't find anything
                if str(realPath) == '/':
                    print('No visit found in the providing path, skipping')
                    exit(1)

        # print(realPath)
        # parents will look up a level if required, 0 is top level
        # p = Path(str(realPath)).parents[0]
        # This will just get the real path and the last directory of that path
        p = str(realPath)
        realDLSName = basename(normpath(p))
        main.visitNameDLS = realDLSName
        main.sessionTypeDLS = realDLSName[:2]

        report_dir = Path(main.report_dir)
        print('Output for report PDF, csv and plots will be: '+str(report_dir))

        # Make report directories
        # Only make report directory, if doesn't exist
        if report_dir.exists() and report_dir.is_dir():
            print('Output directory exists')
        else:
            os.makedirs(report_dir)

        if dep_dir.exists() and dep_dir.is_dir():
            print('Dep directory exists')
        else:
            os.makedirs(dep_dir)

        # If doppio flag declared, you don't need other directories
        if not args.doppio:
            if csv_dir.exists() and csv_dir.is_dir():
                print('csv directory exists')
            else:
                os.makedirs(csv_dir)

            if plot_dir.exists() and plot_dir.is_dir():
                print('plot directory exists')
            else:
                os.makedirs(plot_dir)

            if pdf_dir.exists() and pdf_dir.is_dir():
                print('pdf directory exists')
            else:
                os.makedirs(pdf_dir)

        ## Use function to get all Supervisor directories
        # Search supervisor returns a dataframe with column names [Directory, Epu, Atlas] and thus what their type is
        #This was the old way to find Supervisor directories
        #supervisorList = searchSupervisor(main.dirInput)
        # Now searches directly by looking for top level metadata file
        supervisorList = searchMetadata(main.dirInput)

        supervisorListNo = supervisorList.shape[0]

        # Get directories that are not an atlas, these will be analysed
        directoryList = supervisorList.loc[supervisorList['Atlas'] == 'None']
        directoryListNo = directoryList.shape[0]

        # Now get the atlas and epu directories from the Supervisor dataframe
        # https://stackoverflow.com/questions/17071871/how-do-i-select-rows-from-a-dataframe-based-on-column-values
        atlasList = supervisorList.loc[supervisorList['Atlas'] == 'ScreeningSession.dm']
        atlasListNo = atlasList.shape[0]

        if atlasList.shape[0] > 1:
            print("Warning, multiple Atlas directories, I will use the most recent, hope that's right")
            print()
        # If you need an epuList this is it
        epuList = supervisorList.loc[supervisorList['Epu'] == 'EpuSession.dm']
        epuListNo = epuList.shape[0]

        if supervisorList.empty:
           print('No Supervisor directories found...')
           print('Falling back to count movies and write this out per raw directory...')
           lookupBAG(main.visitNameDLS)
           date = fileTime = datetime.datetime.fromtimestamp(os.path.getmtime(realPath))

           # Find out what microscope this was on from the path
           scope = getScopeFromPath(realPath)

           # Start looking for raw directories
           rawdir = realPath+'/raw*'
           # This will return absolute paths for all raw directories
           rawlist = [f for f in iglob(rawdir, recursive=True) if os.path.isdir(f)]
           # If any raw directory exists
           if rawlist:
              # Loop through raw directories and find movies, reliant on them having standard epu data output name syntax
              for d in rawlist:
                 print('Found: '+d)
                 print()
                 # Count number of raw movies
                 #rawNo = countMovies(realPath, d)
                 rawMovies = glob.glob(d+'/**/FoilHole_*_Data_*', recursive = True)
                 rawMovieName = os.path.splitext(os.path.split(rawMovies[0])[1])[0]
                 tmp1 = rawMovieName.replace('_fractions','')
                 name = tmp1.replace('_EER','')
                 rawNo = len(rawMovies)
                 rawName = os.path.basename(os.path.normpath(d))
                 sessionName = realDLSName+'-'+rawName
                 print('Counted movies in '+sessionName+': '+str(rawNo))
                 print()

                 # DEV DEV DEV Attempt to search raw for movies to get times and then do processed search, not working right now
                 #search_mics('F', d, 'ractions')

                 # See if you can find the processed directory
                 print('Attempting minimal processed check...')
                 print(name)
                 searchProcessed(name, realPath, 'N')

                 print(str(d))
                 # Send the results to standard error output
                 standard_error(sessionName,realPath,date,'Super',d,rawNo,scope,'Unknown')
           else:
               print('No Raw or Supervisor, visit is recorded as Empty')
               print('')
               standard_error(main.visitNameDLS,realPath,date,'Empty','0','0',scope,'Unknown')
        else:
           # Prints out the supervisory directory list, including classification as EPU, Tomo or Atlas
           print(supervisorList[['Supervisor', 'Epu', 'Tomo', 'Atlas']])

        # Now go through each Supervisor directory that is not an atlas in a Loop
        j = 1
        #for index, row in epuList.iterrows():
        for index, row in directoryList.iterrows():
            dir = row.Directory
            sessionName = os.path.basename(os.path.normpath(dir))
            # Date is a best guess based on directory date
            date = fileTime = datetime.datetime.fromtimestamp(os.path.getmtime(dir))
            print('')
            print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
            print('Working on Supervisor directory '+str(j)+'/'+str(directoryListNo)+': '+dir)
            print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

            # Write out a text file for keeping track of script success/failure
            # Can later append to this file to record where script may fail on certain session directories
            #epuDir = basename(dir)
            #with open(main.pdf_dir+'/'+str(epuDir)+'_script.out', 'w') as f:
            #    f.write('epu.xml_summary.py ran on: '+dir)

            # Find out what microscope this was on from the path in case you need to bail out into standard error
            scope = getScopeFromPath(realPath)

            # If no EpuSession.dm found, something wrong or it's tomo, iterate to next directory
            #if row.Epu == 'None' and row.Tomo == 'None' and row.Atlas == 'None':
            #    print('Failed to identify Supervisor directory type, skipping')
            #    print()
            #    standard_error(sessionName,realPath,date,'Error','0','0')
            #    j = j + 1
            #    continue
            if row.Epu == 'None' and row.Tomo == 'None':
                # Get the number of mics
                # Search for dir using findRaw
                # Then count mics using countMovies
                print(dir)
                print(dir+': No EpuSession.dm found, falling back to count mics in epu directory')
                print()
                searchedFiles = glob.glob(dir+"/**/Data/*xml", recursive=True)
                if not searchedFiles:
                   continue
                oneFile = searchedFiles[0]
                name = str(os.path.splitext(os.path.basename(oneFile))[0])
                raw = findRaw(name, dir+'/..')[0]
                # Count movies in appropriate raw folder
                rawNo = countMovies(realPath, raw)

                # See if you can find the processed directory
                print('Attempting minimal processed check...')
                print(name)
                searchProcessed(name, realPath, 'N')

                # This function will attempt to write out what the sirectory is if not EPU
                standard_error(sessionName,realPath,date,'EPU_dm',raw,rawNo,scope,'SPA')
                #standard_error(sessionName,realPath,date,'EPU_dm','0','0',scope)
                j = j + 1
                continue
            elif row.Epu == 'None' and row.Tomo == 'Batch':
                print(dir+': Tomo found, skipping')
                print()
                # This function will attempt to write out what the sirectory is if not EPU
                standard_error(sessionName,realPath,date,'Tomo','0',np.nan,scope,'Tomo')
                j = j + 1
                continue
            else:
                print(dir+': EpuSession.dm found')
                print()
                epuDir = dir

            # If no Atlas found, something is wrong, iterate to the next directory
            if atlasList.empty:
                # Get the number of mics
                # Search for dir using findRaw
                # Then count mics using countMovies
                print(dir)
                print(dir+': No Atlas found, falling back to count mics in epu directory')
                print()
                searchedFiles = glob.glob(dir+"/**/Data/*xml", recursive=True)
                if not searchedFiles:
                   print('Nothing found')
                   continue
                oneFile = searchedFiles[0]
                name = str(os.path.splitext(os.path.basename(oneFile))[0])
                raw = findRaw(name, dir+'/..')[0]
                # Count movies in appropriate raw folder
                rawNo = countMovies(realPath, raw)

                # See if you can find the processed directory
                print('Attempting minimal processed check...')
                print(name)
                searchProcessed(name, realPath, 'N')

                print('No Atlas found, how do you expect me to work under these conditions, skipping')
                print('')
                standard_error(sessionName,realPath,date,'Atl_dm',raw,rawNo,scope,'SPA')
                j = j + 1
                continue
            else:
                #First atlas
                #atlasDir = atlasList.iloc[0]['Directory']
                #Last atlas, i.e. most recent
                atlasDir = atlasList.iloc[-1]['Directory']

            main.epu = epuDir
            main.atlas = atlasDir
            main.visitPath = realPath
            main.epuSessionPath = main.epu
            main.atlasSessionPath = main.atlas

            # Look for EpuSession.dm in defined top level EPU session directory
            path_to_file = str(main.epu+'/EpuSession.dm')
            path = Path(path_to_file)

            if path.is_file():
                print(f'The file {path_to_file} exists')
                # This is the path to the epusession xml
                main.xml = path_to_file
            else:
                print(f'The file {path_to_file} does not exist, exiting')
                exit(1)

            print()
            print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
            print()

            # Before running full eminsight analysis, look for all image files, via xml, mrc or jpg
            # This will declare global searchSupervisorData variables with the file lists, for xml, mrc and jpg
            searchSupervisorAtlas(atlasDir)
            searchSupervisorData(epuDir)

            # Get presets for EPU session xml
            print('')
            print('\033[1m'+'Finding all presets from EPU session:'+'\033[0m')
            print('')

            xml_presets(main.xml)

            # Get main set up parameters from EPU session xml
            print('')
            print('\033[1m'+'Finding main EPU session parameters:'+'\033[0m')
            print('')

            xml_session(main.xml)

            # Find mics via xml for counting
            searchedFiles = find_mics(main.epu, 'xml')
            if searchedFiles == 'exit':
                continue

            count(searchedFiles)
            # Rough and ready data acquisition metadata pull, unsorted mic list
            randomMic = searchedFiles[0]
            xml_presets_data(randomMic)

            # Before running full eminsight analysis run get a basic look at all the essential metadata you would need for an EMDB deposition
            # DEV DEV DEV take variables from xml parsing and put into the deposition file, like defocus range etc
            doppio(main.xml)

            # If doppio flag declared, skip doing full analysis by exiting loop
            if args.doppio == "Y":
                j = j + 1
                continue
            else:
                # Run full eminsight analysis
                run()
                j = j + 1

    print('Finished')

def defineOutputFolders():
    if not args.output:
        # Default locations for file outputs is up one level from EPU session directory
        main.report_dir = args.input+"/EMinsight"
        main.plot_dir = main.report_dir+'/plots'
        main.csv_dir = main.report_dir+'/csv'
        main.pdf_dir = main.report_dir+'/report'
        main.dep_dir = main.report_dir+'/dep'

    else:
        main.report_dir = args.output+"/EMinsight"
        main.plot_dir = main.report_dir+'/plots'
        main.csv_dir = main.report_dir+'/csv'
        main.pdf_dir = main.report_dir+'/report'
        main.dep_dir = main.report_dir+'/dep'

def doppio(xml):

    # Get EPU session name from main EPU xml file, this is a function
    sessionName = xml_sessionName(xml)

    # This is the data xml metadata file already in a dictionary
    data = searchSupervisorData.xmlDataDict["MicroscopeImage"]

    # Get mag
    xmlMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
    xmlMetrePix = data["SpatialScale"]["pixelSize"]["x"]["numericValue"]
    xmlAPix = float(xmlMetrePix) * 1e10
    xmlAPix = roundup(xmlAPix, 1)

    # Get scope and kV
    model = data["microscopeData"]["instrument"]["InstrumentModel"]
    eV = data["microscopeData"]["gun"]["AccelerationVoltage"]

    # Save doppio deposition csv file
    dictHorizontal = [
    {'Microscope': model,
     'epuversion': xml_session.epuVersion,
     'date': xml_session.sessionDate,
     'eV': eV,
     'mag': xmlMag,
     'apix': str(xmlAPix),
     'nominal_defocus_min_microns': xml_session.defocusMin,
     'nominal_defocus_max_microns': xml_session.defocusMax,
     'spot_size': xml_presets.spot,
     'C2_micron': xml_presets.C2,
     'Objective_micron': str(xml_presets_data.objective),
     'Beam_diameter_micron': xml_presets.beamD,
     'collection': xml_session.afisMode,
     'number_of_images': count.mic_count
     }
    ]
    df = pd.DataFrame.from_dict(dictHorizontal)

    ## Deposition file
    depfilepath = main.dep_dir+'/'+sessionName+'_dep.json'
    checksumpath = main.dep_dir+'/'+sessionName+'_dep.checksum'

    # Human readable deposition file
    #df.to_csv (main.dep_dir+'/'+sessionName+'.dep', index = False, header=True)
    df_transpose = df.T
    df_transpose.to_csv(main.dep_dir+'/'+sessionName+'_dep.csv', index = True, header=True)

    # Computer readable deposition file
    df_transpose.to_json(depfilepath, index = True)

    # This can be run before doing full analysis of the session directories
    print("Created full deposition metadata file")

    # Deposition Checksum
    checksum(depfilepath, checksumpath)

def doppio_full(xml):

    # Get EPU session name from main EPU xml file, this is a function
    sessionName = xml_sessionName(xml)

    # This is the data xml metadata file already in a dictionary
    data = searchSupervisorData.xmlDataDict["MicroscopeImage"]

    # Get mag
    xmlMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
    xmlMetrePix = data["SpatialScale"]["pixelSize"]["x"]["numericValue"]
    xmlAPix = float(xmlMetrePix) * 1e10
    xmlAPix = roundup(xmlAPix, 1)

    # Get scope and kV
    model = data["microscopeData"]["instrument"]["InstrumentModel"]
    eV = data["microscopeData"]["gun"]["AccelerationVoltage"]

    # Save doppio deposition csv file
    dictHorizontal = [
    {'Microscope': model,
     'epuversion': xml_session.epuVersion,
     'date': xml_session.sessionDate,
     'eV': eV,
     'mag': xmlMag,
     'apix': str(xmlAPix),
     'nominal_defocus_min_microns': xml_session.defocusMin,
     'nominal_defocus_max_microns': xml_session.defocusMax,
     'spot_size': xml_presets.spot,
     'C2_micron': xml_presets.C2,
     'Objective_micron': str(xml_presets_data.objective),
     'beam_diameter_micron': xml_presets.beamD,
     'collection': xml_session.afisMode,
     'number_of_images': count.mic_count,
     'grid_type': xml_session.gridType,
     'available_squares': xml_targets.totalSquares,
     'collected_squares': xml_targets.collectedSquares,
     'average_foils_per_square': xml_targets.avSqFoilHoleNo,
     'hole_size_micron': xml_targets.holeSizeRegular,
     'hole_space_micron': xml_targets.holeSpaceRegular,
     'shots_per_hole': xml_template_dataframe.shotsPerHole,
     'total_dose_eA2': calc_params.doseTotal,
     'fraction_dose_eA2': calc_params.doseFrameTotal,
     }
    ]
    df = pd.DataFrame.from_dict(dictHorizontal)

    ## Deposition file
    depfilepath = main.dep_dir+'/'+sessionName+'_dep.json'
    checksumpath = main.dep_dir+'/'+sessionName+'_dep.checksum'

    # Human readable deposition file
    #df.to_csv (main.dep_dir+'/'+sessionName+'.dep', index = False, header=True)
    df_transpose = df.T
    df_transpose.to_csv(main.dep_dir+'/'+sessionName+'_dep.csv', index = True, header=True)

    # Computer readable deposition file
    df_transpose.to_json(depfilepath, index = True)

    # This can be run before doing full analysis of the session directories
    print("Created full deposition metadata file")

    # Deposition Checksum
    checksum(depfilepath, checksumpath)

def checksum(path, out):
    # https://www.quickprogrammingtips.com/python/how-to-calculate-sha256-hash-of-a-file-in-python.html
    filename = path
    sha256_hash = hashlib.sha256()
    with open(filename,"rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096),b""):
            sha256_hash.update(byte_block)
        checksum = sha256_hash.hexdigest()

    # Open file for writing
    checkfile = open(out, "w")
    # Write checksum string
    n = checkfile.write(checksum)
    # Close file
    checkfile.close()

    print('Created checksum')
    print()

def run():

    # DEV DEV DEV
    # Look in here, can analysis be seperated from report generation for running wihtout doing full analysis??

    #
    main.apix = args.apix
    main.printXml = args.print

    main.output = args.output

    # Get name of user who ran the script
    main.user = get_username()

    # datetime for when script was run and report generated
    now = datetime.datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S") # Formatted if you require
    # Report
    main.runDate = dt_string

    # Run major functions
    # DEV There are major running order dependancies here, need to tidy up

    if not main.atlas:
        print('No atlas directory provided')
    else:
        atlas_xml(main.atlas)
        screeningSession_xml(main.atlas)

    # Now we know about the session we can define output names
    # Reset the directory outputs before modification
    defineOutputFolders()

    # Modify plots output and make directory
    main.plot_dir = main.plot_dir+'/'+xml_session.sessionName
    try:
        os.mkdir(main.plot_dir)
    except OSError:
        pass

    # Modify PDF output and make directory
    main.pdf_dir = main.pdf_dir+'/'+xml_session.sessionName
    try:
        os.mkdir(main.pdf_dir)
    except OSError:
        pass

    # Output paths and names
    main.report_name = xml_session.sessionName+'_session.pdf'
    main.report_path = main.pdf_dir+'/'+main.report_name

    # Output paths and names
    main.report_processed_name = xml_session.sessionName+'_processed.pdf'
    main.report_processed_path = main.pdf_dir+'/'+main.report_processed_name

    # Output paths and names
    main.report_screening_name = xml_session.sessionName+'_screening.pdf'
    main.report_screening_path = main.pdf_dir+'/'+main.report_screening_name

    # Find BAG
    name, institute, proposal, anonymous = lookupBAG(main.visitNameDLS)

    # Report to terminal
    if not args.silent or args.silent == 'N':
       print('Atlas screening session name: '+str(screeningSession_xml.name))
       print('Screening session start time: '+str(screeningSession_xml.date))
       print('')
       print('Autoloader position: '+str(xml_session.autoSlot))
       print('Autoloader position name: '+str(atlas_xml.atlasName))
       print('')
       print('Atlas xml: '+main.atlasxml)
       print('Original atlas location: '+xml_session.atlasDirOrig)
       print('')
       print('EPU session name: '+xml_session.sessionName)
       print('EPU session created: '+str(xml_session.sessionDate))
       print('')
       print('Clustering mode: '+xml_session.clustering)
       print('Clustering radius: '+str(xml_session.clusteringRadius)+' microns')
       print('')
       print('Focus with: '+xml_session.focusWith)
       print('Focus recurrence: '+xml_session.focusRecurrence)
       print('')
       print('Delay after Image Shift: '+xml_session.delayImageShift+' seconds')
       print('Delay after Stage Shift: '+xml_session.delayStageShift+' seconds')
       print('')
       print('Defocus list (um): ')
       print(xml_session.defocus)
       print('')
       print('BAG name: ')
       print(name)
       print('BAG institute: ')
       print(institute)
       print('BAG proposal: ')
       print(proposal)
       print()

    # Get Supervisor directory size
    if args.sizecalc == 'Y':
        print('Calculating Supervisor directory disk space usage: '+str(main.epu))
        size = get_dir_size(main.epu, 'GB')
        searchSupervisor.size = size
        print(str(size)+' GB')
        print()
    else:
        searchSupervisor.size = np.nan

    # Search for micrographs and report their collection timings, method is args.timemethod X (xml time record), F (file timestamp), M (MRC file header time record)
    # look in epu directory and for xml files
    searchedFiles = find_mics(main.epu, 'xml')
    search_mics(searchedFiles, args.timemethod)
    print('')

    # This will establish the microscope name, but also the angstrom per pixel
    if not main.apix:
        getScopeNameMag(search_mics.firstXml)
    else:
        getScopeName(search_mics.firstXml)

    print('')
    print('\033[1m'+'Finding target statistics:'+'\033[0m')
    print('')
    xml_targets(main.xml)

    xml_presets_data(search_mics.firstXml)

    # Report to terminal
    print('Number of available squares: '+str(xml_targets.totalSquares))
    print('Number of targeted squares: '+str(xml_targets.targetSquares))
    print('Number of collected squares: '+str(xml_targets.collectedSquares))

    print('')
    print('Hole size: '+str(xml_targets.holeSizeRegular)+' microns regularised ('+str(xml_targets.holeSize)+' microns user measured)')
    print('Hole spacing: '+str(xml_targets.holeSpaceRegular)+' microns regularised ('+str(xml_targets.holeSpace)+' microns user measured)')
    print('')

    # Get the shot template coordinates
    xml_template_dataframe(main.xml)
    #Plot the shot template coordinates
    xml_template_dataframe_plot(xml_template_dataframe.dfShot)
    print('')

    #
    calc_params()
    print('')

    # This is dependent on search_mics having been run to get the list of files
    count_mics(search_mics.micList)
    calcTimings()

    # This is not dependent on search_mics
    count_movies(search_mics.micList, main.visitPath)

    # Get size of visit
    #main.size = os.stat(main.visitPath).st_blocks * 512
    #print("Visit directory size: " + str(main.size)+" bytes")

    # Run search processed pipeline routines
    #print(main.visitPath)
    #print(search_mics.lastName)
    #input()
    searchProcessed(search_mics.lastName, main.visitPath, 'Y')

    # Intercept and look for anything strange about the session
    assessSession()
    writeAdvisory()

    #Writwe out data structure
    verifyDataStructure(search_mics.sessionStructure)
    #input('Press Enter to continue')

    ## Create OUTPUTS OUTPUTS OUTPUTS
    print('Writing EPU supervisor directory data structure to csv...')
    print()
    #print(search_mics.sessionStructure)
    search_mics.sessionStructure.to_csv(main.csv_dir+'/'+xml_session.sessionName+'_data_structure.csv', index = False, header=True)

    print('Creating full deposition metadata file')
    print()
    doppio_full(main.xml)

    standard_output('SPA')
    if searchProcessed.ispyb == 'Y':
       processed_output()

    # Session PDF report
    if not args.report:
        print('Not generating PDF session report, will just create standard csv output.')
    else:
        # Create report
        construct_session_report(main.report_path)

    # Processed PDF report
    if not args.report:
        print('Not generating PDF processing report, will just create standard csv output.')
    else:
        if searchProcessed.processedRan == 'Y':
            # Create report
            construct_process_report(main.report_processed_path)
        else:
            print('Not generating PDF processing report, no processed results found.')

    ## Get details of all collected squares and images for screening reporting
    # Relies on data from processed pipeline
    if args.screening == "Y":
        gatherScreening(xml_targets.listSquares, main.report_screening_path)

def getScopeFromPath(path):
    # Expects path /dls/mXX/data/2023/bi23047-106 at DLS but code should handle when mXX is at a different level
    # When it cannot find mXX as might be the case in other facilities the function returns m00

    # Regular expression pattern for finding 'mXX'
    pattern = re.compile(r'/m\d{2}/')

    # Search for the pattern in the path
    match = pattern.search(path)
    if match:
        # Extracting 'mXX' from the match
        mXX = match.group(0).strip('/')
        # Assuming magDfTable is a pre-defined DataFrame with necessary data
        scope = magDfTable.loc[magDfTable['mXX'] == mXX, '#name'].item()
    else:
        # Handle the case where 'mXX' is not found
        #raise ValueError("mXX not found in path")
        scope = 'm00'

    return scope

def lookupBAG(visitID):
    proposal = visitID.split('-')[0]
    #print(proposal)
    #print(bagTable)

    if 'nt' in str(proposal):
       bagName = 'inhouse (nt)'
       bagInstitute = 'eBIC-ext'
       bagAnonymous = '0'
    elif 'nr' in str(proposal):
       bagName = 'inhouse (nr)'
       bagInstitute = 'eBIC'
       bagAnonymous = '0'
    elif 'cm' in str(proposal):
       bagName = 'comission'
       bagInstitute = 'eBIC'
       bagAnonymous = '0'
    else:
       bag = bagTable.loc[bagTable['ID'].str.contains(proposal, case=False)]
       if bag.empty:
           bagName = 'notBAG'
           bagInstitute = 'notBAG'
           bagAnonymous = '1'
       else:
           bagName = bag.to_csv(columns=['name'], sep='\t', index=False, header=False)
           #print(bag(index=False, columns=['name']))
           bagInstitute = bag.to_csv(columns=['institute'], sep='\t', index=False, header=False)
           #print(bag(index=False, columns=['institute']))
           bagAnonymous = bag.to_csv(columns=['anonymous'], sep='\t', index=False, header=False)

    lookupBAG.name = bagName.replace("\n", "")
    lookupBAG.institute = bagInstitute.replace("\n", "")
    lookupBAG.anonymous = bagAnonymous.replace("\n", "")
    lookupBAG.proposal = proposal
    return bagName, bagInstitute, proposal, bagAnonymous

def searchMetadata(path):
    ## EpuSession.dm may be contained in top level Supervisor* directory, or inside raw* data folders (murfey)
    # This will search for SPA collection session files
    #supervisorEpuList = glob.glob(path+'/[Ss]upervisor*/EpuSession.dm') + glob.glob(path+'/raw*/metadata/EpuSession.dm')
    if "EpuSession.dm" in str(path): # If path is direct to EpuSession.dm just confirm its found
        supervisorEpuList = glob.glob(path)
    else: # If path to directory, search directories
        supervisorEpuList = glob.glob(path+'/**/EpuSession.dm') + glob.glob(path+'/EpuSession.dm') + glob.glob(path+'/[Ss]upervisor*/EpuSession.dm') + glob.glob(path+'/raw*/[Mm]etadata*/EpuSession.dm')
    # This will search for Tomo collection session files
    # Are you sure tomography session will always contain this Position_1.mdoc?
    supervisorTomoList = glob.glob(path+'/[Ss]upervisor*/Position_1.mdoc')
    # Then concatonate Data Supervisor directory lists together
    supervisorDataList = supervisorEpuList + supervisorTomoList
    # Perform filtering on list if the user has asked only to analyse in a particular directory
    #list(filter(lambda k: 'ab' in k, lst))

    if not args.atlas:
        # This will search for Atlas screening session files
        supervisorAtlasList = glob.glob(path+'/[Ss]upervisor*/ScreeningSession.dm') + glob.glob(path+'/atlas/[Ss]upervisor*/ScreeningSession.dm')
    elif "ScreeningSession.dm" in str(args.atlas):
        supervisorAtlasList = glob.glob(args.atlas)
    else:
        supervisorAtlasList = glob.glob(args.atlas+'/[Ss]upervisor*/ScreeningSession.dm') + glob.glob(args.atlas+'/atlas/[Ss]upervisor*/ScreeningSession.dm') + glob.glob(args.atlas+'/**/ScreeningSession.dm') + glob.glob(args.atlas+'/ScreeningSession.dm')

    # Then concatonate Data and Atlas Supervisor directory lists together
    supervisorList = supervisorAtlasList + supervisorDataList

    DirectoryList = pd.DataFrame({"Directory":[],"Epu":[],"Tomo":[],"Atlas":[]})
    #epuDirectoryList.append(path) # This could be a way to get the Visit
    # Add dir to directorylist if it contains EpuSession.dm files
    for f in supervisorList:

        # Supervisor folder name
        supervisor = os.path.basename(os.path.normpath(f))

        # Is it a Tomo directory
        #tomoSession = glob.glob(f+'/Batch')
        if 'mdoc' in f:
            #print(tomoSession)
            Tomo = 'Batch'
        else:
            Tomo = 'None'

        # Is it an EPU directory with metadata
        # If EPU session then grab EpuSession name from the dm file
        if 'EpuSession' in f:
            #print(epuSession)
            Epu = 'EpuSession.dm'

            with open(f, "r") as xml:
                for_parsing = xml.read()
                data = xmltodict.parse(for_parsing)
            data = data["EpuSessionXml"]

            # EPU session name
            supervisor = data["Name"]["#text"]
        else:
            Epu = 'None'

        # Is it an Atlas directory with metadata
        if '/ScreeningSession.dm' in f:
            #print(epuSession)
            Atlas = 'ScreeningSession.dm'

            with open(f, "r") as xml:
                for_parsing = xml.read()
                data = xmltodict.parse(for_parsing)
            data = data["ScreeningSessionXml"]

            # EPU session name
            supervisor = data["Name"]["#text"]
        else:
            Atlas = 'None'

        # Get full path and the basename directory
        full = os.path.realpath(f)
        dir = os.path.dirname(full)

        # Setup table describing Supervisor directories
        #DirectoryList = DirectoryList.append({'Directory':f,'Epu':Epu,'Atlas':Atlas}, ignore_index=True)
        new_row = pd.DataFrame({'Supervisor':supervisor,'Full path':full,'Directory':dir,'Epu':Epu,'Tomo':Tomo,'Atlas':Atlas}, index=[0])
        DirectoryList = pd.concat([new_row,DirectoryList.loc[:]]).reset_index(drop=True)

    # Return list of Supervisor directories that contain EpuSession.dm xmls
    return DirectoryList

def searchSupervisor(path):
    # Get list of Supervisor directories
    # Make glob insensitive to Supervisor or supervisor
    # DEV DEV DEV Concatonate path for manually transferred Supervisor directories, and murfey transferred directories into raw*/metadata
    # Sort by modified time
    supervisorList = glob.glob(path+'/[Ss]upervisor*') + glob.glob(path+'/raw*/metadata') + glob.glob(path+'/atlas/[Ss]upervisor*')
    print(supervisorList)
    supervisorList.sort(key=os.path.getmtime)

    DirectoryList = pd.DataFrame({"Directory":[],"Epu":[],"Tomo":[],"Atlas":[]})
    #epuDirectoryList.append(path) # This could be a way to get the Visit
    # Add dir to directorylist if it contains EpuSession.dm files
    for f in supervisorList:

        # Supervisor folder name
        supervisor = os.path.basename(os.path.normpath(f))

        # Is it a Tomo directory
        #tomoSession = glob.glob(f+'/Batch')
        tomoSession = glob.glob(f+'/*mdoc')
        if tomoSession:
            #print(tomoSession)
            Tomo = 'Batch'
        else:
            Tomo = 'None'

        # Is it an EPU directory with metadata
        # If EPU session then grab EpuSession name from the dm file
        epuSession = glob.glob(f+'/EpuSession.dm')
        if epuSession:
            #print(epuSession)
            Epu = 'EpuSession.dm'

            with open(f+'/EpuSession.dm', "r") as xml:
                for_parsing = xml.read()
                data = xmltodict.parse(for_parsing)
            data = data["EpuSessionXml"]

            # EPU session name
            supervisor = data["Name"]["#text"]
        else:
            Epu = 'None'

        # Is it an Atlas directory with metadata
        atlasSession = glob.glob(f+'/ScreeningSession.dm')
        if atlasSession:
            #print(epuSession)
            Atlas = 'ScreeningSession.dm'
        else:
            Atlas = 'None'

        # Setup table describing Supervisor directories
        #DirectoryList = DirectoryList.append({'Directory':f,'Epu':Epu,'Atlas':Atlas}, ignore_index=True)
        new_row = pd.DataFrame({'Supervisor':supervisor,'Directory':f,'Epu':Epu,'Tomo':Tomo,'Atlas':Atlas}, index=[0])
        DirectoryList = pd.concat([new_row,DirectoryList.loc[:]]).reset_index(drop=True)

    # Return list of Supervisor directories that contain EpuSession.dm xmls
    return DirectoryList

def searchSupervisorAtlas(path):
    print('Searching Supervisor Atlas directory for xmls, mrc, and jpg')
    print()

    xmlAtlasList = findpattern('Atlas*.xml', path) # You then need to remove any item containing *Data*
    try:
      xmlAtlas = xmlAtlasList[0]
      #print(xml)
    except:
      xmlAtlas = 'None'

    searchSupervisorAtlas.xmlAtlasList = xmlAtlasList
    searchSupervisorAtlas.xmlAtlas = xmlAtlas

    xmlAtlasTileList = findpattern('Tile*.xml', path) # You then need to remove any item containing *Data*
    try:
      xmlAtlasTile = xmlAtlasTileList[0]
      #print(xml)
    except:
      xmlAtlasTile = 'None'

    searchSupervisorAtlas.xmlAtlasTileList = xmlAtlasTileList
    searchSupervisorAtlas.xmlAtlasTile = xmlAtlasTile

    print('Found representative xml file for pulling meta data about Atlas session')
    print('Atlas: '+searchSupervisorAtlas.xmlAtlas)
    print('Atlas tile: '+searchSupervisorAtlas.xmlAtlasTile)
    print()

    # Store representative xml as global dictionary for reference anywhere in script (reduce I/O)
    with open(xmlAtlas, "r") as xml:
        for_parsing = xml.read()
        searchSupervisorAtlas.xmlAtlasDict = xmltodict.parse(for_parsing)

    with open(xmlAtlasTile, "r") as xml:
        for_parsing = xml.read()
        searchSupervisorAtlas.xmlAtlasTileDict = xmltodict.parse(for_parsing)

def searchSupervisorData(path):
    print('Searching Supervisor Data directory for xmls, mrc, and jpg')
    print()

    print('Finding GridSquare xml')
    xmlSquareList = findpattern('GridSquare*.xml', path)
    try:
      xmlSquare = xmlSquareList[0]
      print('Done')
      #print(xml)
    except:
      xmlSquare = 'None'
      print('None found')

    print('Finding FoilHole xml')
    xmlHoleList = findpattern('FoilHole*.xml', path)
    xmlHoleList = [ x for x in xmlHoleList if "Data" not in x ] # This will remove items in list containing *Data*, i.e. DataAcquisition xml files
    try:
      xmlHole = xmlHoleList[0]
      print('Done')
      #print(xml)
    except:
      xmlHole = 'None'
      print('None found')

    print('Finding AcquisitionData xml')
    xmlDataList = findpattern('FoilHole*Data*.xml', path)
    try:
      xmlData = xmlDataList[0]
      print('Done')
      #print(xml)
    except:
      xmlData = 'None'
      print('None found')

    print('Finding AquisitionData mrc')
    mrcDataList = findpattern('FoilHole*Data*.mrc', path)
    try:
      mrc = mrcDataList[0]
      print('Done')
      #print(mrc)
    except:
      mrc = 'None'
      print('None found')

    print('Finding AquisitionData jpg')
    jpgDataList = findpattern('FoilHole*Data*.jp*g', path)
    try:
      jpg = jpgDataList[0]
      print('Done')
      #print(jpg)
    except:
      jpg = 'None'
      print('None found')

    searchSupervisorData.xmlSquareList = xmlSquareList
    searchSupervisorData.xmlHoleList = xmlHoleList
    searchSupervisorData.xmlDataList = xmlDataList

    searchSupervisorData.xmlSquare = xmlSquare
    searchSupervisorData.xmlHole = xmlHole
    searchSupervisorData.xmlData = xmlData

    print('Found representative xml file for pulling meta data about EPU session')
    print('Square: '+searchSupervisorData.xmlSquare)
    print('Hole: '+searchSupervisorData.xmlHole)
    print('Acquisition: '+searchSupervisorData.xmlData)
    print()

    # Store representative xml as global dictionary for reference anywhere in script (reduce I/O)
    try:
        with open(xmlSquare, "r") as xml:
            for_parsing = xml.read()
            searchSupervisorData.xmlSquareDict = xmltodict.parse(for_parsing)
    except:
        print('searchSupervisorData error')

    try:
        with open(xmlHole, "r") as xml:
            for_parsing = xml.read()
            searchSupervisorData.xmlHoleDict = xmltodict.parse(for_parsing)
    except:
        print('searchSupervisorData error')

    try:
        with open(xmlData, "r") as xml:
            for_parsing = xml.read()
            searchSupervisorData.xmlDataDict = xmltodict.parse(for_parsing)
    except:
        print('searchSupervisorData error')

    searchSupervisorData.mrcDataList = mrcDataList

    searchSupervisorData.jpgDataList = jpgDataList

def getXmlMag(xml_path: Path) -> Dict[str, Any]:
    try:
        with open(xml_path, "r") as xml:
            for_parsing = xml.read()
            data = xmltodict.parse(for_parsing)
        data = data["MicroscopeImage"]
    except:
        xmlMag = 0
        xmlAPix = 0
    else:
        xmlMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
        xmlMetrePix = data["SpatialScale"]["pixelSize"]["x"]["numericValue"]

        xmlAPix = float(xmlMetrePix) * 1e10
        xmlAPix = roundup(xmlAPix, 1)

    return xmlMag, str(xmlAPix)

def searchSupervisorXml(path):
    # Is it an EPU directory with GridSquare directories (may not be inside Images-Disc1)
    # Here you can look for images and return the path so even without metadata you know how to find them
    xmlList = findpattern('FoilHole*Data*.xml', path)
    try:
      xml = xmlList[0]
      print(xml)
    except:
      xml = 'None'

    return xmlList

def searchSupervisorMrc(path):
    # Is it an EPU directory with GridSquare directories (may not be inside Images-Disc1)
    # Here you can look for images and return the path so even without metadata you know how to find them
    xmlList = findpattern('FoilHole*Data*.mrc', path)
    try:
      xml = mrcList[0]
      print(xml)
    except:
      xml = 'None'

    return mrcList

def searchSupervisorJpg(path):
    # Is it an EPU directory with GridSquare directories (may not be inside Images-Disc1)
    # Here you can look for images and return the path so even without metadata you know how to find them
    xmlList = findpattern('FoilHole*Data*.jp*g', path)
    try:
      xml = jpgList[0]
      print(xml)
    except:
      xml = 'None'

    return jpgList

def searchSupervisorXml(path):
    # Is it an EPU directory with GridSquare directories (may not be inside Images-Disc1)
    # Here you can look for images and return the path so even without metadata you know how to find them
    xmlList = findpattern('FoilHole*Data*.xml', path)
    try:
      xml = xmlList[0]
      print(xml)
    except:
      xml = 'None'

    return xmlList

#very very lazy stack exchanging in these find functions
#https://stackoverflow.com/questions/1724693/find-a-file-in-python

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result

# Add third input to return relative or full path, the use case above requires relative path, not full
# example - findpattern('*.txt', '/path/to/dir')
def findpattern(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    #result = glob.glob('**/'+str(pattern), recursive=True)
    return result

def findpattern_dev(pattern, path):
    files = []
    for entry in scandir.walk(path):
        for file_entry in entry.files():
            if file_entry.name.endswith('.txt'):
                txt_files.append(file_entry.path)
    return txt_files

def get_dir_size(dir, format):
    # DEV DEV DEV, this is very slow, needs to be sorted out

    # Gets whole directory size in bytes
    # Multiply by 9.53674e-7 for MB
    # Multiply by 9.31322e-10 for GB
    # Multiply by 9.09495e-13 for TB

    sizeformat = {'MB': 9.53674e-7, 'GB': 9.31322e-10, 'TB': 9.09495e-13}
    factor = sizeformat[format]

    #total = 0
    #with os.scandir(dir) as it:
    #    for entry in it:
    #        if entry.is_file():
    #            total += entry.stat().st_size
    #        elif entry.is_dir():
    #            total += get_dir_size(entry.path, format)
    #size = roundup((total * factor),3)
    #return size

    #size = 0
    #for path, dirs, files in os.walk(dir):
    #    for f in files:
    #        fp = os.path.join(path, f)
    #        size += os.path.getsize(fp)
    #for ele in os.scandir(dir):
    #    size+=os.stat(ele).st_size
    #size = roundup((size * factor),3)
    #return size

    dir_size = 0
    for f in Path(dir).glob("**/*"):
        if f.is_dir() or f.is_symlink():
            continue
        dir_size += f.stat().st_size

    size = roundup((dir_size * factor),3)
    return size

def standard_error(sessionName,realPath,date,error,raw,micNo,scope,mode):
    # Do some date conversion for reporting
    sessionDate = date.strftime("%Y-%m-%d")
    # Do some date conversion for reporting
    epuDate = date.strftime("%Y-%m-%d %H:%M:%S")

    # Get BAG name still
    lookupBAG(main.visitNameDLS)

    # These aren't working in all cases so have a bail out
    try:
       searchProcessed.ispyb
    except:
       searchProcessed.ispyb = np.nan
    try:
       searchProcessed.path
    except:
       searchProcessed.path = np.nan
    try:
       searchProcessed.number
    except:
       searchProcessed.number = np.nan

    # Session performance data
    dictHorizontal = [
    {'Proposal': lookupBAG.proposal,
     'Visit': main.visitNameDLS,
     'Visit type': main.sessionTypeDLS,
     'BAG': lookupBAG.name,
     'Institute': lookupBAG.institute,
     'Anonymous': lookupBAG.anonymous,
     'Session date': sessionDate,
     'Directory': realPath,
     'Error': error,
     'Type': mode,
     'Advisory': np.nan,
     'Instrument': scope,
     'EPU Session': sessionName,
     'EPU session date': epuDate,
     'Supervisor size (GB)': np.nan,
     'Grid_type': np.nan,
     'I0 set': np.nan,
     'I0 min int': np.nan,
     'I0 max int': np.nan,
     'Grids': np.nan,
     'Atlas': np.nan,
     'Setup time (hrs)': np.nan,
     'Total squares': np.nan,
     'Targeted squares': np.nan,
     'Collected squares': np.nan,
     'Average Foils per Square': np.nan,
     'Movie_dir': raw,
     'Movie_format': np.nan,
     'Movie_ext': np.nan,
     'Total_movies': micNo,
     'superres bin': np.nan,
     'Size (TB)': np.nan,
     'Tilt': np.nan,
     'Collection time (hrs)': np.nan,
     'Rate (mic/hr)': np.nan,
     'Total_EPU_mics': np.nan,
     'Rate (um^2/hr)': np.nan,
     'Total area (um^2)': np.nan,
     'Mag (X)': np.nan,
     'apix (A/px)': np.nan,
     'Hole (um)': np.nan,
     'Space (um)': np.nan,
     'Hole_measured (um)': np.nan,
     'Space_measured (um)': np.nan,
     'Beam (um)': np.nan,
     'Shots per hole': np.nan,
     'ShotID': np.nan,
     'Stage wait (s)': np.nan,
     'AFIS/Accurate': np.nan,
     'AFIS wait (s)': np.nan,
     'AFIS radius (um)': np.nan,
     'Exposure time (s)': np.nan,
     'Dose rate (e/A2/s)': np.nan,
     'Dose': np.nan,
     'Total dose (e/A2)': np.nan,
     'Fractions': np.nan,
     'Fraction dose (e/A2)': np.nan,
     'Defocus max (um)': np.nan,
     'Defocus min (um)': np.nan,
     'Processed': searchProcessed.ispyb,
     'Relion': searchProcessed.path,
     'Jobs': searchProcessed.number,
     'Processed size (TB)': np.nan,
     'Report path': np.nan,
     'Session report': np.nan,
     'epuVersion': np.nan}
     #'First mic': count_mics.first,
     # 'Last mic': count_mics.last,
     #'Screen time (hrs)': roundup(calcTimings.screenTimeHr,1),
    ]

    df = pd.DataFrame.from_dict(dictHorizontal)
    print('Saving data to '+main.csv_dir+'/'+sessionName+'_session.csv')
    df.to_csv (main.csv_dir+'/'+sessionName+'_session.csv', index = False, header=True)

    # Microscope optical set up
    dictHorizontal = [
    {'Proposal': lookupBAG.proposal,
     'Visit': main.visitNameDLS,
     'Visit type': main.sessionTypeDLS,
     'BAG': lookupBAG.name,
     'Institute': lookupBAG.institute,
     'Anonymous': lookupBAG.anonymous,
     'Grid_type': np.nan,
     'I0 set': np.nan,
     'I0 min int': np.nan,
     'I0 max int': np.nan,
     'Hole (um)': np.nan,
     'Space (um)': np.nan,
     'Hole_measured (um)': np.nan,
     'Space_measured (um)': np.nan,
     'Instrument': scope,
     'Session date': sessionDate,
     'Directory': realPath,
     'Error': error,
     'Type': mode,
     'Advisory': np.nan,
     'EPU Session': sessionName,
     'EPU session date': epuDate,
     'Tilt': np.nan,
     'Defocus max (um)': np.nan,
     'Defocus min (um)': np.nan,
     'Spot': np.nan,
     'C2': np.nan,
     'Beam_um': np.nan,
     'Mag_X': np.nan,
     'apix (A/px)': np.nan,
     'Slit': np.nan,
     'Slit width': np.nan,
     'Atlas_mag': np.nan,
     'Atlas_apix': np.nan,
     'Atlas_probe': np.nan,
     'Atlas_spot': np.nan,
     'Atlas_C2': np.nan,
     'Atlas_beam': np.nan,
     'Atlas_DF': np.nan,
     'Sq_mag': np.nan,
     'Sq_apix': np.nan,
     'Sq_probe': np.nan,
     'Sq_spot': np.nan,
     'Sq_C2': np.nan,
     'Sq_beam': np.nan,
     'Sq_DF': np.nan,
     'Hole_mag': np.nan,
     'Hole_apix': np.nan,
     'Hole_probe': np.nan,
     'Hole_spot': np.nan,
     'Hole_C2': np.nan,
     'Hole_beam': np.nan,
     'Hole_DF': np.nan,
     'Data_mag': np.nan,
     'Data_apix': np.nan,
     'Data_probe': np.nan,
     'Data_spot': np.nan,
     'Data_C2': np.nan,
     'Data_beam': np.nan,
     'Data_DF': np.nan,
     'epuVersion': np.nan}
    ]

    df = pd.DataFrame.from_dict(dictHorizontal)
    print('Saving data to '+main.csv_dir+'/'+sessionName+'_optics.csv')
    df.to_csv (main.csv_dir+'/'+sessionName+'_optics.csv', index = False, header=True)

def standard_output(mode):
    # It is very important that these data structures mimic the data structures in the standard_error function

    # Do some date conversion for reporting
    sessionDate = screeningSession_xml.date.strftime("%Y-%m-%d")

    # DEV DEV not great that this is done here but calculate average and stdev holes per square from data
    # Need a real function for looking through squares and calculating foilholes per squares
    #foilHoleNo = count_mics.number / xml_targets.holeSpace
    #foilHolePerSqAv = foilHoleNo / xml_targets.collectedSquares

    # Session performance data
    dictHorizontal = [
    {'Proposal': lookupBAG.proposal,
     'Visit': main.visitNameDLS,
     'Visit type': main.sessionTypeDLS,
     'BAG': lookupBAG.name,
     'Institute': lookupBAG.institute,
     'Anonymous': lookupBAG.anonymous,
     'Session date': sessionDate,
     'Directory': main.visitPath,
     'Error': 'None',
     'Type': mode,
     'Advisory': assessSession.advisoryCount,
     'Instrument': getScopeNameMag.scope,
     'EPU Session': xml_session.sessionName,
     'EPU session date': xml_session.sessionDate,
     'Supervisor size (GB)': searchSupervisor.size,
     'Grid_type': xml_session.gridType,
     'I0 set': xml_session.I0set,
     'I0 min int': xml_session.I0MinInt,
     'I0 max int': xml_session.I0MaxInt,
     'Grids': screeningSession_xml.atlasNo,
     'Atlas': 'Sample'+str(round(xml_session.autoSlot)),
     'Setup time (hrs)': roundup(calcTimings.setupTimeHr,1),
     'Total squares': xml_targets.totalSquares,
     'Targeted squares': xml_targets.targetSquares,
     'Collected squares': xml_targets.collectedSquares,
     'Average Foils per Square': xml_targets.avSqFoilHoleNo,
     'Movie_dir': findRaw.dir,
     'Movie_format': xml_session.doseFractionOutputFormat,
     'Movie_ext': count_movies.format,
     'Total_movies': count_movies.number,
     'superres bin': xml_presets_data.superResBin,
     'Size (TB)': findRaw.size,
     'Tilt': xml_presets_data.stageAlpha,
     'Collection time (hrs)': roundup(calcTimings.timeHr,1),
     'Rate (mic/hr)': calcTimings.rate,
     'Total_EPU_mics': count_mics.number,
     'Rate (um^2/hr)': calcTimings.rateArea,
     'Total area (um^2)': str(round(float(getScopeNameMag.area)*float(count_mics.number))),
     'Mag (X)': getScopeNameMag.mag,
     'apix (A/px)': getScopeNameMag.apix,
     'Hole (um)': xml_targets.holeSizeRegular,
     'Space (um)': xml_targets.holeSpaceRegular,
     'Hole_measured (um)': xml_targets.holeSize,
     'Space_measured (um)': xml_targets.holeSpace,
     'Beam (um)': xml_presets.beamD,
     'Shots per hole': xml_template_dataframe.shotsPerHole,
     'ShotID': np.nan, # DEV DEV DEV Not yet used
     'Stage wait (s)': xml_session.delayStageShift,
     'AFIS/Accurate': xml_session.afisMode,
     'AFIS wait (s)': xml_session.delayImageShift,
     'AFIS radius (um)': xml_session.afisRadius,
     'Exposure time (s)': roundup(float(xml_presets.time),3),
     'Dose rate (e/A2/s)': calc_params.doseRate,
     'Dose': calc_params.doseMessage,
     'Total dose (e/A2)': calc_params.doseTotal,
     'Fractions': calc_params.fractionNo,
     'Fraction dose (e/A2)': calc_params.doseFrameTotal,
     'Defocus max (um)': xml_session.defocusMax,
     'Defocus min (um)': xml_session.defocusMin,
     'Processed': searchProcessed.ispyb,
     'Relion': searchProcessed.path,
     'Jobs': searchProcessed.number,
     'Processed size (TB)': searchProcessed.size,
     'Report path': main.pdf_dir,
     'Session report': main.report_name,
     'epuVersion': xml_session.epuVersion}
     #'First mic': count_mics.first,
     # 'Last mic': count_mics.last,
     #'Screen time (hrs)': roundup(calcTimings.screenTimeHr,1),
    ]

    df = pd.DataFrame.from_dict(dictHorizontal)
    df.to_csv (main.csv_dir+'/'+xml_session.sessionName+'_session.csv', index = False, header=True)

    plotTimings(main.csv_dir+'/'+xml_session.sessionName+'_session.csv')

    # Microscope optical set up
    dictHorizontal = [
    {'Proposal': lookupBAG.proposal,
     'Visit': main.visitNameDLS,
     'Visit type': main.sessionTypeDLS,
     'BAG': lookupBAG.name,
     'Institute': lookupBAG.institute,
     'Anonymous': lookupBAG.anonymous,
     'Grid_type': xml_session.gridType,
     'I0 set': xml_session.I0set,
     'I0 min int': xml_session.I0MinInt,
     'I0 max int': xml_session.I0MaxInt,
     'Hole (um)': xml_targets.holeSizeRegular,
     'Space (um)': xml_targets.holeSpaceRegular,
     'Hole_measured (um)': xml_targets.holeSize,
     'Space_measured (um)': xml_targets.holeSpace,
     'Instrument': getScopeNameMag.scope,
     'Session date': sessionDate,
     'Directory': main.visitPath,
     'Error': 'None',
     'Type': mode,
     'Advisory': assessSession.advisoryCount,
     'EPU Session': xml_session.sessionName,
     'EPU session date': xml_session.sessionDate,
     'Tilt': xml_presets_data.stageAlpha,
     'Defocus max (um)': xml_session.defocusMax,
     'Defocus min (um)': xml_session.defocusMin,
     'Spot': xml_presets.spot,
     'C2': xml_presets.C2,
     'Beam_um': xml_presets.beamD,
     'Mag_X': getScopeNameMag.mag,
     'apix (A/px)': getScopeNameMag.apix,
     'Slit': xml_presets_data.filterSlit,
     'Slit width': xml_presets_data.filterSlitWidth,
     'Atlas_mag': xml_presets.magPresetList[0],
     'Atlas_apix': xml_presets.apixPresetList[0],
     'Atlas_probe': xml_presets.probePresetList[0],
     'Atlas_spot': xml_presets.spotPresetList[0],
     'Atlas_C2': xml_presets.c2PresetList[0],
     'Atlas_beam': xml_presets.beamDPresetList[0],
     'Atlas_DF': xml_presets.defocusPresetList[0],
     'Sq_mag': xml_presets.magPresetList[1],
     'Sq_apix': xml_presets.apixPresetList[1],
     'Sq_probe': xml_presets.probePresetList[1],
     'Sq_spot': xml_presets.spotPresetList[1],
     'Sq_C2': xml_presets.c2PresetList[1],
     'Sq_beam': xml_presets.beamDPresetList[1],
     'Sq_DF': xml_presets.defocusPresetList[1],
     'Hole_mag': xml_presets.magPresetList[2],
     'Hole_apix': xml_presets.apixPresetList[2],
     'Hole_probe': xml_presets.probePresetList[2],
     'Hole_spot': xml_presets.spotPresetList[2],
     'Hole_C2': xml_presets.c2PresetList[2],
     'Hole_beam': xml_presets.beamDPresetList[2],
     'Hole_DF': xml_presets.defocusPresetList[2],
     'Data_mag': xml_presets.magPresetList[3],
     'Data_apix': xml_presets.apixPresetList[3],
     'Data_probe': xml_presets.probePresetList[3],
     'Data_spot': xml_presets.spotPresetList[3],
     'Data_C2': xml_presets.c2PresetList[3],
     'Data_beam': xml_presets.beamDPresetList[3],
     'Data_DF': xml_presets.defocusPresetList[3],
     'epuVersion': xml_session.epuVersion}
    ]

    df = pd.DataFrame.from_dict(dictHorizontal)
    df.to_csv (main.csv_dir+'/'+xml_session.sessionName+'_optics.csv', index = False, header=True)

def searchProcessed(micName, sessionPath, analyse):
    ## This now needs major development
    ## Use search_mics.firstXml which note is xml and search the main directory for the raw folder  contains this xml as mrc
    ## That is then the name of the raw folder to build the processed path
    ## Search for that processed path and this is the yes/no as to whether the pipeline was used
    ## Then knowing this, pass the path to analyseProcessed function to assess relion directory structure

    # Will use the most recent (last) micrograph to find out which raw, sometimes first micrographs can be confusing, i.e. user collected in superres first and then changed parameters
    rawName = findRaw(micName, sessionPath)[0]

    if rawName == 'None':
       searchProcessed.ispyb = 'N'
       searchProcessed.path = 'None'
       searchProcessed.number = '0'
       searchProcessed.size = '0'
       print('No raw directory, skipping processed search')
    else:
       processedPath = sessionPath+'/processed/'+rawName+'/relion'

       print('Searching for processed output...')
       print('Looking in: processed/'+rawName+'/relion')
       # Need a bail out here if processed doesnt exist
       # see bi28654-10
       if os.path.exists(processedPath):
            #print(processed)
            realProcessed = os.path.realpath(processedPath)
            #print(realProcessed)
            # Count number of processing directory, relion symbolic link won't be included in this
            # You have to move up two directories, look at symbolic link structure
            number = len(os.listdir(processedPath+'/../..'))
            #print(number)
            searchProcessed.ispyb = 'Y'
            searchProcessed.path = realProcessed
            searchProcessed.number = number
            print('Processed subdirectories found, will record pipeline as having run')
            print('')

            # Global variable to delcare processed found
            searchProcessed.processedRan = 'Y'

            if analyse == 'Y':
                analyseProcessed(processedPath)
            elif analyse == 'N':
                print('Skipping full processed analysis')
                print('')
            else:
                print('Skipping full processed analysis')
                print('')

            # Get Supervisor directory size
            if args.sizecalc == 'Y':
                print('Calculating processed directory size (TB)')
                size = get_dir_size(processedPath, 'TB')
                searchProcessed.size = size
                print(str(size)+' TB')
                print()
            else:
                searchProcessed.size = np.nan
       else:
            # Global variable to delcare processed not found
            searchProcessed.processedRan = 'N'

            searchProcessed.ispyb = 'N'
            searchProcessed.path = 'None'
            searchProcessed.number = '0'
            searchProcessed.size = '0'

            analyseProcessed.motionCor = 'N'
            analyseProcessed.ctfFind = 'N'
            analyseProcessed.extract = 'N'
            analyseProcessed.class2D = 'N'
            analyseProcessed.class3D = 'N'

            print('Processed subdirectories not found, will record pipeline as not having run')
            print('')

def analyseProcessed(path):

    # Standard pipeline paths
    # Use search functions in here in case job names are not consistent across different pipeline versions
    analyseProcessed.relionPath = path
    analyseProcessed.dirStruct = 'Movies/GridSquare*/Data'
    analyseProcessed.relion_it = path+'/relion_it_options.py'
    analyseProcessed.motioncorrPath = glob.glob(path+'/MotionCorr/job0*', recursive = True)[0]
    analyseProcessed.ctffindPath = glob.glob(path+'/CtfFind/job0*', recursive = True)[0]
    #analyseProcessed.extractPath = glob.glob(path+'/Extract/job0*', recursive = True)[0] # picks may not be done so need this to behave more like 2D and 3D analysis functions
    analyseProcessed.extractPath = path+'/Extract'
    # path+'/Extract/job008' dev dev dev glob glob
    analyseProcessed.class2dPath = path+'/Class2D'
    analyseProcessed.class3dPath = path+'/Class3D'

    # Some metadata is more robust to search directly in the relion directory structure but can also pull some from relion_it.py
    if os.path.isfile(analyseProcessed.relion_it):
        analyseProcessed.params = analyseProcessed.relion_it

        # Function for pulling certain params from relion_it.py
        analyseprocessedRelionIt(analyseProcessed.params)
        print('relion_it_options.py found, pulling parameters')

        # Run these as seperate functions to (hypothesis) avoid memory issues in running altogether
        analyseProcessedMovies(analyseProcessed.motioncorrPath)
        analyseProcessedCtf(analyseProcessed.ctffindPath)

        # The next is important
        # Analyse the picks that have been extracted by relion and probably subject to filtering by cryolo thresholding, or analyse the cryolo picks directly
        # You can run both but there is an unresolved error when later reading from the main dataframe containing ideal_picking column
        #analyseProcessedPicking.type = 'relion' # Can be relion (extracted coordinates) or cryolo (all picks made by cryolo)
        analyseProcessedPicking.type = 'cryolo' # Can be relion (extracted coordinates) or cryolo (all picks made by cryolo)
        analyseProcessedPicking(analyseProcessed.extractPath, analyseProcessedPicking.type)

        analyseProcessedClassification(analyseProcessed.class2dPath)
    else:
        print('relion_it_options.py not found, uh oh')

    ## Analyse processing parameters on a per shot basis

    # Need to keep track of what per shot parameters have been analysed
    shotAnalysisList = []

    # Each analysis function above delcares a flag to say whether preprocessing results were found so base this on those flags
    # DEV DEV DEV Nedd to pass plotting color over
    if analyseProcessed.motionCor == 'Y':
        plotdir = str(main.plot_dir)+'/shotAnalysis'
        try:
            os.makedirs(plotdir)
        except:
            pass
        # Get motion per mic per shot area and plot
        analyseProcessedPerShot(xml_template_dataframe.dfShot, search_mics.sessionStructure, 'rlnAccumMotionEarly', 'tomato', plotdir, True)
        analyseProcessedPerShot(xml_template_dataframe.dfShot, search_mics.sessionStructure, 'rlnAccumMotionTotal', 'green', plotdir, True)
        shotAnalysisList.append('rlnAccumMotionEarly')
        shotAnalysisList.append('rlnAccumMotionTotal')

        # After the shot stats have been added to xml_template_dataframe.dfShot, plot these as a new template
        xml_template_dataframe_plot_stats(xml_template_dataframe.dfShot, 'rlnAccumMotionEarly', 'bwr', 'Angstrom')
        xml_template_dataframe_plot_stats(xml_template_dataframe.dfShot, 'rlnAccumMotionTotal', 'bwr', 'Angstrom')

        # Drop normalised and cmap_color columns
        xml_template_dataframe.dfShot = dropCmapDf(xml_template_dataframe.dfShot)

    if analyseProcessed.ctfFind == 'Y':
        plotdir = str(main.plot_dir)+'/shotAnalysis'
        try:
            os.makedirs(plotdir)
        except:
            pass
        # Get ctf statistics per mic per shot area and plot
        analyseProcessedPerShot(xml_template_dataframe.dfShot, search_mics.sessionStructure, 'rlnCtfMaxResolution', 'gold', plotdir, True)
        shotAnalysisList.append('rlnCtfMaxResolution')

        # After the shot stats have been added to xml_template_dataframe.dfShot, plot these as a new template
        xml_template_dataframe_plot_stats(xml_template_dataframe.dfShot, 'rlnCtfMaxResolution', 'bwr', 'Angstrom')

        # Drop normalised and cmap_color columns
        xml_template_dataframe.dfShot = dropCmapDf(xml_template_dataframe.dfShot)

    if analyseProcessed.extract == 'Y':
        plotdir = str(main.plot_dir)+'/shotAnalysis'
        try:
            os.makedirs(plotdir)
        except:
            pass
        # Get particles per mic per shot area and plot
        #analyseProcessedPerShot(xml_template_dataframe.dfShot, search_mics.sessionStructure, 'ptcls_per_mic')
        analyseProcessedPerShot(xml_template_dataframe.dfShot, search_mics.sessionStructure, 'ideal_picking', 'violet', plotdir, False)
        analyseProcessedPerShot(xml_template_dataframe.dfShot, search_mics.sessionStructure, 'pct_clustered', 'lightgreen', plotdir, False)
        shotAnalysisList.append('ideal_picking')
        shotAnalysisList.append('pct_clustered')

        # After the shot stats have been added to xml_template_dataframe.dfShot, plot these as a new template
        xml_template_dataframe_plot_stats(xml_template_dataframe.dfShot, 'ideal_picking', 'bwr_r', 'Normalised particle density')
        xml_template_dataframe_plot_stats(xml_template_dataframe.dfShot, 'pct_clustered', 'bwr', 'Percentage of clustered particles (%)')

        # Drop normalised and cmap_color columns
        xml_template_dataframe.dfShot = dropCmapDf(xml_template_dataframe.dfShot)

    if analyseProcessed.class2D == 'Y':
        plotdir = str(main.plot_dir)+'/shotAnalysis'
        try:
            os.makedirs(plotdir)
        except:
            pass
        # Get particles per mic per shot area and plot
        analyseProcessedPerShot(xml_template_dataframe.dfShot, search_mics.sessionStructure, 'Class2dRes_mean', 'deeppink', plotdir, False)
        #analyseProcessedPerShot(xml_template_dataframe.dfShot, search_mics.sessionStructure, 'Class2dRes_mode', 'deeppink', plotdir, False)
        shotAnalysisList.append('Class2dRes_mean')
        #shotAnalysisList.append('Class2dRes_mode')

        # After the shot stats have been added to xml_template_dataframe.dfShot, plot these as a new template
        xml_template_dataframe_plot_stats(xml_template_dataframe.dfShot, 'Class2dRes_mean', 'bwr', 'Average resolution of Class2D particle membership')
        #xml_template_dataframe_plot_stats(xml_template_dataframe.dfShot, 'Class2dRes_mode', 'bwr', 'Modal resolution of Class2D particle membership')

        # Drop normalised and cmap_color columns
        xml_template_dataframe.dfShot = dropCmapDf(xml_template_dataframe.dfShot)

def analyseProcessedPerShot(dfShot, dfSessionDataStructure, column, color, outdir, filter):

    print('Analysing statistics within each AFIS or Accurate exposure for: '+str(column))

    # columns for dataframe
    colColumnMean = column#+'_mean'
    colColumnCount = column#+'_count'

    # Master dataframe to put all presets into
    #dfShotStatsAll = pd.DataFrame(columns=['shotID',colColumnMean,colColumnCount])
    dfShotStatsAll = pd.DataFrame(columns=['shotID',colColumnMean])

    # Drop autofocus from the shot dataframe
    # Generally expect xml_template_dataframe.dfShot to be passed to dfShot
    dfShot = dfShot[dfShot["shotID"].str.contains("autofocus") == False]

    # Apply filtering so histograms are not dominated by processing outliers, get minimum and maximum values
    # ctfMaxReolsution often has failed fits that report back as 25 Angstroms skewing the data
    if filter == True:
        # Sigma outlier filter
        sigma = 3

        # Filter outliers beyond n sigma
        # DEV DEV DEV
        # Put this into a function for general filtering of dataframe outliers
        # It should report what the range in, and range after filtering, how many data points were removed
        dfSessionDataStructure = dfSessionDataStructure[(np.abs(stats.zscore(dfSessionDataStructure[column])) < sigma)]

        # Get min and max values from the queried column for plotting
        dfSessionDataStructure[column].replace([np.inf, -np.inf], np.nan, inplace=True)
        minimum = dfSessionDataStructure[column].dropna().min()
        maximum = dfSessionDataStructure[column].dropna().max()
        #print(minimum)
        #print(maximum)
        message = 'filtered'
    elif filter == False:
        minimum = dfSessionDataStructure[column].min()
        maximum = dfSessionDataStructure[column].max()
        message = 'unfiltered'

    print('Range: '+str(minimum)+'-'+str(maximum)+' ('+str(message)+')')
    print()

    # Box plots of statistics for all shots in the hole
    name = column.strip("_") # Strip out '_' prefix for relion parameters
    filename = str(xml_session.sessionName)+'_shotAnalysis_'+str(name)
    analyseProcessedPerShotBoxPlot(dfSessionDataStructure, column, outdir, filename)

    # Loop through shot dataframe by the shotID to plot statistics for individual shots
    for index, row in dfShot.iterrows():
        shotID = row['shotID']
        #print('Analysing '+str(shotID))

        # Filter the micrograph data structure dataframe by the shotID being worked on
        # Generally expect search_mics.sessionStructure to be passed to dfSessionDataStructure
        dfSessionDataStructureReduced = dfSessionDataStructure[dfSessionDataStructure['shotID'].str.match(str(shotID))]

        # Pull out all of the values to plot histograms of the data at each shot location
        name = colColumnMean.strip("_") # Strip out '_' prefix for relion parameters
        title = 'Micrographs with shotID '+str(shotID)+': '+colColumnMean
        yaxis = 'Number of micrographs'
        xaxis = name
        outdir = outdir # consider going in new directory for these plots
        filename = str(xml_session.sessionName)+'_shotAnalysis_'+str(name)+'_shotID_'+str(shotID)
        analyseProcessedPerShotHistogram(dfSessionDataStructureReduced, column, 40, color, color, title, yaxis, xaxis, outdir, filename, minimum, maximum)

        # Calculate mean of values in the column with column passed in 'column' for those mics
        # i.e. Average of the ptcls_per_mic for the micrographs taken with specific shotID
        # DEV DEV DEV can you do standard dev too?
        dfShotStats = dfSessionDataStructureReduced.groupby('shotID')[column].aggregate(['mean']).round(2)
        #dfShotStats = dfSessionDataStructureReduced.groupby('shotID')[column].aggregate(['mean','count']).round(2)

        dfShotStats.rename(columns = {'mean':colColumnMean}, inplace = True)
        #dfShotStats.rename(columns = {'count':colColumnCount}, inplace = True)

        # Add dataframes to master dataframe
        #dfShotStatsAll = pd.merge(dfShotStatsAll, dfShotStats, on=['shotID', colColumnMean, colColumnCount], how='outer')
        dfShotStatsAll = pd.merge(dfShotStatsAll, dfShotStats, on=['shotID', colColumnMean], how='outer')

    # Add new per shot statistics to the shot statistics dataframe
    xml_template_dataframe.dfShot = pd.merge(xml_template_dataframe.dfShot, dfShotStatsAll, on="shotID", how="left")

    # Append to saved csv
    xml_template_dataframe.dfShot.to_csv(main.csv_dir+'/'+xml_session.sessionName+'_shots.csv', index = False, header=True)

def analyseProcessedPerShotHistogram(df, column, bins, fill, edge, title, yaxis, xaxis, outdir, filename, min, max):

    #print('Plotting histogram analysis of column: '+str(column))

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    suffix = str(filename)
    yaxis = str(yaxis)
    xaxis = str(xaxis)
    title = str(title)

    # Numeric
    df_plot = df
    #df_plot[col] = pd.to_numeric(df_plot[col])

    # Plot histogram
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=column, kind='hist', kde='False',
            height=6, aspect=1.4, bins=bins, color=fill,
            edgecolor=edge, linewidth=0.5, binrange=(min, max))
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.title(title, fontsize=10)
    plt.xlim(min, max)
    plt.tight_layout()
    ax10.figure.savefig(str(outdir)+'/'+str(suffix)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

def analyseProcessedPerShotBoxPlot(df, column, outdir, filename):

    #print('Plotting boxplot analysis of column: '+str(column))

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    suffix = str(filename)

    # Plot box plot
    fig11 = plt.figure(11)
    ax11 = df.boxplot(by ='shotID', column =[column], grid = False, figsize=(5,4))
    ax11.set_ylabel(column)
    ax11.yaxis.set_label_position("left")
    ax11.yaxis.tick_left()
    ax11.tick_params(axis='x', labelrotation = 45)
    title_boxplot = str(column)
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    fig11 = ax11.get_figure()
    fig11.savefig(str(outdir)+'/'+str(suffix)+'_boxplot.png', dpi=300)
    plt.figure(11).clear()
    plt.close(11)
    #plt.show()

def analyseprocessedRelionIt(path):
    lines = open(path, "r").readlines()
    for line in lines:
        if re.search(r"do_class2d =", line, re.IGNORECASE):
            result = line.split()[2]
            if result == '1':
                request2D = 'Y'
            else:
                request2D = 'N'
        if re.search(r"do_class3d =", line, re.IGNORECASE):
            result = line.split()[2]
            if result == '1':
                request3D = 'Y'
            else:
                request3D = 'N'

    analyseProcessed.request2D = request2D
    analyseProcessed.request3D = request3D

def analyseProcessedMovies(path):

    print('Analysing preprocessing')

    # Get number of motion corrected micrographs
    if os.path.isdir(path):
        analyseProcessed.motionCor = 'Y'

        # DEV DEV known error in bi22238-40
        corStar = glob.glob(path+'/corrected_micrographs.star', recursive = True)[0]
        print('Analysing motion in '+str(corStar))

        try:
            #corMics, metadata = fileparser.getparticles(corStar)
            corMics = starToDataFrame(corStar, 'micrographs') # Alistair Burt star parser
        except:
            print('Failed to read star file')
            numFiles = '0'
            totMotionMean = '0'
            totMotionMax = '0'
            totMotionMin = '0'
            earlyMotionMean = '0'
            earlyMotionMax = '0'
            earlyMotionMin = '0'
            lateMotionMean = '0'
            lateMotionMax = '0'
            lateMotionMin = '0'
        else:
            numFiles = len(corMics)
            # Motion min and max
            totMotionMean = roundup(corMics['rlnAccumMotionTotal'].astype(float).mean(),2)
            totMotionMax = roundup(corMics['rlnAccumMotionTotal'].astype(float).max(),2)
            totMotionMin = roundup(corMics['rlnAccumMotionTotal'].astype(float).min(),2)

            earlyMotionMean = roundup(corMics['rlnAccumMotionEarly'].astype(float).mean(),2)
            earlyMotionMax = roundup(corMics['rlnAccumMotionEarly'].astype(float).max(),2)
            earlyMotionMin = roundup(corMics['rlnAccumMotionEarly'].astype(float).min(),2)

            lateMotionMean = roundup(corMics['rlnAccumMotionLate'].astype(float).mean(),2)
            lateMotionMax = roundup(corMics['rlnAccumMotionLate'].astype(float).max(),2)
            lateMotionMin = roundup(corMics['rlnAccumMotionLate'].astype(float).min(),2)

            # Plot motion for this session
            processed_plotMotion(corMics)

            ## Get motioncor values per micrograph
            # Merge to get motion values per micrograph in main session structure dataframe
            ## Add motion per mic and add to main data structure dataframe

            corMics = addDataAcquisitionColumn(corMics)

            ## DEV DEV DEV This merging is failing on certain visits - see bi28713-40
            # Now merge dataframes on basename column and search_mics.sessionStructure to get pick counts into main dataframe
            corPerMic = corMics[['DataAcquisition_name', 'rlnAccumMotionEarly', 'rlnAccumMotionLate', 'rlnAccumMotionTotal']]

            search_mics.sessionStructure = search_mics.sessionStructure.merge(corPerMic, how = 'inner', on = ['DataAcquisition_name'])

        print('MotionCorr mrc: '+str(numFiles))
        analyseProcessed.motionCorNo = numFiles
        print('Total motion range: '+str(totMotionMin)+'-'+str(totMotionMax))
        print()

        analyseProcessed.totMotionMean = totMotionMean
        analyseProcessed.totMotionMax = totMotionMax
        analyseProcessed.totMotionMin = totMotionMin

        analyseProcessed.earlyMotionMean = earlyMotionMean
        analyseProcessed.earlyMotionMax = earlyMotionMax
        analyseProcessed.earlyMotionMin = earlyMotionMin

        analyseProcessed.lateMotionMean = lateMotionMean
        analyseProcessed.lateMotionMax = lateMotionMax
        analyseProcessed.lateMotionMin = lateMotionMin

        # Plot motioncor timestamps versus acquisition timestamps
        # Should have this as an option, as it is slow
        # Also, not accurate in conjunction with fast timestamp reading of exposures
        if args.corcount == 'Y':
           countMotionCor(analyseProcessed.motioncorrPath)
           plotMotionCorTimings(analyseProcessed.motioncorrPath)
        #else:
           #countMotionCor(analyseProcessed.motioncorrPath)
    else:
        print('No motioncor found')
        analyseProcessed.motionCor = 'N'
        analyseProcessed.motionCorNo = '0'

        analyseProcessed.totMotionMean = '0'
        analyseProcessed.totMotionMax = '0'
        analyseProcessed.totMotionMin = '0'

        analyseProcessed.earlyMotionMean = '0'
        analyseProcessed.earlyMotionMax = '0'
        analyseProcessed.earlyMotionMin = '0'

        analyseProcessed.lateMotionMean = '0'
        analyseProcessed.lateMotionMax = '0'
        analyseProcessed.lateMotionMin = '0'

    # Deal with variables that might be empty
    try:
       countMotionCor.firstCollectedTime
    except:
       countMotionCor.firstCollectedTime = 'ND'
    try:
       countMotionCor.firstCorrectedTime
    except:
       countMotionCor.firstCorrectedTime = 'ND'
    try:
       countMotionCor.lastCollectedTime
    except:
       countMotionCor.lastCollectedTime = 'ND'
    try:
       countMotionCor.lastCorrectedTime
    except:
       countMotionCor.lastCorrectedTime = 'ND'

def analyseProcessedCtf(path):

    # Get number of ctffind log files and ctf params
    if os.path.isdir(path):
        analyseProcessed.ctfFind = 'Y'

        # Find ctf micrograph star file
        ctfStar = glob.glob(path+'/micrographs_ctf.star', recursive = True)[0]

        print('Analysing ctf in '+str(ctfStar))

        # Star file into dataframe
        # https://pypi.org/project/starparser/
        try:
            #ctfMics, metadata = fileparser.getparticles(ctfStar)
            ctfMics = starToDataFrame(ctfStar, 'micrographs') # Alistair Burt star parser
        except:
            print('Failed to read star file')
            numFiles = '0'
            ctfMax = '0'
            ctfMin = '0'
        else:
            # Number of star file lines
            numFiles = len(ctfMics)
            # Do some sigma value thresholding to remove outliers, look at graph histogram plotting
            # DEV DEV DEV

            # CTF max and min
            ctfMax = roundup(ctfMics['rlnCtfMaxResolution'].astype(float).max(),2)
            ctfMin = roundup(ctfMics['rlnCtfMaxResolution'].astype(float).min(),2)
            # Plot CTF for this session
            #plotPerMic(df, column, bins, fill, edge, title, yaxis, xaxis, filename)
            plotPerMic(ctfMics, 'rlnCtfMaxResolution', 'filter', 40, 'gold', 'darkorange', 'CTF max resolution per micrograph (Angstrom)', 'Micrograph count', 'CTF max resolution', '_ctf')

        ## Get ctf values per micrograph
        # Merge to get ctfmax values per micrograph in main session structure dataframe
        ## Add ptcls per mic and ideal packing score to main data structure dataframe

        # Convert MicrographName in base data acquisition name
        ctfMics = addDataAcquisitionColumn(ctfMics)

        #verifyDataStructure(search_mics.sessionStructure)
        #input('Press Enter to continue')

        # Now merge dataframes on basename column and search_mics.sessionStructure to get pick counts into main dataframe
        ctfPerMic = ctfMics[['DataAcquisition_name', 'rlnCtfMaxResolution']]
        search_mics.sessionStructure = search_mics.sessionStructure.merge(ctfPerMic, how = 'inner', on = ['DataAcquisition_name'])
        #print(search_mics.sessionStructure[['DataAcquisition_name', 'rlnCtfMaxResolution']].iloc[0])
        #input('Enter')

        print('CtfFind logs: '+str(numFiles))
        analyseProcessed.ctfFindNo = numFiles
        print('Ctf best fit resolution range: '+str(ctfMin)+'-'+str(ctfMax))
        print()
        analyseProcessed.ctfMax = ctfMax
        analyseProcessed.ctfMin = ctfMin
    else:
        print('No ctf data found')
        analyseProcessed.ctfFind = 'N'
        analyseProcessed.ctfFindNo = '0'
        analyseProcessed.ctfMax = '0'
        analyseProcessed.ctfMin = '0'

def addDataAcquisitionColumn(df):
    # Add a column which is just the data acquisition name, for use elsewhere
    # Create column with data acquisition basename, turn to path, basename to file with extension, remove _fractions suffix
    df["rlnMicrographName"] = df["rlnMicrographName"].apply(Path)
    df['DataAcquisition_name'] = df["rlnMicrographName"].apply(lambda x: os.path.basename(x))
    df['DataAcquisition_name'] = df["DataAcquisition_name"].apply(lambda x: os.path.splitext(x)[0])
    # Deal with potential suffixes on raw data
    df['DataAcquisition_name'] = df["DataAcquisition_name"].apply(lambda x: x.replace('_fractions', ''))
    df['DataAcquisition_name'] = df["DataAcquisition_name"].apply(lambda x: x.replace('_Fractions', ''))
    df['DataAcquisition_name'] = df["DataAcquisition_name"].apply(lambda x: x.replace('_EER', ''))

    return df

def analyseProcessedPicking(path, picktype):

    print('Analysing picking in :'+str(path))

    # Particle and picking info
    if os.path.isdir(path):

        # Set flag
        analyseProcessed.extract = 'Y'

        # Path is expected to be a relion path
        rlnpath = analyseProcessed.relionPath

        # Relion extract path
        extract = glob.glob(rlnpath+'/Extract/job0*', recursive = True)[0]
        # Relion particles.star path
        particleStar = glob.glob(extract+'/particles.star', recursive = True)[0]

        # Relion_it_options.py
        try:
            relion_it = rlnpath+'/relion_it_options.py'
        except:
            print('Relion_it.py required in Relion project directory')
            exit()

        # Cryolo autopick path # note two paths for old and new cryolo job structure
        #cryolo = glob.glob(rlnpath+'/AutoPick/job0*', recursive = True)[0]
        cryolo = glob.glob(rlnpath+'/AutoPick/job0*', recursive = True) + glob.glob(rlnpath+'/External/crYOLO_AutoPick/gen_pick', recursive = True)
        cryolo = cryolo[0]
        print(cryolo)
        # Cryolo cbox file path
        cboxPath = glob.glob(cryolo+'/CBOX', recursive = True)[0]
        cboxes = glob.glob(cboxPath+'/*.cbox', recursive = True)
        # Cryolo star path
        cstarPath = glob.glob(cryolo+'/STAR', recursive = True)[0]
        cstar = glob.glob(cstarPath+'/*.star', recursive = True)

        # Get particle diameter/size for clustering analysis
        #ptclDPx, analyseProcessed.allowedPx = getParticleDiameterRelion(relion_it,overlap)
        maskA, angpix = getMaskDiameterRelion(relion_it)
        analyseProcessed.angpix = angpix

        # Calculate particle parameters diameter, pixel, Angstrom and pixel size
        samplesize = 250 # This is the sample size for how many mics to read ptcl diameter from
        if len(cboxes) < samplesize:
            samplesize = len(cboxes)

        print('Number of cbox files: '+str(len(cboxes)))
        print('Sampling to get particle diameter: '+str(samplesize))

        ptclDPx = getParticleDiameterCryolo(cboxes, samplesize) # override diameter with cryolo measured diameters
        analyseProcessed.ptclDPx = roundup(ptclDPx,3)
        analyseProcessed.ptclRPx = analyseProcessed.ptclDPx/2

        ###################################
        ## Analyse particle pick clustering
        ###################################

        # Make a directory for the particles
        try:
            os.makedirs(main.plot_dir+'/particles/')
        except:
            pass

        # CLUSTERING ANALYSIS METHOD
        #method = 'dbscan'
        method = 'nn'

        # Calculate optimal packing for clustering analysis
        overlap = 0.8 # Permitted overlap as decimal
        overlapPct = overlap*100

        # Calculate from diameter what particle overlap is permitted
        analyseProcessed.ptclDAng = roundup(float(ptclDPx)*float(angpix),3)
        analyseProcessed.allowedAng = analyseProcessed.ptclDAng*overlap
        analyseProcessed.allowedPx = roundup(analyseProcessed.allowedAng/(float(angpix)),3)

        # Report
        print('Found cryolo cbox files in: ')
        print(str(cboxPath))
        print('Found cryolo star files in: ')
        print(str(cstarPath))
        print('Found relion particles files in: ')
        print(str(extract))
        print('Found relion_it_options.py files in: ')
        print(str(relion_it))
        print()
        print('Pixel size (Ang/px): '+str(angpix))
        print('Particle diameter (px): '+str(analyseProcessed.ptclDPx))
        print('Particle diameter (Ang): '+str(analyseProcessed.ptclDAng))
        print('Relion mask diameter (Ang): '+str(maskA))
        print()
        print('Permitted particle overlap ptcl diameter (%): '+str(overlapPct))
        print('Permitted particle distance (Ang): '+str(analyseProcessed.allowedAng))
        print('Permitted particle distance (px): '+str(analyseProcessed.allowedPx))
        print()

        # Read particles star file into dataframe
        # This will need to come under a conditional to state whether using relion particles.star or cryolo files
        if picktype == 'relion':
            particlesDF = particlesDataFrame(particleStar, 'relion')
            particlesDF = addDataAcquisitionColumn(particlesDF)
            analyseProcessedPicking.picktype = 'extracted'
        elif picktype == 'cryolo':
            particlesDF = particlesDataFrame(cboxPath, 'cbox')
            particlesDF = addDataAcquisitionColumn(particlesDF)
            analyseProcessedPicking.picktype = 'cryolo'
        else:
            particlesDF = particlesDataFrame(particleStar, 'relion')
            particlesDF = addDataAcquisitionColumn(particlesDF)
            analyseProcessedPicking.picktype = 'extracted'


        # If running dbscan do an optimisation to establish optimal parameter for eps, Author Jaehoon Cha SciML
        if method == 'dbscan':
            # Run eps optimisation prior to DBSCAN
            analyseProcessedPickingClustering_optimisation(particleStar, particlesDF)
        else:
            eps_opt.eps = np.nan

        # Run cluster analysis DBSCAN or NN
        particleClustering = analyseProcessedPickingClustering(particlesDF, method)

        # Reporting to terminal
        n_mics = particleClustering.shape[0]
        mean_ptcls = particleClustering[["n_particles"]].mean().round()
        mean_clustered = particleClustering[["n_particles_clustered"]].mean().round(2)
        mean_sparse = particleClustering[["n_particles_sparse"]].mean().round(2)
        mean_clustered_pct = particleClustering[["pct_clustered"]].mean().round(2)*100
        mean_sparse_pct = particleClustering[["pct_sparse"]].mean().round(2)*100
        print()
        print('Cluster analysis particle coordinates (datapoints) over micrographs (datasets)')
        print('Number of micrographs (datasets) analysed: '+str(n_mics))
        print('Average number of coordinates per micrograph: '+str(mean_ptcls['n_particles']))
        print('Average number of clustered coordinates per micrograph: '+str(mean_clustered['n_particles_clustered'])+' ('+str(mean_clustered_pct['pct_clustered'])+' %)')
        print('Average number of sparse coordinates per micrograph: '+str(mean_sparse['n_particles_sparse'])+' ('+str(mean_sparse_pct['pct_sparse'])+' %)')
        print()

        # Plot histograms of analysis
        processed_plotClustering(particleClustering)

        # Save analysis to data frame
        #df.to_csv(str(main.output)+'/particles_analysis.csv', index = False, header=True)

        # Merge particle clustering analysis with particle packing analysis

        # Convert MicrographName in base data acquisition name
        particleClustering = addDataAcquisitionColumn(particleClustering)
        # If working with a star file
        # Drop rlnMicrographName from star file dataframe as other star file merges will duplicate this column
        particleClustering.drop('rlnMicrographName', inplace=True, axis=1)

        # Now merge dataframes on basename column and search_mics.sessionStructure to get pick counts into main dataframe
        search_mics.sessionStructure = search_mics.sessionStructure.merge(particleClustering, how = 'inner', on = ['DataAcquisition_name'])

        #verifyDataStructure(search_mics.sessionStructure)
        #input('Press Enter to continue')

        # Calculate averages and stdev to put in csv report
        #print(particleClustering['pct_clustered'])
        #input()
        analyseProcessed.particlesPerMicClusterAve = roundup(particleClustering['pct_clustered'].mean(),2)
        analyseProcessed.particlesPerMicClusterStd = roundup(particleClustering['pct_clustered'].std(),2)
        print('Average percentage of particle clustering over all micrographs: '+str(analyseProcessed.particlesPerMicClusterAve)+' (+/- '+str(analyseProcessed.particlesPerMicClusterStd)+')')
        print()
        #input()

        ## DEV DEV DEV Can you use the same star file reader to look into the CRYOLO .CBOX file for average particle diameter
        ## cryolo CBOX files are star file formatted if you look!!
        ## DEV DEV DEV This is the next thing you need to do asap

        ## First step is to identify where the cbox files are

        ###################################
        ## Analyse particle packing density
        ###################################

        print('Analysing picking density in '+str(particleStar))

        # Particles star file already read into dataframe
        # Get number of particle picks
        numFiles = len(particlesDF)
        analyseProcessed.particleNo = numFiles
        print('Extracted particles: '+str(numFiles))

        ## Get picking params
        try:
            analyseProcessed.params
        except NameError:
            analyseProcessed.picking = '0'
            analyseProcessed.pickMin = '0'
            analyseProcessed.pickMax = '0'
        else:
            lines = open(analyseProcessed.params, "r").readlines()
            for line in lines:
               if re.search(r"autopick_do_cryolo", line, re.IGNORECASE):
                  if line.split()[2] == '1' or line.split()[2] == 'True':
                     analyseProcessed.picking = 'cryolo'
                  elif line.split()[2] == '0' or line.split()[2] == 'False':
                     analyseProcessed.picking = 'relion'
               if re.search(r"autopick_LoG_diam_min", line, re.IGNORECASE):
                  analyseProcessed.pickMin = line.split()[2]
               if re.search(r"autopick_LoG_diam_max", line, re.IGNORECASE):
                  analyseProcessed.pickMax = line.split()[2]

        ## Get number of particles picked per micrograph
        # Grouping to get number of particles per micrograph
        ### DEV DEV DEV error in bi29255-4
        particlesByMic = particlesDF.groupby(["rlnMicrographName"]).size().to_frame('ptcls_per_mic')
        particlesByMic = particlesByMic.reset_index()
        #print(particlesByMic['rlnMicrographName'].iloc[0])
        #input('Enter')

        # Check how many micrographs this comes from
        analyseProcessed.particleMicNo = particlesByMic.shape[0]
        print('Extracted from '+str(analyseProcessed.particleMicNo)+' micrographs')
        print()

        # Get box size
        jobStar = glob.glob(analyseProcessed.extractPath+'/**/job.star', recursive = True)[0]
        lines = open(jobStar, "r").readlines()
        for line in lines:
            if re.search(r"extract_size", line, re.IGNORECASE):
                analyseProcessed.extractPx = line.split()[1]
                analyseProcessed.extractAng = float(analyseProcessed.extractPx) * float(getScopeNameMag.apix)

        # Get size of particle by 2D mask (if not default to some value?)
        # Get mask diameter from relion_it_options.py
        #lines = open(analyseProcessed.params, "r").readlines()
        #for line in lines:
        #   if re.search(r"mask_diameter", line, re.IGNORECASE):
        #      analyseProcessed.ptclDAng = line.split()[2]

        # Calclate particle area assuming it's a square for calculating optimal packing
        particleAreaAng = float(analyseProcessed.ptclDAng) * float(analyseProcessed.ptclDAng)

        ## Calculate optimal ptcl number for idealised packing
        # Calculate the sensor size on the template image
        cameraXY = getCameraDim(search_mics.firstXml)
        cameraXang = float(cameraXY[0])*10000
        cameraYang = float(cameraXY[1])*10000
        cameraAreaang = float(cameraXang) * float(cameraYang)
        print('Camera area (Angstroms): '+str(cameraAreaang))

        print()
        print('Particle diameter (A): '+str(analyseProcessed.ptclDAng))
        print('Particle diameter (px): '+str(analyseProcessed.ptclDPx))
        print('Pixel size (A/px): '+str(analyseProcessed.angpix))
        print('Particle area as square (A^2): '+str(particleAreaAng))
        print()

        # Calculate idealised particle count for particle of n diameter
        ptclPerfectPacking = cameraAreaang / particleAreaAng
        print('Idealised pick number: '+str(ptclPerfectPacking))
        print()

        # Add a column to the particle star file dataframe for mic picks normalised to idealised ptcl density
        print('Calculating picking value normalised to ideal particle density for diameter of '+str(analyseProcessed.ptclDAng)+' Angstroms')

        # DEV DEV DEV This is a hack idea to get cryolo, this code will fail if both relion and cryolo cboxes are not both analysed
        if picktype == 'relion':
            particlesByMic['ideal_picking'] = particlesByMic['ptcls_per_mic']/ptclPerfectPacking
            particlesByMic['ideal_picking_relion'] = particlesByMic['ptcls_per_mic']/ptclPerfectPacking
        elif picktype == 'cryolo':
            particlesByMic['ideal_picking'] = particlesByMic['ptcls_per_mic']/ptclPerfectPacking
            particlesByMic['ideal_picking_cryolo'] = particlesByMic['ptcls_per_mic']/ptclPerfectPacking
        else:
            particlesByMic['ideal_picking'] = particlesByMic['ptcls_per_mic']/ptclPerfectPacking
            particlesByMic['ideal_picking_relion'] = particlesByMic['ptcls_per_mic']/ptclPerfectPacking

        # Calculate averages and stdev to put in csv report
        analyseProcessed.particlesPerMicAve = roundup(particlesByMic['ptcls_per_mic'].mean(),0)
        analyseProcessed.particlesPerMicStd = roundup(particlesByMic['ptcls_per_mic'].std(),0)
        print('Average number of particles per micrograph: '+str(analyseProcessed.particlesPerMicAve)+' (+/- '+str(analyseProcessed.particlesPerMicStd)+')')
        analyseProcessed.particlesPerMicAveNorm = roundup(particlesByMic['ideal_picking'].mean(),2)
        analyseProcessed.particlesPerMicStdNorm = roundup(particlesByMic['ideal_picking'].std(),2)
        print('Average number of particles per micrograph (idealised): '+str(analyseProcessed.particlesPerMicAveNorm)+' (+/- '+str(analyseProcessed.particlesPerMicStdNorm)+')')
        print()

        ## Add ptcls per mic and ideal packing score to main data structure dataframe

        # Convert MicrographName in base data acquisition name
        particlesByMic = addDataAcquisitionColumn(particlesByMic)
        # If working with a star file
        # Drop rlnMicrographName from star file dataframe as other star file merges will duplicate this column
        particlesByMic.drop('rlnMicrographName', inplace=True, axis=1)
        #print(particlesByMic)

        #verifyDataStructure(search_mics.sessionStructure)
        #input('Press Enter to continue')

        # Now merge dataframes on basename column and search_mics.sessionStructure to get pick counts into main dataframe
        search_mics.sessionStructure = search_mics.sessionStructure.merge(particlesByMic, how = 'inner', on = ['DataAcquisition_name'])
        #print(search_mics.sessionStructure)
        #print(search_mics.sessionStructure[['rlnMicrographName', 'ptcls_per_mic', 'ideal_picking']])
        #input()

        #verifyDataStructure(search_mics.sessionStructure)
        #input('Press Enter to continue')

        ## Plot histogram of particles per mic and idealised packing score
        processed_plotPtclsPerMic(particlesByMic)

    else:
        print('No picks found...')
        analyseProcessed.extract = 'N'
        analyseProcessed.particleNo = '0'
        analyseProcessed.extractPx = '0'
        analyseProcessed.extractAng = '0'
        analyseProcessed.picking = '0'
        analyseProcessed.pickMin = '0'
        analyseProcessed.pickMax = '0'
        analyseProcessed.particleMicNo = '0'
        analyseProcessed.particlesPerMicAve = '0'
        analyseProcessed.particlesPerMicStd = '0'
        analyseProcessed.particlesPerMicAveNorm = '0'
        analyseProcessed.particlesPerMicStdNorm = '0'
        analyseProcessed.particlesPerMicClusterAve = '0'
        analyseProcessed.particlesPerMicClusterStd = '0'
        eps_opt.eps = '0'
        analyseProcessed.ptclDAng = '0'

def getParticleDiameterRelion(relion_it, overlap):
    # Find relion_it.py file for particle diameter inferred from mask used
    # Very much subject to user error

    lines = open(relion_it, "r").readlines()
    for line in lines:
        if re.search(r"mask_diameter", line, re.IGNORECASE):
            maskA = line.split()[2]

    lines = open(relion_it, "r").readlines()
    for line in lines:
        if re.search(r"angpix", line, re.IGNORECASE):
            angpix = line.split()[2]
            break # Break because angpix appears on other lines in the relion_it_options.py file

    ptclDPx = float(maskA)/float(angpix)

    #############################################################################
    # DEV DEV DEV This parameter is key to what is considered clustered or sparse
    #############################################################################
    # Calculate permiited overlap in distance based on nearest neighbour
    allowedPx = ptclDPx*overlap

    return ptclDPx, allowedPx

def getMaskDiameterRelion(relion_it):
    # Find relion_it.py file for particle diameter inferred from mask used
    # Very much subject to user error

    lines = open(relion_it, "r").readlines()
    for line in lines:
        if re.search(r"mask_diameter", line, re.IGNORECASE):
            maskA = line.split()[2]

    lines = open(relion_it, "r").readlines()
    for line in lines:
        if re.search(r"angpix", line, re.IGNORECASE):
            angpix = line.split()[2]
            break # Break because angpix appears on other lines in the relion_it_options.py file

    return maskA, angpix

def getParticleDiameterCryolo(cboxes, sampling):

    print('Analysing cryolo cbox files to determine particle size')
    # Author Stephen Riggs, adapted by Kyle Morris
    # https://github.com/DiamondLightSource/python-relion/blob/8b0b6e0d0655790e50a0f83000ead99d52df2b8b/src/relion/zocalo/cryolo.py#L170

    # Cryolo threshold
    threshold = 0.3

    if sampling == 'all':
        cboxes_sample = cboxes
    else:
        # Get random sample of micrographs to check cryolo sizes in cbox files
        cboxes_sample = random.sample(cboxes, sampling)

    averageSizes = []

    for cbox in tqdm(cboxes_sample):

        # Find the diameters of the particles
        cryolo_particle_sizes = np.array([])
        cbox_block = cif.read_file(cbox).find_block("cryolo")

        cbox_sizes = np.append(
                np.array(cbox_block.find_loop("_EstWidth"), dtype=float),
                np.array(cbox_block.find_loop("_EstHeight"), dtype=float),
        )
        cbox_confidence = np.append(
            np.array(cbox_block.find_loop("_Confidence"), dtype=float),
            np.array(cbox_block.find_loop("_Confidence"), dtype=float),
        )
        cryolo_particle_sizes = np.append(
            cryolo_particle_sizes,
            cbox_sizes[cbox_confidence > threshold],
        )

        averageSize = np.average(cryolo_particle_sizes)
        averageSizes.append(averageSize)

        #print(averageSizes)

    # Get average size of particle in list
    #averageSize = sum(averageSizes) / len(averageSizes)
    #print('Average particle diameter: '+str(averageSize))

    # Get average size of particle in list using np to filter out any nan values arising from micrographs with no picks
    averageSize = np.nanmean(averageSizes)
    print('Average particle diameter (px): '+str(roundup(averageSize,2))+' ('+str(analyseProcessed.angpix)+' apix)')

    ptclDPx = roundup(averageSize,3)

    return ptclDPx

def getParticleCoordsCryolo(cboxes):

    # Allow all cryolo particle picks through for unbiased assessment using 0, filter out false positives using 0.3
    #threshold = 0
    threshold = 0.3

    coordinates_list = []  # List to store particle coordinates

    for cbox in tqdm(cboxes):

        # Find the coords of the particles
        cryolo_coords = np.array([])
        cbox_block = cif.read_file(cbox).find_block("cryolo")

        # In pixels
        cbox_coords_x = np.array(cbox_block.find_loop("_CoordinateX"), dtype=float)
        cbox_coords_y = np.array(cbox_block.find_loop("_CoordinateY"), dtype=float)
        cbox_box_size_x = np.array(cbox_block.find_loop("_Width"), dtype=float)
        cbox_box_size_y = np.array(cbox_block.find_loop("_Height"), dtype=float)
        cbox_ptcl_size_x = np.array(cbox_block.find_loop("_EstWidth"), dtype=float)
        cbox_ptcl_size_y = np.array(cbox_block.find_loop("_EstHeight"), dtype=float)
        cbox_confidence = np.array(cbox_block.find_loop("_Confidence"), dtype=float)

        # Filter coordinates based on confidence threshold
        filtered_coords_x = cbox_coords_x[cbox_confidence > threshold]
        filtered_coords_y = cbox_coords_y[cbox_confidence > threshold]
        filtered_box_size_x = cbox_box_size_x[cbox_confidence > threshold]
        filtered_box_size_y = cbox_box_size_y[cbox_confidence > threshold]
        filtered_ptcl_size_x = cbox_ptcl_size_x[cbox_confidence > threshold]
        filtered_ptcl_size_y = cbox_ptcl_size_y[cbox_confidence > threshold]

        # Append coordinates to the list
        for x, y, bx, by, px, py in zip(filtered_coords_x, filtered_coords_y, filtered_box_size_x, filtered_box_size_y, filtered_ptcl_size_x, filtered_ptcl_size_y):
            coordinates_list.append({'rlnCoordinateX': x, 'rlnCoordinateY': y, 'rlnMicrographName': os.path.basename(cbox), 'boxX': bx, 'boxY': by, 'ptclX': px, 'ptclY': py})

    # Create a DataFrame from the list of dictionaries
    df = pd.DataFrame(coordinates_list)

    # Report in cryolo picks extracted and made into star like dataframe
    num_rows = df.shape[0]
    print('Particles in cryolo cboxes above threshhold ('+str(threshold)+'): '+str(num_rows))

    return df

def starToDataFrame(path, datablock):
    #print('Reading star file: '+str(os.path.basename(Path(path))))
    starDF = starfile.read(path) # Star file is returned as ordered dictionary continaining dataframes
    datablockDF = starDF[datablock]
    return datablockDF

def particlesDataFrame(path, method):
    # This retrieves the particle co-ordinates into a dataframe for cluster analysis
    if method == 'relion':
        #analyseProcessed.extract = 'Y'
        print('Analysing picking in relion file: '+str(path))

        # Get total number of extracted particles
        # https://pypi.org/project/starparser/
        try:
            # This will be a dataframe of the particles.star file
            # Includes all fields you'd expect, mic name, x,y coords
            #particlesDF, metadata = fileparser.getparticles(path)
            particlesDF = starToDataFrame(path, 'particles') # Alistair Burt star parser
        except:
            print('Fail')
            numFiles = np.nan
        else:
            numFiles = len(particlesDF)
        #analyseProcessed.particleNo = numFiles

        print('Extracted particles: '+str(numFiles))
        print()

        ## Get number of particles picked per micrograph
        # Grouping to get number of particles per micrograph
        ### DEV DEV DEV error in bi29255-4
        particlesByMic = particlesDF.groupby(["rlnMicrographName"]).size().to_frame('ptcls_per_mic')
        particlesByMic = particlesByMic.reset_index()
        #print(particlesByMic['rlnMicrographName'].iloc[0])
        #input('Enter')

        #print(particles["rlnMicrographName"])

        return particlesDF

    elif method == 'cbox':
        # DEV DEV DEV not finished this
        print('Analysing picking in cryolo cbox files in: '+str(path))
        cboxPath = path
        cboxes = glob.glob(cboxPath+'/*.cbox', recursive = True)

        particlesDF = getParticleCoordsCryolo(cboxes)

        # Get box cryolo box size for adjustment of coordinates
        cryoloBoxSize = particlesDF['boxX'].iloc[0]
        cryoloBoxHalf = cryoloBoxSize / 2

        # Adjust the coordinates to be centre of image rather than lower left quandrant
        particlesDF['rlnCoordinateX'] = particlesDF['rlnCoordinateX']+cryoloBoxHalf
        particlesDF['rlnCoordinateY'] = particlesDF['rlnCoordinateY']+cryoloBoxHalf

        return particlesDF

def analyseProcessedPickingClustering(particles, method):

    # Copy random mics into plot folder for comparison to the clustering analysis
    shutil.copy(str(search_mics.randomName0Path)+'.jpg', main.plot_dir+'/particles/')
    shutil.copy(str(search_mics.randomName1Path)+'.jpg', main.plot_dir+'/particles/')
    shutil.copy(str(search_mics.randomName2Path)+'.jpg', main.plot_dir+'/particles/')
    shutil.copy(str(search_mics.randomName3Path)+'.jpg', main.plot_dir+'/particles/')

    ## Pull out the x,y coordinates and calculate sparsity or aggregation
    # Get unique micrograph names
    micrographNames = particles.groupby(["rlnMicrographName"]).agg(['unique'])
    micrographNames = micrographNames.reset_index()

    # Empty df for gathering
    analyseProcessedPickingClustering.micname_list = []
    analyseProcessedPickingClustering.n_particles_list = []
    analyseProcessedPickingClustering.n_particles_clustered_list = []
    analyseProcessedPickingClustering.n_particles_sparse_list = []
    analyseProcessedPickingClustering.pct_clustered_list = []
    analyseProcessedPickingClustering.pct_sparse_list = []
    analyseProcessedPickingClustering.n_clusters__list = []

    # This will need to come under a conditional input argument as to whether to use nearest neighbours or dbscan or something else
    #micrographNames.progress_apply(lambda row : print(row), axis = 1)
    if method == 'dbscan':
        print('Analysing co-ordinate clustering using dbscan')
        micrographNames.progress_apply(lambda row : process_dbscan(row), axis = 1)
    elif method == 'nn':
        print('Analysing co-ordinate clustering using nearest neighbour analysis')
        micrographNames.progress_apply(lambda row : process_nn(row), axis = 1)

    df = pd.DataFrame({'rlnMicrographName': analyseProcessedPickingClustering.micname_list, 'n_particles': analyseProcessedPickingClustering.n_particles_list, 'n_particles_clustered': analyseProcessedPickingClustering.n_particles_clustered_list, 'n_particles_sparse': analyseProcessedPickingClustering.n_particles_sparse_list, 'pct_clustered': analyseProcessedPickingClustering.pct_clustered_list, 'pct_sparse': analyseProcessedPickingClustering.pct_sparse_list, 'n_clusters_': analyseProcessedPickingClustering.n_clusters__list})
    del analyseProcessedPickingClustering.micname_list, analyseProcessedPickingClustering.n_particles_list, analyseProcessedPickingClustering.n_particles_clustered_list, analyseProcessedPickingClustering.n_particles_sparse_list, analyseProcessedPickingClustering.pct_clustered_list, analyseProcessedPickingClustering.pct_sparse_list, analyseProcessedPickingClustering.n_clusters__list

    # DEV DEV DEV Now that dataframe needs to be added to the dataframe containing info on every micrograph
    #input()
    return(df)

def analyseProcessedPickingClustering_optimisation(particleStar, particles):
    ### get star files paths
    starfiles = particleStar
    print(starfiles)
    time_start = time()
    particles_dic = {}

    dataset_path = ('/').join(starfiles.split('/')[:2])
    particles_dic[dataset_path] = {}

    ### get df from star files
    #particles, metadata = fileparser.getparticles(starfiles)
    particles = starToDataFrame(starfiles, 'particles') # Alistair Burt star parser
    particles_dic[dataset_path]['particles'] = particles
    particles_dic[dataset_path]['particlesByMic']  = particles_dic[dataset_path]['particles'].groupby(["rlnMicrographName"]).size().to_frame('ptcls_per_mic')
    particles_dic[dataset_path]['particlesByMic']  = particles_dic[dataset_path]['particlesByMic'] .reset_index()

    ### get XY coordinate in the dataset
    particles_dic[dataset_path]['XY_Micrograph'] = particles_dic[dataset_path]['particles'][['rlnCoordinateX', 'rlnCoordinateY']]

    ### get the mean of samples per micrograph
    n_rnd_sample = int(particles_dic[dataset_path]['particles'].groupby(["rlnMicrographName"]).count()['rlnCoordinateX'].mean())
    rnd_idx = np.random.choice(len(particles_dic[dataset_path]['XY_Micrograph']), n_rnd_sample, replace = False)

    ### make eps vs # of clusters image save folder per dataset
    folder_path = os.path.join(main.plot_dir, 'eps')
    ptcl_plot_path = os.path.join(main.plot_dir, 'particles')
    try:
        os.mkdir(folder_path)
    except OSError:
        pass

    try:
        os.mkdir(ptcl_plot_path)
    except OSError:
        pass

    ### obtain the peak eps and images
    df = particles_dic[dataset_path]['XY_Micrograph'].iloc[rnd_idx]
    particles_dic[dataset_path]['eps'] = eps_opt(df, folder_path)

    time_end = time()
    print('Optimising eps value took {}'.format(time_end - time_start))
    print('Optimised eps value: '+str(eps_opt.eps))

def eps_opt(df, folder_path):
    xy = df.values.astype(np.float32)
    xy = StandardScaler().fit_transform(xy)

    ### get the number of clusters w.r.t. eps values from 0.01 to 1.
    n_clusters_dic = {}
    for eps in np.linspace(0.01, 1., num = 100):
        db = DBSCAN(eps=eps).fit(xy)
        labels = db.labels_
        n_clusters_dic[eps] = len(np.unique(labels)) - 1

    eps = np.array(list(n_clusters_dic.keys()))
    n_clusters = np.array(list(n_clusters_dic.values()))
    peak_eps = eps[np.argmax(n_clusters)]
    plt.figure(figsize = (10, 4))
    plt.plot(eps, n_clusters)
    plt.axvline(x = peak_eps, c = 'red')
    plt.text(peak_eps, 2.0, 'EPS at the peak is {}'.format(peak_eps))
    plt.title('EPS vs # of clusters')
    plt.xlabel('EPS')
    plt.ylabel('# of clusters')
    filename = str(xml_session.sessionName)+str('_rnd_sample_eps_opt.png')
    img_path = os.path.join(folder_path, str(filename))
    plt.savefig(img_path)
    plt.close()
    eps_opt.eps = peak_eps
    main.epsilon = peak_eps
    return peak_eps

def process_dbscan(row):
    #print(row)
    #print(type(row))
    #input()
    # row is returned as a pandas.core.series.Series, some values are simpler to extract values from a dataframe so convert here
    #df = row.to_frame().T
    #print(df)
    #print(type(df))
    #input()
    # Pull mic name, row is dtype series object so [0] pulls out as string for use later
    #micname = row[["rlnMicrographName"]][0]
    micname = row['rlnMicrographName'].iloc[0]
    #print(micname)
    #print(type(micname))
    #input()
    #micname_tolist = micname.tolist()
    #print(micname_tolist)
    #print(type(micname_tolist))
    #input()

    dataname = row['DataAcquisition_name']#.iloc[0]
    #print(dataname)
    #print(type(dataname))
    #input()

    # These return as list
    micCoordsX = row["rlnCoordinateX"][0]
    micCoordsY = row["rlnCoordinateY"][0]
    #print(len(micCoordsX))
    #print(len(micCoordsY))
    #print(micCoordsX)
    #print(micCoordsY)

    # Sometimes there are different numbers of x and y coordinates???
    if len(micCoordsX) == len(micCoordsY):
        # Put back into sensible dataframe
        micCoords = pd.DataFrame({'rlnCoordinateX': micCoordsX, 'rlnCoordinateY': micCoordsY})

        #print(micCoords)
        #input()

        n_particles, n_particles_clustered, n_particles_sparse, pct_clustered, pct_sparse, n_clusters_ = dbscan(micCoords, micname)

    else:

        n_particles = np.nan
        n_particles_clustered = np.nan
        n_particles_sparse = np.nan
        pct_clustered = np.nan
        pct_sparse = np.nan
        n_clusters_ = np.nan

    analyseProcessedPickingClustering.micname_list.append(micname)
    analyseProcessedPickingClustering.n_particles_list.append(n_particles)
    analyseProcessedPickingClustering.n_particles_clustered_list.append(n_particles_clustered)
    analyseProcessedPickingClustering.n_particles_sparse_list.append(n_particles_sparse)
    analyseProcessedPickingClustering.pct_clustered_list.append(pct_clustered)
    analyseProcessedPickingClustering.pct_sparse_list.append(pct_sparse)
    analyseProcessedPickingClustering.n_clusters__list.append(n_clusters_)

def dbscan(df, micname):

    # Readthedocs
    # https://scikit-learn.org/stable/modules/clustering.html
    # https://scikit-learn.ru/example/demo-of-dbscan-clustering-algorithm/

    # This pulls the dataframe into an np array
    X = df[['rlnCoordinateX', 'rlnCoordinateY']].values
    X = StandardScaler().fit_transform(X)

    # #############################################################################
    # Compute DBSCAN
    db = DBSCAN(eps=eps_opt.eps, min_samples=10).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    n_particles = len(X)
    n_particles_clustered = n_particles-n_noise_
    n_particles_sparse = n_noise_
    pct_clustered = n_particles_clustered/n_particles
    pct_sparse = n_noise_/n_particles

    # Report on analysis
    #print("Number of particles coordinates: "+str(n_particles))
    #print('Estimated number of clustered particles: %d' % n_particles_clustered+ " (" + str(round(pct_clustered*100)) + " %)")
    #print('Estimated number of sparse particles: %d' % n_particles_sparse+ " (" + str(round(pct_sparse*100)) + " %)")
    #print()
    #print('Estimated number of clusters: %d' % n_clusters_)
    #print('Estimated number of noise points: %d' % n_noise_)
    #print()

    # DEV DEV DEV here the micname needs to be compared to the random mic names and if they match to set plot to yes
    # This way you plot particle aggregation for the micrographs that are also shown in the report
    plot = 'n'

    # Do a string comparison to the randomly selected micrographs made earlier, if analysing one of the random mics, plot the clustering analysis
    name = Path(micname).stem
    currentmic = name.replace('_fractions', '')
    #currentmic = name.replace('_EER','')

    if currentmic == search_mics.random0Name or currentmic == search_mics.random1Name or currentmic == search_mics.random2Name or currentmic == search_mics.random3Name:
        plot = 'y'

    if plot == 'y':
        #print('Plotting particle clustering')
        dbscanplot(X, core_samples_mask, labels, n_clusters_, micname)

    # Return data from function
    return n_particles, n_particles_clustered, n_particles_sparse, pct_clustered, pct_sparse, n_clusters_

def dbscanplot(X, core_samples_mask, labels, n_clusters_, micname):

    # #############################################################################
    # Plot result
    import matplotlib.pyplot as plt

    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
            for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_member_mask = (labels == k)

        xy = X[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                markeredgecolor='k', markersize=8)

        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                markeredgecolor='k', markersize=8)

    plt.title('Estimated number of clusters by DBSCAN: ' + str(n_clusters_) + ' (eps: ' +str(eps_opt.eps)+ ')')

    filename = os.path.basename(Path(micname))

    plt.savefig(main.plot_dir+'/particles/'+str(filename)+'_eps_'+str(eps_opt.eps)+'_particles.png', dpi=300)
    plt.clf()
    #plt.show()

def process_nn(row):
    micname = row['rlnMicrographName'].iloc[0]

    dataname = row['DataAcquisition_name']#.iloc[0]

    # These return as list
    micCoordsX = row["rlnCoordinateX"][0]
    micCoordsY = row["rlnCoordinateY"][0]

    # Sometimes there are different numbers of x and y coordinates???
    if len(micCoordsX) == len(micCoordsY):
        # Put back into sensible dataframe
        micCoords = pd.DataFrame({'rlnCoordinateX': micCoordsX, 'rlnCoordinateY': micCoordsY})

        radius = float(analyseProcessed.allowedPx)
        n_particles, n_particles_clustered, n_particles_sparse, pct_clustered, pct_sparse, n_clusters_ = nn_labels(micCoords, micname, radius)

    else:

        n_particles = np.nan
        n_particles_clustered = np.nan
        n_particles_sparse = np.nan
        pct_clustered = np.nan
        pct_sparse = np.nan
        n_clusters_ = np.nan

    analyseProcessedPickingClustering.micname_list.append(micname)
    analyseProcessedPickingClustering.n_particles_list.append(n_particles)
    analyseProcessedPickingClustering.n_particles_clustered_list.append(n_particles_clustered)
    analyseProcessedPickingClustering.n_particles_sparse_list.append(n_particles_sparse)
    analyseProcessedPickingClustering.pct_clustered_list.append(pct_clustered)
    analyseProcessedPickingClustering.pct_sparse_list.append(pct_sparse)
    analyseProcessedPickingClustering.n_clusters__list.append(n_clusters_)

def nn_labels(df, micname, radius):
    """
    Calculates labels for the coordinates based on the minimum nearest neighbor distance.

    Args:
        df (pandas.DataFrame): A DataFrame containing 'x' and 'y' columns representing the coordinates.
        radius (float): The radius to determine if a coordinate is 'close' or 'sparse'.

    Returns:
        pandas.DataFrame: A new DataFrame with added 'label' column.
    """
    df_numeric = df[['rlnCoordinateX', 'rlnCoordinateY']].copy()

    if len(df_numeric) > 1:
        nbrs = NearestNeighbors(n_neighbors=min(2, len(df_numeric))).fit(df_numeric)
        distances, indices = nbrs.kneighbors(df_numeric)

        df_new = df.copy()
        df_new['_min_distance'] = distances[:, 1]  # Distance to the second nearest neighbor

        df_new['label'] = df_new.apply(lambda row: 'close' if row['_min_distance'] <= radius else 'sparse', axis=1)
    else:
        df_new = df.copy()
        df_new['label'] = 'sparse'  # Label all points as 'sparse' if only 1 data point

    #clusters = estimate_clusters(df)

    n_particles = df_new.shape[0]
    n_particles_clustered = len(df_new[df_new['label'].str.contains('close')])
    n_particles_sparse = len(df_new[df_new['label'].str.contains('sparse')])
    pct_clustered = n_particles_clustered/n_particles
    pct_sparse = n_particles_sparse/n_particles
    n_clusters_ = np.nan

    # DEV DEV DEV here the micname needs to be compared to the random mic names and if they match to set plot to yes
    # This way you plot particle aggregation for the micrographs that are also shown in the report
    plot = 'n'
    micpath = None

    # Do a string comparison to the randomly selected micrographs made earlier, if analysing one of the random mics, plot the clustering analysis
    name = Path(micname).stem
    currentmic = name.replace('_fractions', '')
    currentmic = currentmic.replace('_EER','')

    if currentmic == search_mics.random0Name:
        plot = 'y'
        micpath = str(search_mics.randomName0Path)+'.jpg'

    if currentmic == search_mics.random1Name:
        plot = 'y'
        micpath = str(search_mics.randomName1Path)+'.jpg'

    if currentmic == search_mics.random2Name:
        plot = 'y'
        micpath = str(search_mics.randomName2Path)+'.jpg'

    if currentmic == search_mics.random3Name:
        plot = 'y'
        micpath = str(search_mics.randomName3Path)+'.jpg'

    # Plot for individual micrograph
    if plot == 'y':
        #plot_colored_nn_coordinates_old(df, micname)
        plot_colored_nn_coordinates(df_new, analyseProcessed.ptclDPx, 0.3, 0.7, 'black', micpath, micname)
        #plot_colored_nn_coordinates_new(df, marker_size=analyseProcessed.ptclDPx, center_alpha=0.1, edge_alpha=0.9, edge_color='black', background_image_path=None, micname=micname)

    #return df
    return n_particles, n_particles_clustered, n_particles_sparse, pct_clustered, pct_sparse, n_clusters_

def plot_colored_nn_coordinates(df, marker_size, center_alpha, edge_alpha, edge_color, background_image_path, micname):

    # Values are pulled as strings and need to be numeric
    df_plot = df
    df_plot["rlnCoordinateX"] = pd.to_numeric(df_plot["rlnCoordinateX"])
    df_plot["rlnCoordinateY"] = pd.to_numeric(df_plot["rlnCoordinateY"])

    # Make a plot of the hole, using the measured hole diameter to input radius of circle
    # https://www.codespeedy.com/how-to-draw-shapes-in-matplotlib-with-python/
    fig2 = plt.figure(2)
    ax1 = fig2.add_subplot()

    # Create legend elements
    red_circle = Circle((0.5, 0.5), radius=marker_size/2, facecolor='red', edgecolor=edge_color, alpha=center_alpha, linewidth=0.5, label='clustered')
    black_circle = Circle((0.5, 0.5), radius=marker_size/2, facecolor='black', edgecolor=edge_color, alpha=center_alpha, linewidth=0.5, label='sparse')

    # Calculate the aspect ratio to preserve the circles' shape
    aspect_ratio = float(getCameraDim.cameraXpx) / float(getCameraDim.cameraYpx)

    # Adjust the figure size to extend the Y-axis for the legend
    plotSizeX = 5
    plotSizeY = plotSizeX / aspect_ratio

    plt.gcf().set_size_inches(plotSizeX, plotSizeY)

    # Set the aspect ratio of the plot to be equal
    ax1.set_aspect('equal')

    # Assign colors to the clustering labels
    colors = {'close': 'red', 'sparse': 'black'}

    # Iterate through coordinates in dataframe and plot as circles on plot, color by nearest neighbour label close or sparse
    for index, row in df_plot.iterrows():
        # Load and set background image if provided
        #if background_image_path:
            #print('Merge image for this one')
            #img = plt.imread(background_image_path)
            #plt.imshow(img, extent=[df_plot['rlnCoordinateX'].min(), df_plot['rlnCoordinateX'].max(), df_plot['rlnCoordinateY'].min(), df_plot['rlnCoordinateY'].max()], aspect='auto', alpha=0.5)
            #ax1 = plt.gca()


        #circle = plt.Circle((row['rlnCoordinateX'], row['rlnCoordinateY']), radius=analyseProcessed.ptclDPx/2, color="black", alpha = 0.3)
        circleFace = plt.Circle(
            (row['rlnCoordinateX'], row['rlnCoordinateY']),
            radius=marker_size/2,
            facecolor=colors[row['label']],
            alpha=center_alpha)
        ax1 = plt.gca()
        ax1.add_patch(circleFace)

        circleEdge = plt.Circle(
            (row['rlnCoordinateX'], row['rlnCoordinateY']),
            radius=marker_size/2,
            facecolor='none',
            edgecolor=edge_color,
            alpha=edge_alpha,
            linewidth=0.5)
        ax1 = plt.gca()
        ax1.add_patch(circleEdge)

    # Add a rectangle patch to create a border
    rect = Rectangle((0, 0), float(getCameraDim.cameraXpx), float(getCameraDim.cameraYpx), linewidth=1, edgecolor='black', facecolor='none', alpha=0.5)
    #rect = Rectangle((df_plot['rlnCoordinateX'].min(), df_plot['rlnCoordinateY'].min()), df_plot['rlnCoordinateX'].max() - df_plot['rlnCoordinateX'].min(), df_plot['rlnCoordinateY'].max() - df_plot['rlnCoordinateY'].min(), linewidth=1, edgecolor='black', facecolor='none', alpha=0.5)
    plt.gca().add_patch(rect)

    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)

    plt.xlim(0,float(getCameraDim.cameraXpx))
    plt.ylim(0,float(getCameraDim.cameraYpx))

    # Add the legend to the plot
    plt.legend(handles=[red_circle, black_circle], loc='lower center', bbox_to_anchor=(0.5, -0.15), fancybox=False, shadow=False, ncol=2)

    #ax1.axis('scaled')
    plt.tight_layout()

    #Save or display the plot
    filename = os.path.basename(Path(micname))
    tmp1 = filename.replace('_fractions.mrc','')
    tmp2 = tmp1.replace('_fractions.cbox','')
    tmp3 = tmp2.replace('_EER.mrc','')
    name = tmp3.replace('_EER.cbox','')

    name = name+'_'+str(analyseProcessedPicking.picktype)
    #plt.title('Nearest neighbour clustering analysis')
    #plt.show()
    fig2.savefig(main.plot_dir+'/particles/'+str(name)+'_nnclustering_particles_trans.png', dpi=300, transparent=True, bbox_inches='tight')
    fig2.savefig(main.plot_dir+'/particles/'+str(name)+'_nnclustering_particles.png', dpi=300, transparent=False, bbox_inches='tight')
    plt.clf()
    plt.figure(2).clear()
    plt.close('all')

def processed_plotClustering(df):

    picktype = analyseProcessedPicking.picktype

    print('Plotting particle cluster analysis: percentage of clustered particles')
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = 'pct_clustered'
    #filename = str(xml_session.sessionName)+'_particle_cluster'
    suffix = '_ptcls_per_mic_clustered_'+str(picktype)
    filename = str(xml_session.sessionName)+str(suffix)
    xaxis = 'Percentage of particles per micrograph in cluster'
    title = 'Particle clustering, closer than '+str(roundup(analyseProcessed.allowedAng,1))+' A'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    #print(filtered)
    # Quartile outliers
    #q_low = df_plot[col].quantile(0.01)
    #q_hi  = df_plot[col].quantile(0.99)
    #df_plot = df_plot[(df_plot[col] < q_hi) & (df_plot[col] > q_low)]
    # Need to implement subplotting to show totalmotion, initialmotion and latemotion
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='lightgreen',
            binrange=(0,1),
            edgecolor="navy", linewidth=0.5).set(title=title)

    ax10.set(xlabel=xaxis, ylabel='Micrograph count')
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    print('Plotting particle cluster analysis: percentage of unclustered particles')
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = 'pct_sparse'
    #filename = str(xml_session.sessionName)+'_particle_cluster'
    suffix = '_ptcls_per_mic_unclustered_'+str(picktype)
    filename = str(xml_session.sessionName)+str(suffix)
    xaxis = 'Percentage of particles per micrograph not in cluster'
    title = 'Particle sparsity, further than '+str(analyseProcessed.allowedAng)+' A'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    #print(filtered)
    # Quartile outliers
    #q_low = df_plot[col].quantile(0.01)
    #q_hi  = df_plot[col].quantile(0.99)
    #df_plot = df_plot[(df_plot[col] < q_hi) & (df_plot[col] > q_low)]
    # Need to implement subplotting to show totalmotion, initialmotion and latemotion
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='lightgreen',
            binrange=(0,1),
            edgecolor="navy", linewidth=0.5).set(title=title)

    ax10.set(xlabel=xaxis, ylabel='Micrograph count')
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    print('Plotting particle cluster analysis: number of clusters')
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = 'n_clusters_'
    #filename = str(xml_session.sessionName)+'_particle_cluster'
    suffix = '_ptcls_per_mic_cluster_number_'+str(picktype)
    filename = str(xml_session.sessionName)+str(suffix)
    xaxis = 'Particle clustering: number of clusters'
    title = 'Particle clustering per micrograph: number of clusters'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    #print(filtered)
    # Quartile outliers
    #q_low = df_plot[col].quantile(0.01)
    #q_hi  = df_plot[col].quantile(0.99)
    #df_plot = df_plot[(df_plot[col] < q_hi) & (df_plot[col] > q_low)]
    # Need to implement subplotting to show totalmotion, initialmotion and latemotion
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='lightgreen',
            edgecolor="navy", linewidth=0.5).set(title=title)

    ax10.set(xlabel=xaxis, ylabel='Micrograph count')
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Multi motion plot
    suffix = '_ptcls_per_mic_cluster_analysis_'+str(picktype)
    filename = str(xml_session.sessionName)+str(suffix)

    col0 = 'pct_clustered'
    col1 = 'pct_sparse'
    #col2 = 'n_clusterd'
    xaxis0 = 'Clustered (%)'
    xaxis1 = 'Unclustered (%)'
    #xaxis2 = 'Clusters (n)'

    df_plot[col0] = pd.to_numeric(df_plot[col0])
    df_plot[col1] = pd.to_numeric(df_plot[col1])
    #df_plot[col2] = pd.to_numeric(df_plot[col2])

    filtered0 = df_plot[(np.abs(stats.zscore(df_plot[col0])) < 3.0)]
    filtered1 = df_plot[(np.abs(stats.zscore(df_plot[col1])) < 3.0)]
    #filtered2 = df_plot[(np.abs(stats.zscore(df_plot[col2])) < 3.0)]

    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(6, 8))
    figure.suptitle('Particle clustering analysis')
    axes[0].set_title(xaxis0)
    axes[1].set_title(xaxis1)
    #axes[2].set_title(xaxis2)

    axes[0].set_ylabel('Micrographs')
    axes[1].set_ylabel('Micrographs')
    #axes[2].set_ylabel('Micrographs')

    axes[1].set_xlabel('Percentage of particle picks per micrograph')

    sns.histplot(data=df_plot, x=col0, kde='False',
            bins=40, color='tomato',
            binrange=(0,1),
            edgecolor="navy", linewidth=0.5, ax = axes[0])
    sns.histplot(data=df_plot, x=col1, kde='False',
            bins=40, color='royalblue',
            binrange=(0,1),
            edgecolor="navy", linewidth=0.5, ax = axes[1])
    #sns.histplot(data=filtered2, x=col2, kde='False',
    #        bins=40, color='lightgreen',
    #        binrange=(0,1),
    #        edgecolor="navy", linewidth=0.5, ax = axes[2])

    figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)

def analyseProcessedClassification(path):

    # Get Class2D directories
    class2djobs = sorted(glob.glob(path+'/job**', recursive = True))

    if not class2djobs:
        print('Class2D jobs not found')
        analyseProcessed.class2D = 'N'
        analyseProcessed.ptclDAng = '0'
        analyseProcessed.maxProbMax = '0'
        analyseProcessed.maxProbMin = '0'
        analyseProcessed.class2dresMax = '0'
        analyseProcessed.class2dresMin = '0'
        analyseProcessed.class2dresMean = '0'
        analyseProcessed.class2dresMode = '0'
    else:
        analyseProcessed.class2D = 'Y'

        print('Analysing 2D classification in: '+str(path))
        print()

        ## Assign Class2D resolutions and occupancies to their respectve particles
        # This dataframe could be used elsewhere to plot average Class2D res in shot areas
        class2DPtclsDf = assignClass2dParam2Ptcl(class2djobs, 'rlnEstimatedResolution', 'rlnEstimatedResolutionClass2D')
        class2DPtclsDf2 = assignClass2dParam2Ptcl(class2djobs, 'rlnClassDistribution', 'rlnClassDistributionClass2D')
        # Merge
        class2DPtclsDf['rlnClassDistributionClass2D'] = class2DPtclsDf2['rlnClassDistributionClass2D']

        # Plot the 2D class resolution each particle is assigned to
        print('Coloring and plotting coordinates of particles with Class2D class property assignments')
        print()

        # Get only the randomly selected mics so a color mapping can be applied that is relvant to them all in comparison
        search_strings = [search_mics.random0Name, search_mics.random1Name, search_mics.random2Name, search_mics.random3Name]
        df = class2DPtclsDf[class2DPtclsDf['rlnMicrographName'].str.contains('|'.join(search_strings))]
        # Replace inf values that interfere with color mapping
        df['rlnEstimatedResolutionClass2D'].replace([np.inf, -np.inf], np.nan, inplace=True)
        df['rlnClassDistributionClass2D'].replace([np.inf, -np.inf], np.nan, inplace=True)

        # Make plots of particles with Class2D property assignments
        plotClass2Dptcls(df, 'rlnEstimatedResolutionClass2D')
        plotClass2Dptcls(df, 'rlnClassDistributionClass2D')
        plotClass2Dptcls(df, 'rlnClassNumber')

        ## Assess Class2D
        # Find rlnMaxValueProbDist
        maxProbMax = roundup(class2DPtclsDf['rlnMaxValueProbDistribution'].astype(float).max(),2)
        maxProbMin = roundup(class2DPtclsDf['rlnMaxValueProbDistribution'].astype(float).min(),2)
        print('Max probability particle range: '+str(maxProbMin)+'-'+str(maxProbMax))
        analyseProcessed.maxProbMax = maxProbMax
        analyseProcessed.maxProbMin = maxProbMin

        # Replace infinite values in the class2D with np.nan so classes with infinite resolution are not counted in max range
        # This means some particles will be assigned np.nan as their class resolution, is that what you want?
        class2DPtclsDf.replace([np.inf, -np.inf], np.nan, inplace=True)
        #class2DPtclsDf.dropna(subset=["rlnEstimatedResolutionClass2D"], how="all", inplace=True)

        # Find rlnEstimatedResolution based
        class2dresMax = class2DPtclsDf['rlnEstimatedResolutionClass2D'].astype(float).max()
        class2dresMin = class2DPtclsDf['rlnEstimatedResolutionClass2D'].astype(float).min()
        class2dresMean = class2DPtclsDf['rlnEstimatedResolutionClass2D'].astype(float).mean()
        class2dresMode = class2DPtclsDf['rlnEstimatedResolutionClass2D'].astype(float).mode()
        print('Resolution range of 2D classes: '+str(class2dresMin)+'-'+str(class2dresMax))
        analyseProcessed.class2dresMax = round(class2dresMax, 2)
        analyseProcessed.class2dresMin = round(class2dresMin, 2)
        analyseProcessed.class2dresMean = round(class2dresMean)
        analyseProcessed.class2dresMode = round(class2dresMode[0])

        # Return dataframe of class2dPtcls for use elsewhere
        # DEV DEV DEV should this be appended to the main particle dataframe?

        # Group by micrograph name and get average class2D resolution for particles on respective micrographs
        particlesByMic = class2DPtclsDf.groupby('rlnMicrographName')['rlnEstimatedResolutionClass2D'].aggregate(['mean']).round(2)
        particlesByMic = particlesByMic.reset_index()
        particlesByMic.rename(columns = {'mean':'Class2dRes_mean'}, inplace = True)
        # Convert MicrographName in base data acquisition name
        particlesByMic = addDataAcquisitionColumn(particlesByMic)
        # If working with a star file
        # Drop rlnMicrographName from star file dataframe as other star file merges will duplicate this column
        particlesByMic.drop('rlnMicrographName', inplace=True, axis=1)

        # Verify data structure
        #verifyDataStructure(search_mics.sessionStructure)
        #input('Press Enter to continue')
        # Now merge dataframes on basename column and search_mics.sessionStructure to get pick counts into main dataframe
        search_mics.sessionStructure = search_mics.sessionStructure.merge(particlesByMic, how = 'inner', on = ['DataAcquisition_name'])

        # Group by micrograph name and get modal class2D resolution for particles on respective micrographs
        #particlesByMic = class2DPtclsDf.groupby('rlnMicrographName')['rlnEstimatedResolutionClass2D'].agg(lambda x: pd.Series.mode(x)[0])
        #particlesByMic = particlesByMic.reset_index()
        #particlesByMic.rename(columns = {'rlnEstimatedResolutionClass2D':'Class2dRes_mode'}, inplace = True)
        # Convert MicrographName in base data acquisition name
        #particlesByMic = addDataAcquisitionColumn(particlesByMic)
        # If working with a star file
        # Drop rlnMicrographName from star file dataframe as other star file merges will duplicate this column
        #particlesByMic.drop('rlnMicrographName', inplace=True, axis=1)

        # Verify data structure
        #verifyDataStructure(search_mics.sessionStructure)
        #input('Press Enter to continue')
        # Now merge dataframes on basename column and search_mics.sessionStructure to get pick counts into main dataframe
        #search_mics.sessionStructure = search_mics.sessionStructure.merge(particlesByMic, how = 'inner', on = ['DataAcquisition_name'])

        # Plot 2D assignments per micrograph for this session
        plotPerMic(search_mics.sessionStructure, 'Class2dRes_mean', None, 40, 'pink', 'purple', 'Average Class2D resolution per micrograph (Angstrom)', 'Micrograph count', 'Average Class2D resolution', '_class2D_mean')
        #plotPerMic(search_mics.sessionStructure, 'Class2dRes_mode', None, 40, 'pink', 'purple', 'Modal Class2D resolution per micrograph (Angstrom)', 'Micrograph count', 'Modal Class2D resolution', '_class2D_mode')

    # Analyse Class3D in dev
    if os.path.isdir(analyseProcessed.class3dPath):
        analyseProcessed.class3D = 'Y'
    else:
        analyseProcessed.class3D = 'N'

    print('')

def plotClass2Dptcls(df, column):
    print('Plotting particles with Class2D assignments of: '+str(column))
    # Apply color map (cmap) to values in 'column' and return color values in df['cmap_color']
    #column = 'rlnEstimatedResolutionClass2D'
    cmapname = 'coolwarm'
    class2DPtclsDf_subset = applyCmapDf(df, column, cmapname)

    # Min max values for all mics in subset to be used for scale bar in plotitng
    min_val = class2DPtclsDf_subset[column].min()
    max_val = class2DPtclsDf_subset[column].max()

    #### random0
    name = search_mics.random0Name
    path = search_mics.randomName0Path
    # Filter dataframe by string
    df = class2DPtclsDf_subset[class2DPtclsDf_subset['rlnMicrographName'].str.contains(str(name))]

    if df.shape[0] > 0:
        micname = df['rlnMicrographName'].iloc[0]
        micpath = str(path)+'.jpg'
        #analyseProcessed_plot_colored_coordinates(df, analyseProcessed.ptclDPx, 0.15, 0.7, 'black', micpath, micname, True, column, cmapname, min_val, max_val)
        analyseProcessed_plot_colored_coordinates(df, analyseProcessed.ptclDPx, 0.66, 1, 'black', micpath, micname, False, column, cmapname, min_val, max_val)

    #### random1
    name = search_mics.random1Name
    path = search_mics.randomName1Path
    # Filter dataframe by string
    df = class2DPtclsDf_subset[class2DPtclsDf_subset['rlnMicrographName'].str.contains(str(name))]

    if df.shape[0] > 0:
        micname = df['rlnMicrographName'].iloc[0]
        micpath = str(path)+'.jpg'
        #analyseProcessed_plot_colored_coordinates(df, analyseProcessed.ptclDPx, 0.15, 0.7, 'black', micpath, micname, True, column, cmapname, min_val, max_val)
        analyseProcessed_plot_colored_coordinates(df, analyseProcessed.ptclDPx, 0.66, 1, 'black', micpath, micname, False, column, cmapname, min_val, max_val)

    #### random2
    name = search_mics.random2Name
    path = search_mics.randomName2Path
    # Filter dataframe by string
    df = class2DPtclsDf_subset[class2DPtclsDf_subset['rlnMicrographName'].str.contains(str(name))]

    if df.shape[0] > 0:
        micname = df['rlnMicrographName'].iloc[0]
        micpath = str(path)+'.jpg'
        #analyseProcessed_plot_colored_coordinates(df, analyseProcessed.ptclDPx, 0.15, 0.7, 'black', micpath, micname, True, column, cmapname, min_val, max_val)
        analyseProcessed_plot_colored_coordinates(df, analyseProcessed.ptclDPx, 0.66, 1, 'black', micpath, micname, False, column, cmapname, min_val, max_val)

    #### random3
    name = search_mics.random3Name
    path = search_mics.randomName3Path
    # Filter dataframe by string
    df = class2DPtclsDf_subset[class2DPtclsDf_subset['rlnMicrographName'].str.contains(str(name))]

    if df.shape[0] > 0:
        micname = df['rlnMicrographName'].iloc[0]
        micpath = str(path)+'.jpg'
        #analyseProcessed_plot_colored_coordinates(df, analyseProcessed.ptclDPx, 0.15, 0.7, 'black', micpath, micname, True, column, cmapname, min_val, max_val)
        analyseProcessed_plot_colored_coordinates(df, analyseProcessed.ptclDPx, 0.66, 1, 'black', micpath, micname, False, column, cmapname, min_val, max_val)

def assignClass2dParam2Ptcl(class2djobs, modelColumn, newColumnName):
    # Initialise dataframe for particles
    result_df = None

    print('Assigning Class2D model label to individual particles: '+str(modelColumn))
    for dir in tqdm(class2djobs):

        # Find _data and _model
        class2dstars = sorted(glob.glob(dir+'/run_it*_data.star', recursive = True))
        class2dmodels = sorted(glob.glob(dir+'/run_it*_model.star', recursive = True))

        if class2dstars and class2dmodels:
            class2dstar = class2dstars[-1]
            class2dmodel = class2dmodels[-1]

            # Read data from files into dataframes
            dataDf = starToDataFrame(class2dstar, 'particles') # Alistair Burt star parser
            modelDf = starToDataFrame(class2dmodel, 'model_classes') # Alistair Burt star parser

            # Create an index/class number column in the second dataframe
            modelDf['class'] = modelDf.index

            # Create a new column in the class2D particles dataframe
            dataDf[newColumnName] = None

            ## Assign the resolution of the 2D class a particle belongs to, to that particle
            # New vectorised code
            # Create a dictionary mapping 'rlnClassNumber' to 'rlnEstimatedResolution'
            resolution_mapping = modelDf.set_index('class')[modelColumn].to_dict()
            # Use the map function to create a new column 'rlnEstimatedResolutionClass2D' based on 'rlnClassNumber'
            dataDf[newColumnName] = dataDf['rlnClassNumber'].map(resolution_mapping)

            # Old loop code
            # Iterate through rows of the particles and model dataframes
            #for index, row in dataDf.iterrows():
            #    matching_row = modelDf.loc[modelDf['class'] == row['rlnClassNumber']]
            #    if not matching_row.empty:
            #        dataDf.at[index, "rlnEstimatedResolutionClass2D"] = matching_row.iloc[0]["rlnEstimatedResolution"]

            # Append the updated first dataframe to the result dataframe
            if result_df is None:
                result_df = dataDf
            else:
                result_df = pd.concat([result_df, dataDf], ignore_index=True)

    print()
    class2DPtclsDf = result_df

    return class2DPtclsDf

def analyseProcessed_plot_colored_coordinates(df, marker_size, center_alpha, edge_alpha, edge_color, background_image_path, micname, transparent, column, cmapname, min_val, max_val):

    print('Plotting Class2D on particles for: '+str(column))
    print(str(micname))
    print()

    # Values are pulled as strings and need to be numeric
    df_plot = df
    df_plot["rlnCoordinateX"] = pd.to_numeric(df_plot["rlnCoordinateX"])
    df_plot["rlnCoordinateY"] = pd.to_numeric(df_plot["rlnCoordinateY"])

    # Make a plot of the hole, using the measured hole diameter to input radius of circle
    # https://www.codespeedy.com/how-to-draw-shapes-in-matplotlib-with-python/
    fig2 = plt.figure(2)
    ax1 = fig2.add_subplot()

    # Dimensions for plot based on mic size
    plotSizeX = 5
    plotSizeY = (plotSizeX/(float(getCameraDim.cameraXpx)/float(getCameraDim.cameraYpx)))+1
    plt.gcf().set_size_inches(plotSizeX, plotSizeY)

    # Iterate through coordinates in dataframe and plot as circles on plot, color by nearest neighbour label close or sparse
    for index, row in df_plot.iterrows():
        circleFace = plt.Circle(
            (row['rlnCoordinateX'], row['rlnCoordinateY']),
            radius=marker_size/2,
            facecolor=row['cmap_color'])
            #alpha=center_alpha) # Some cmap colors are assinged as black but alpha 0, this face alpha will override that and show them up
        ax1 = plt.gca()
        ax1.add_patch(circleFace)

        circleEdge = plt.Circle(
            (row['rlnCoordinateX'], row['rlnCoordinateY']),
            radius=marker_size/2,
            facecolor='none',
            edgecolor=edge_color,
            alpha=edge_alpha,
            linewidth=0.5)
        ax1 = plt.gca()
        ax1.add_patch(circleEdge)

    # Add a rectangle patch to create a border
    rect = Rectangle((0, 0), float(getCameraDim.cameraXpx), float(getCameraDim.cameraYpx), linewidth=1, edgecolor='black', facecolor='none', alpha=0.5)
    #rect = Rectangle((df_plot['rlnCoordinateX'].min(), df_plot['rlnCoordinateY'].min()), df_plot['rlnCoordinateX'].max() - df_plot['rlnCoordinateX'].min(), df_plot['rlnCoordinateY'].max() - df_plot['rlnCoordinateY'].min(), linewidth=1, edgecolor='black', facecolor='none', alpha=0.5)
    plt.gca().add_patch(rect)

    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)

    plt.xlim(0,float(getCameraDim.cameraXpx))
    plt.ylim(0,float(getCameraDim.cameraYpx))

    ## Create a color map legend
    # Find min and max to assign these to min and max of cmap
    #min_val = df[column].min()
    #max_val = df[column].max()
    sm = plt.cm.ScalarMappable(cmap=cmapname, norm=plt.Normalize(vmin=min_val, vmax=max_val))
    sm._A = []  # Empty array to fool the colorbar into thinking there's data
    # Create a divider for existing axes instance
    divider = make_axes_locatable(ax1)
    # Define the size and position of the colorbar legend
    cax = divider.append_axes("bottom", size="5%", pad=0.5)  # Adjust 'size' and 'pad' as needed

    # Add colorbar legend
    cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label(str('Assigned 2D class property to particle '+str(column)), size=6)  # Set the label for the colorbar

    # Set the font size of the colorbar tick labels
    cbar.ax.tick_params(labelsize=6)

    #ax1.axis('scaled')
    plt.tight_layout()
    filename = os.path.basename(Path(micname))
    tmp1 = filename.replace('_fractions.mrc','')
    name = tmp1.replace('_EER.mrc','')
    #plt.title('Nearest neighbour clustering analysis')
    #plt.show()
    if transparent == True:
        fig2.savefig(main.plot_dir+'/particles/'+str(name)+'_'+str(column)+'_trans.png', dpi=300, transparent=True, bbox_inches='tight')
    else:
        fig2.savefig(main.plot_dir+'/particles/'+str(name)+'_'+str(column)+'.png', dpi=300, transparent=False, bbox_inches='tight')
    plt.clf()
    plt.figure(2).clear()
    plt.close('all')

def processed_plotMotion(df):

    print('Plotting Motion...')
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = 'rlnAccumMotionTotal'
    filename = str(xml_session.sessionName)+'_motion_total'
    xaxis = 'Accumulated motion (Angstroms)'
    title = 'Accumulated motion (Angstroms) with total '+ str(calc_params.doseTotal)+' e/A^2'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    #print(filtered)
    # Quartile outliers
    #q_low = df_plot[col].quantile(0.01)
    #q_hi  = df_plot[col].quantile(0.99)
    #df_plot = df_plot[(df_plot[col] < q_hi) & (df_plot[col] > q_low)]
    # Need to implement subplotting to show totalmotion, initialmotion and latemotion
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='green',
            binrange=(0,30),
            edgecolor="navy", linewidth=0.5).set(title=title)

    ax10.set(xlabel=xaxis, ylabel='Micrograph count')
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = 'rlnAccumMotionEarly'
    filename = str(xml_session.sessionName)+'_motion_early'
    xaxis = 'Early motion (Angstroms)'
    title = 'Early motion (Angstroms) with total '+ str(calc_params.doseTotal)+' e/A^2'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    #print(filtered)
    # Quartile outliers
    #q_low = df_plot[col].quantile(0.01)
    #q_hi  = df_plot[col].quantile(0.99)
    #df_plot = df_plot[(df_plot[col] < q_hi) & (df_plot[col] > q_low)]
    # Need to implement subplotting to show totalmotion, initialmotion and latemotion
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='tomato',
            binrange=(0,30),
            edgecolor="navy", linewidth=0.5).set(title=title)

    ax10.set(xlabel=xaxis, ylabel='Micrograph count')
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = 'rlnAccumMotionLate'
    filename = str(xml_session.sessionName)+'_motion_late'
    xaxis = 'Late motion (Angstroms)'
    title = 'Late motion (Angstroms) with total '+ str(calc_params.doseTotal)+' e/A^2'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    #print(filtered)
    # Quartile outliers
    #q_low = df_plot[col].quantile(0.01)
    #q_hi  = df_plot[col].quantile(0.99)
    #df_plot = df_plot[(df_plot[col] < q_hi) & (df_plot[col] > q_low)]
    # Need to implement subplotting to show totalmotion, initialmotion and latemotion
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='royalblue',
            binrange=(0,30),
            edgecolor="navy", linewidth=0.5).set(title=title)

    ax10.set(xlabel=xaxis, ylabel='Micrograph count')
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    # Multi motion plot
    col0 = 'rlnAccumMotionEarly'
    col1 = 'rlnAccumMotionLate'
    col2 = 'rlnAccumMotionTotal'
    xaxis2 = 'Total'
    xaxis0 = 'Early'
    xaxis1 = 'Late'
    df_plot[col2] = pd.to_numeric(df_plot[col2])
    df_plot[col0] = pd.to_numeric(df_plot[col0])
    df_plot[col1] = pd.to_numeric(df_plot[col1])
    filtered2 = df_plot[(np.abs(stats.zscore(df_plot[col2])) < 3.0)]
    filtered0 = df_plot[(np.abs(stats.zscore(df_plot[col0])) < 3.0)]
    filtered1 = df_plot[(np.abs(stats.zscore(df_plot[col1])) < 3.0)]

    figure, axes = plt.subplots(3, 1, sharex=True, figsize=(6, 8))
    figure.suptitle('Accumulated motion (Angstroms) with total '+ str(calc_params.doseTotal)+' e/A^2')
    axes[0].set_title(xaxis0)
    axes[1].set_title(xaxis1)
    axes[2].set_title(xaxis2)

    axes[0].set_ylabel('Micrographs')
    axes[1].set_ylabel('Micrographs')
    axes[2].set_ylabel('Micrographs')

    axes[2].set_xlabel('Angstroms')

    sns.histplot(data=filtered0, x=col0, kde='False',
            bins=40, color='tomato',
            binrange=(0,30),
            edgecolor="navy", linewidth=0.5, ax = axes[0])
    sns.histplot(data=filtered1, x=col1, kde='False',
            bins=40, color='royalblue',
            binrange=(0,30),
            edgecolor="navy", linewidth=0.5, ax = axes[1])
    sns.histplot(data=filtered2, x=col2, kde='False',
            bins=40, color='green',
            binrange=(0,30),
            edgecolor="navy", linewidth=0.5, ax = axes[2])

    figure.savefig(str(main.plot_dir)+'/'+str(filename)+'_multi.png', dpi=300)

def plotPerMic(df, column, filter, bins, fill, edge, title, yaxis, xaxis, filename):
    # plotPerMic(df, 'ptcls_per_mic', 40, deeppink, darkviolet, 'Particles per mic ('+str(GridSquare[0])+')', Micrograph count, Particle per micrograph, '_ptcls_per_mic_'+str(GridSquare[0]))
    print('Plotting histogram analysis of column: '+str(column))
    #print()
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = column
    suffix = str(xml_session.sessionName)+str(filename)
    yaxis = 'Micrograph count'
    xaxis = str(xaxis)
    title = str(title)
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    if filter == "filter":
        filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    elif filter == "nofilter":
        filtered = df_plot
    else:
        filtered = df_plot
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=bins, color=fill,
            edgecolor=edge, linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(suffix)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

def plotPtclsPerMic(df, column, title, name):
    #print('Plotting picking analysis')
    #print()
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = column
    filename = str(xml_session.sessionName)+str(name)
    yaxis = 'Micrograph count'
    xaxis = 'Particles per micrograph'
    title = str(title)
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='deeppink',
            edgecolor="darkviolet", linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

# DEV DEV DEV These plotting per mic functions could be combined into a single function
def plotCtfPerMic(df, column, title, name):
    #print('Plotting picking analysis')
    #print()
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = column
    filename = str(xml_session.sessionName)+str(name)
    yaxis = 'Micrograph count'
    xaxis = 'Ctf per micrograph'
    title = str(title)
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='deeppink',
            edgecolor="darkviolet", linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)
    plt.close(10)

def plotCorPerMic(df, column, title, name):
    #print('Plotting picking analysis')
    #print()
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = column
    filename = str(xml_session.sessionName)+str(name)
    yaxis = 'Micrograph count'
    xaxis = 'Motion per micrograph'
    title = str(title)
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='deeppink',
            edgecolor="darkviolet", linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

def processed_plotPtclsPerMic(df):
    print('Plotting picking analysis')
    print()
    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = 'ptcls_per_mic'
    filename = str(xml_session.sessionName)+'_ptcls_per_mic'
    yaxis = 'Micrograph count'
    xaxis = 'Particles per micrograph'
    title = 'Particles per micrograph'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='deeppink',
            edgecolor="darkviolet", linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib plots DEV DEV DEV
    # Data DEV DEV DEV
    col = 'ideal_picking'
    filename = str(xml_session.sessionName)+'_ptcls_per_mic_normalised'
    yaxis = 'Micrograph count'
    xaxis = 'Idealised particle number per micrograph'
    title = 'Picks normalised to idealised density for particle of '+str(analyseProcessed.ptclDAng)+' Ang'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    #Filtering to remove motion outliers
    # DEV DEV DEV
    # Put this into a function for general filtering of dataframe outliers
    # It should report what the range in, and range after filtering, how many data points were removed
    filtered = df_plot[(np.abs(stats.zscore(df_plot[col])) < 3.0)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=filtered, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='violet',
            edgecolor="indigo", linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(main.plot_dir)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)
    # Necessary to clear, otherwise other plots pick up this one
    plt.clf()

def processed_output():
    # Do some date conversion for reporting
    sessionDate = screeningSession_xml.date.strftime("%Y-%m-%d")

    dictHorizontal = [
    {'Proposal': lookupBAG.proposal,
     'Visit': main.visitNameDLS,
     'BAG': lookupBAG.name,
     'Institute': lookupBAG.institute,
     'Anonymous': lookupBAG.anonymous,
     'Session_date': sessionDate,
     'Directory': main.visitPath,
     'Error': 'None',
     'Advisory': assessSession.advisoryCount,
     'EPU session': xml_session.sessionName,
     'EPU session date': xml_session.sessionDate,
     'Atlas': 'Sample'+str(round(xml_session.autoSlot)),
     'Instrument': getScopeNameMag.scope,
     'AFIS/Accurate': xml_session.afisMode,
     'Grid_type': xml_session.gridType,
     'I0 set': xml_session.I0set,
     'I0 min int': xml_session.I0MinInt,
     'I0 max int': xml_session.I0MaxInt,
     'Hole (um)': xml_targets.holeSizeRegular,
     'Space (um)':xml_targets.holeSpaceRegular,
     'Beam (um)': xml_presets.beamD,
     'Shots per hole': xml_template_dataframe.shotsPerHole,
     'ShotID': xml_template_dataframe.templateNameList,
     'Movie_dir': findRaw.dir,
     'Movie_format': xml_session.doseFractionOutputFormat,
     'Movie_ext': count_movies.format,
     'Relion': searchProcessed.path,
     'Jobs': searchProcessed.number ,
     'Processed': searchProcessed.ispyb,
     'apix': getScopeNameMag.apix,
     'Targeted_squares': xml_targets.targetSquares,
     'Collected_squares': xml_targets.collectedSquares,
     'Average Foils per Square': xml_targets.avSqFoilHoleNo,
     'Tilt': xml_presets_data.stageAlpha,
     'Total_EPU_mics': count_mics.number,
     'Total_movies': count_movies.number,
     'MotionCorr': analyseProcessed.motionCorNo,

     'Early_motion_mean': analyseProcessed.earlyMotionMean,
     'Early_motion_min': analyseProcessed.earlyMotionMin,
     'Early_motion_max': analyseProcessed.earlyMotionMax,
     'Late_motion_mean': analyseProcessed.lateMotionMean,
     'Late_motion_min': analyseProcessed.lateMotionMin,
     'Late_motion_max': analyseProcessed.lateMotionMax,
     'Total_motion_mean': analyseProcessed.totMotionMean,
     'Total_motion_min': analyseProcessed.totMotionMin,
     'Total_motion_max': analyseProcessed.totMotionMax,

     'First_mic': str(countMotionCor.firstCollectedTime),
     'First_cor': str(countMotionCor.firstCorrectedTime),
     'Last_mic': str(countMotionCor.lastCollectedTime),
     'Last_cor': str(countMotionCor.lastCorrectedTime),

     'CtfFind': analyseProcessed.ctfFindNo,
     'Ctf_min_A': analyseProcessed.ctfMin,
     'Ctf_max_A': analyseProcessed.ctfMax,
     'Picking': analyseProcessed.picking,
     'Pick_min': analyseProcessed.pickMin,
     'Pick_max': analyseProcessed.pickMax,
     'Ptcl_d_A': analyseProcessed.ptclDAng,
     'Particles': analyseProcessed.particleNo,
     'Extracted_mics': analyseProcessed.particleMicNo,

     'Box_px': analyseProcessed.extractPx,
     'Box_ang': analyseProcessed.extractAng,
     'Av_ptcls_per_mic': analyseProcessed.particlesPerMicAve,
     'Std_ptcls_per_mic': analyseProcessed.particlesPerMicStd,
     'Norm_av_ptcls_per_mic': analyseProcessed.particlesPerMicAveNorm,
     'Norm_std_ptcls_per_mic': analyseProcessed.particlesPerMicStdNorm,
     'DBSCAN_eps:': eps_opt.eps,
     'Mean_ptcls_pct_clustered': analyseProcessed.particlesPerMicClusterAve,
     'Std_ptcls_pct_clustered': analyseProcessed.particlesPerMicClusterStd,
     'Class2D_requested': analyseProcessed.request2D,
     'Class2D_results': analyseProcessed.class2D,
     'Mask_ang': analyseProcessed.ptclDAng,
     'Class2D_prob_max': analyseProcessed.maxProbMax,
     'Class2D_prob_min': analyseProcessed.maxProbMin,
     'Class2D_res_max': analyseProcessed.class2dresMax,
     'Class2D_res_min': analyseProcessed.class2dresMin,
     'Class2D_res_mean': analyseProcessed.class2dresMean,
     'Class2D_res_mode': analyseProcessed.class2dresMode,
     'Class3D_requested': analyseProcessed.request3D,
     'Class3D_results': analyseProcessed.class3D,
     'Report path': main.pdf_dir,
     'Session report': main.report_name,
     'Processed report': main.report_processed_name,
     'epuVersion': xml_session.epuVersion}
    ]

    df = pd.DataFrame.from_dict(dictHorizontal)
    df.to_csv (main.csv_dir+'/'+xml_session.sessionName+'_processed.csv', index = False, header=True)

# This function will check for strange things in the session and create/add to an advisory report
def assessSession():
   print('Assessing session setup for common errors (advisories)')

   # This will keep track of how many warnings are generated for the session
   assessSession.advisoryCount = 0

   ## The following are the advisory checks
   # Frame dose
   if type(calc_params.doseFrameTotal) == float:
      threshold = 1.3
      if calc_params.doseFrameTotal > threshold:
         print('Frame dose >'+str(threshold)+' e/A^2: '+str(calc_params.doseFrameTotal))
         print()
         advisory('dose','Frame dose of '+str(calc_params.doseFrameTotal)+' ('+str(calc_params.doseMessage)+'), is greater than '+str(threshold)+' e/A^2')

   # Report advisory number
   print('Number of advisories for this session: '+str(assessSession.advisoryCount))
   print()

# This function will add advisories to a dataframe which can later be written out
def advisory(type, advice):

   # Do some date conversion for reporting
   sessionDate = screeningSession_xml.date.strftime("%Y-%m-%d")

   # Advance advisory count
   assessSession.advisoryCount = assessSession.advisoryCount+1

   # Add advisory into dataframe for reporting
   try:
      # Look for exisitng advisory dataframe
      print(advisory.df.head())
   # create advisory.df when it is None
   except AttributeError:
      data = [[lookupBAG.proposal, main.visitNameDLS, main.sessionTypeDLS, \
               lookupBAG.name, lookupBAG.institute, sessionDate, main.visitPath, \
               xml_session.sessionName, xml_session.sessionDate, getScopeNameMag.scope, \
               type, advice, main.pdf_dir, main.report_name]]
      advisory.df = pd.DataFrame(data, columns=['Proposal', 'Visit', 'Session type', \
               'BAG', 'Institute', 'Session date', 'Directory', 'EPU Session', \
               'EPU session date', 'Instrument', 'Advisory_type', 'Advisory', 'Report path', 'Session report'])
   # create advisory.df when it hasn't even been defined
   except NameError:
      data = [[lookupBAG.proposal, main.visitNameDLS, main.sessionTypeDLS, \
               lookupBAG.name, lookupBAG.institute, sessionDate, main.visitPath, \
               xml_session.sessionName, xml_session.sessionDate, getScopeNameMag.scope, \
               type, advice, main.pdf_dir, main.report_name]]
      advisory.df = pd.DataFrame(data, columns=['Proposal', 'Visit', 'Session type', \
               'BAG', 'Institute', 'Session date', 'Directory', 'EPU Session', \
               'EPU session date', 'Instrument', 'Advisory_type', 'Advisory', 'Report path', 'Session report'])
   # Otherwise append to existing advisory.df dataframe
   else:
      data = [[lookupBAG.proposal, main.visitNameDLS, main.sessionTypeDLS, \
               lookupBAG.name, lookupBAG.institute, sessionDate, main.visitPath, \
               xml_session.sessionName, xml_session.sessionDate, getScopeNameMag.scope, \
               type, advice, main.pdf_dir, main.report_name]]
      df2 = pd.DataFrame(data, columns=['Proposal', 'Visit', 'Session type', \
               'BAG', 'Institute', 'Session date', 'Directory', 'EPU Session', \
               'EPU session date', 'Instrument', 'Advisory_type', 'Advisory', 'Report path', 'Session report'])
      advisory.df.append(df2, ignore_index = True)

def writeAdvisory():

   try:
      # Look for exisitng advisory dataframe
      print(advisory.df.head())
   # catch when advisory.df is None
   except AttributeError:
      pass
   # catch when it hasn't even been defined
   except NameError:
      pass
   # Otherwise write existing dataframe to file
   else:
      advisory.df.to_csv (main.csv_dir+'/'+xml_session.sessionName+'_advisory.csv', index = False, header=True)

def countMotionCor(motioncorrPath):
        print('')
        # DEV DEV DEV Use motioncocrrected MRC timestamp, can you compare to file time of micrograph acquisition?
        print('Finding motion corrected micrographs (known to be slow)')
        corMics = glob.glob(motioncorrPath+'/Movies/**/*mrc', recursive = True)
        sortedFiles = sorted(corMics, key = lambda file: os.path.getmtime(file))
        fileTimes = []
        fileNumber = []
        # Get data acquisition time by file stamp - fast but system timestamp dependent
        j = 1
        for file in tqdm(sortedFiles):
            # Get motion corrected file write time by file timestamp
            fileTimeStamp = datetime.datetime.fromtimestamp(os.path.getmtime(file))
            fileTime = datetime.datetime.strftime(fileTimeStamp, '%Y-%m-%d %H:%M:%S')
            fileTimes.append(fileTime)
            fileNumber.append(j)
            j = j+1

        # Pull in the micrograph timestamps
        df = plotExposures.micTimes
        countMotionCor.firstCollectedTime = df['Time'].iloc[0]
        countMotionCor.lastCollectedTime = df['Time'].iloc[-1]
        #print(df)
        #df.to_csv(main.plot_dir+'/'+xml_session.sessionName+'_motioncorr_timings_df_mics.csv', index = False, header=True)

        # Pull in the corrected micrograph timestamps
        data = {'Count': fileNumber, 'corTime': fileTimes}
        df2 = pd.DataFrame(data)
        countMotionCor.micTimes = df2
        countMotionCor.firstCorrectedTime = df2['corTime'].iloc[0]
        countMotionCor.lastCorrectedTime = df2['corTime'].iloc[-1]
        #print(df2)
        #df2.to_csv(main.plot_dir+'/'+xml_session.sessionName+'_motioncorr_timings_df_cor2.csv', index = False, header=True)

        # Report micrograph and corrected timestamps
        print('')
        print('Number of micrograph times (EPU directory): '+str(len(df.index)))
        print('Number of raw movies (raw* directory): '+str(count_movies.number))
        print('Format of raw movies (raw* directory): '+str(count_movies.format))
        print('Number of corrected times (Relion directory): '+str(len(df2.index)))
        print('')
        print('Time of first collected micrograph in EPU directory: '+str(countMotionCor.firstCollectedTime))
        print('Time of first corrected micrograph in Relion directory: '+str(countMotionCor.firstCorrectedTime))
        print('')
        print('Time of final collected micrograph in EPU directory: '+str(countMotionCor.lastCollectedTime))
        print('Time of final corrected micrograph in Relion directory: '+str(countMotionCor.lastCorrectedTime))

def plotMotionCorTimings(motioncorrPath):

        df = plotExposures.micTimes
        df2 = countMotionCor.micTimes

        #corTimes = df2['corTime']
        df3 = pd.merge(pd.DataFrame(df), pd.DataFrame(df2), left_on=['Count'],
             right_on= ['Count'], how='left')
        #df3.to_csv(main.plot_dir+'/'+xml_session.sessionName+'_motioncorr_timings_df_merge-filt.csv', index = False, header=True)
        #print(df3)

        #fig98 = plt.figure(figsize=(8,4))
        #ax98 = fig98.add_subplot(111)
        #ax98 = df.plot.scatter(x="Time" , y="Count", s=25, figsize=(8,4), color='blue', label="Mic", marker='s')
        #ax99 = df.plot.scatter(x="corTime" , y="Count", s=20, figsize=(8,4), color='none', edgecolors='cyan', label="Corrected", marker='o', ax=ax98)
        #ax98.plot(df['Time'], df['Count'])
        #ax98.plot(df2['corTime'], df['Count'])
        #ax98.axes.xaxis.set_visible(False)
        #ax98.set_ylabel("Micrograph count")
        #plt.tight_layout()
        #fig98 = ax98.get_figure()
        #fig98.savefig(main.plot_dir+'/'+xml_session.sessionName+'_motioncorr_timings.png', dpi=300)
        #plt.figure(98).clear()
        #plt.close(98)

        #fig99,ax = plt.subplots()
        #ax.axes.xaxis.set_visible(False)
        #df3.plot.scatter(x="Time" , y="Count", s=30, figsize=(8,4), color='blue', label="Acquisition", marker='.', ax=ax)
        #df3.plot.scatter(x="corTime" , y="Count", s=30, figsize=(8,4), color='orange', label="MotionCorr", marker='.', ax=ax)
        #ax.set_xlabel("Time")
        #ax.set_ylabel("Micrograph count")
        #plt.xlim([0, str(len(df.index))])
        #plt.tight_layout()
        #fig99.savefig(main.plot_dir+'/'+xml_session.sessionName+'_motioncorr_timings.png', dpi=300)
        #plt.figure(99).clear()
        #plt.close(99)

        # This is working but I think something is wrong with it, compare to the _exposures plot, there is a difference
        ax0 = df.plot.scatter(x="Time" , y="Count", figsize=(7,4), color='blue', label="Acquisition")
        df3.plot.scatter(x="corTime" , y="Count", figsize=(7,4), color='white', edgecolors='orange', label="Corrected", marker='o', ax=ax0)
        #df['Time'] = pd.to_datetime(df['Time'])
        #df['Time'] = df['Time'].dt.strftime("%Y-%m-%d %H:%M:%S")

        #import matplotlib.dates as mdates
        #myFmt = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
        #ax0.xaxis.set_major_formatter(myFmt)

        plt.title("Acquisition time")
        plt.xlabel("Time")
        ax0.set_ylabel("Cumulative exposures & corrections", color='black')
        ax0.xaxis.set_major_locator(plt.MaxNLocator(10))
        plt.xticks( rotation=25 )
        plt.xlim([countMotionCor.firstCollectedTime,countMotionCor.lastCollectedTime])
        #ax0.set_xticklabels(labels=df['Time'], rotation=25, rotation_mode="anchor", ha="right", fontsize=6)

        plt.tight_layout()
        fig1 = ax0.get_figure()
        fig1.savefig(main.plot_dir+'/'+xml_session.sessionName+'_motioncorr_timings.png', dpi=300)
        plt.figure(1).clear()
        plt.close()
        plt.cla()
        plt.clf()

def countProcessed(path, ext):
    print('Scanning processed files '+ext+' in '+path)
    # This will search all subdirectories, does not care about directort structure
    # Sometimes very slow, abandon for now and just read from star files
    list = glob.glob(path+'/**/*'+ext, recursive = True)
    print('Scan complete...')
    return len(list)

def get_username():
    return pwd.getpwuid(os.getuid())[0]

def find_mics(path, search):
    # Need to have an independent function to find the mics, then move into search_mics to sort them out
    # So find mics can be used independently

    print('Looking for micrograph data in EPU directory using extension: '+search)
    print('')
    # Just get file names
    #files = glob.glob("./Images-Disc1/GridSquare*/Data/*xml")
    #files.sort(key=os.path.getmtime)
    #print("\n".join(files))

    # Old method of finding xml files
    searchedFiles = glob.glob(path+"/**/GridSquare*/Data/*"+search+'*')
    #searchedFiles = glob.iglob(main.epu+"/Images-Disc1/GridSquare*/Data/*"+search)
    if searchedFiles:
       print('Found micrograph data: '+str(len(searchedFiles)))
    else:
       print('No micrographs found with search term: '+search)
       searchedFiles = 'exit'

    # New method of finding xml files
    #searchedFiles = searchSupervisorData.xmlList

    return searchedFiles

def count(searchedFiles):
    count.mic_count = len(searchedFiles)

def search_mics(searchedFiles, method):
    fileTimes = []
    fileNumber = []

    if method == 'F':
        print('Reading times from file timestamps (modification), will be fast but could be inaccurate')
        print('')
        # https://stackoverflow.com/questions/23430395/glob-search-files-in-date-order
        # Get all data acquisition xml file names with time stamps
        print('Sorting files by modifications time')
        sortedFiles = sorted( searchedFiles, key = lambda file: os.path.getmtime(file))
        search_mics.micList = sortedFiles

        print('Iterating through files and gathering sorted collection times')
        # Get data acquisition time by file stamp - fast but system timestamp dependent
        j = 1
        for file in tqdm(sortedFiles):
            # Get data acquisiton time by file timestamp
            fileTimeStamp = datetime.datetime.fromtimestamp(os.path.getmtime(file))
            fileTime = datetime.datetime.strftime(fileTimeStamp, '%Y-%m-%d %H:%M:%S')
            # Returned as <class 'datetime.datetime'>
            #print(type(fileTime))

            # Get data acquisition time by time recorded in xml, very slow for all files
            fileTimes.append(fileTime)
            fileNumber.append(j)
            j = j+1

    elif method == 'X':
        print('Reading file times from xml, will be accurate but likely slow due to I/O')
        print('')

        # Create a dataframe for aggregating GridSquare and Acquisition names
        search_mics.sessionStructure = pd.DataFrame()

        print('Iterating through files and gathering collection times')
        # Get data acquisition time by time written in xml file - more likely accurate
        j = 1
        start_time = time()

        # Precompile regular expressions
        grid_square_regex = re.compile(r"GridSquare")
        data_name_split_regex = re.compile(r'[^_]+')

        # Initialize j and prepare lists to collect data
        j = 0
        fileTimes = []
        fileNumber = []
        data_frames = []

        for file in tqdm(searchedFiles):
            try:
                with open(file) as xml:
                    fileTimeString = get_mic_date_xmltree(file)

                    fileTime = datetime.datetime.strptime(fileTimeString, '%Y-%m-%d %H:%M:%S')
                    fileTimes.append(fileTime)
                    fileNumber.append(j)
                    j += 1

                    # Get full path to DataAcquisition
                    path = os.path.normpath(file)
                    pathList = path.split(os.sep)

                    # GridSquare name and path
                    grid_square_match_path = [match for match in pathList if grid_square_regex.search(match)]
                    findPath = ['/'.join(group) for f, group in groupby(pathList, lambda x: x == grid_square_match_path[0]) if not f]
                    GridSquarePath = str(findPath[0] + '/' + grid_square_match_path[0])
                    GridSquareName = grid_square_match_path[0]

                    # GridSquare image
                    GridSquareImages = glob.glob(os.path.join(GridSquarePath, 'GridSquare_*.jpg'))
                    GridSquareImagePath = os.path.basename(GridSquareImages[0]) if GridSquareImages else ""

                    # DataAcquisition
                    DataPath = os.path.dirname(path)
                    DataName = os.path.splitext(os.path.basename(path))[0]
                    split_data_name = data_name_split_regex.findall(DataName)
                    shotID, imgshiftID = split_data_name[3:5]

                    temp = {
                        'Atlas_path': main.atlasdir,
                        'Atlas_img': main.atlasimg,
                        'GridSquare_path': GridSquarePath,
                        'GridSquare_name': GridSquareName,
                        'GridSquare_img': GridSquareImagePath,
                        'DataAcquisition_path': DataPath,
                        'DataAcquisition_name': DataName,
                        'shotID': shotID,
                        'imgshiftID': imgshiftID,
                        'TimeStamp': fileTimeString,
                    }
                    data_frames.append(temp)

            except IOError as err:
                continue

        # Combine all data frames into a single DataFrame
        search_mics.sessionStructure = pd.DataFrame(data_frames)
        print(search_mics.sessionStructure)

        end_time = time()
        elapsed_time = end_time - start_time
        print("Total time:", elapsed_time, "seconds")

        # Sort dataframe on datetime
        search_mics.sessionStructure['TimeStamp'] = pd.to_datetime(search_mics.sessionStructure['TimeStamp'])
        #print(search_mics.sessionStructure.dtypes)
        search_mics.sessionStructure.sort_values(by='TimeStamp', inplace = True)

        #verifyDataStructure(search_mics.sessionStructure)
        #input('Press Enter to continue')

        # The final Atlas, GridSquare, DataAcquisition data structure
        #print(search_mics.sessionStructure)
        #print(search_mics.sessionStructure['Atlas_path'].iloc[0])
        #print(search_mics.sessionStructure['Atlas_img'].iloc[0])
        #print(search_mics.sessionStructure['GridSquare_path'].iloc[0])
        #print(search_mics.sessionStructure['GridSquare_name'].iloc[0])
        #print(search_mics.sessionStructure['GridSquare_img'].iloc[0])
        #print(search_mics.sessionStructure['DataAcquisition_path'].iloc[0])
        #print(search_mics.sessionStructure['DataAcquisition_name'].iloc[0])
        #print(search_mics.sessionStructure['Atlas_path'].iloc[0])
        #print(search_mics.sessionStructure['Atlas_path'].iloc[0])
        #input('Press Enter to continue')

        print('Sorting files by times read from xml files')
        # The following sorting is definitely working because identifying first and last mic later is working
        sortedFiles = [x for _, x in sorted(zip(fileTimes, searchedFiles), key=lambda pair: pair[0])]
        search_mics.micList = sortedFiles
        # But then you need to get the datetime formatted fileTimes sorted for plotting
        fileTimes = sorted(fileTimes)

    # Message
    print('')
    print('Finished reading in files and sorting by time')

    # Get first and last file in DataAcquisition xml file paths
    # Remember, this is the path into the EPU directory
    first = sortedFiles[0]
    search_mics.firstXml = first
    last = sortedFiles[-1]
    search_mics.lastXml = last
    # DEV DEV DEV, it woould be better to create a dataframe with the random choices and then extract names and paths out of that dataframe
    # Get some random xml files of mics for display
    random = secrets.choice(sortedFiles)
    search_mics.random0Xml = random
    random1 = secrets.choice(sortedFiles)
    search_mics.random1Xml = random1
    random2 = secrets.choice(sortedFiles)
    search_mics.random2Xml = random2
    random3 = secrets.choice(sortedFiles)
    search_mics.random3Xml = random3

    # Transform to get first and last name, without path
    search_mics.firstName = os.path.basename(Path(first).with_suffix(''))
    search_mics.lastName = os.path.basename(Path(last).with_suffix(''))
    search_mics.random0Name = os.path.basename(Path(random).with_suffix(''))
    search_mics.random1Name = os.path.basename(Path(random1).with_suffix(''))
    search_mics.random2Name = os.path.basename(Path(random2).with_suffix(''))
    search_mics.random3Name = os.path.basename(Path(random3).with_suffix(''))

    search_mics.random0JPG = search_mics.random0Name+'.jpg'
    search_mics.random1JPG = search_mics.random1Name+'.jpg'
    search_mics.random2JPG = search_mics.random2Name+'.jpg'
    search_mics.random3JPG = search_mics.random3Name+'.jpg'

    search_mics.firstNamePath = Path(first).with_suffix('')
    search_mics.lastNamePath = Path(last).with_suffix('')
    search_mics.randomName0Path = Path(random).with_suffix('')
    search_mics.randomName1Path = Path(random1).with_suffix('')
    search_mics.randomName2Path = Path(random2).with_suffix('')
    search_mics.randomName3Path = Path(random3).with_suffix('')

    # Code to get square image and foil hole image for insertion
    squareDir = Path(search_mics.randomName0Path).parents[1]
    squarePath = glob.glob(str(squareDir)+'/GridSquare_*.jpg')
    search_mics.randomSquarePath = str(squarePath[0])
    search_mics.randomSquareName = os.path.basename(Path(squarePath[0]).with_suffix(''))
    search_mics.randomMicPath = str(search_mics.randomName0Path)+'.jpg'

    squareDir = Path(search_mics.randomName1Path).parents[1]
    squarePath = glob.glob(str(squareDir)+'/GridSquare_*.jpg')
    search_mics.randomSquare1Path = str(squarePath[0])
    search_mics.randomSquare1Name = os.path.basename(Path(squarePath[0]).with_suffix(''))
    search_mics.randomMic1Path = str(search_mics.randomName1Path)+'.jpg'

    squareDir = Path(search_mics.randomName2Path).parents[1]
    squarePath = glob.glob(str(squareDir)+'/GridSquare_*.jpg')
    search_mics.randomSquare2Path = str(squarePath[0])
    search_mics.randomSquare2Name = os.path.basename(Path(squarePath[0]).with_suffix(''))
    search_mics.randomMic2Path = str(search_mics.randomName2Path)+'.jpg'

    squareDir = Path(search_mics.randomName3Path).parents[1]
    squarePath = glob.glob(str(squareDir)+'/GridSquare_*.jpg')
    search_mics.randomSquare3Path = str(squarePath[0])
    search_mics.randomSquare3Name = os.path.basename(Path(squarePath[0]).with_suffix(''))
    search_mics.randomMic3Path = str(search_mics.randomName3Path)+'.jpg'

    print('Found first and last micrographs')

    # initialize data of lists.
    data = {'Count': fileNumber, 'Time': fileTimes}

    # Create DataFrame
    df = pd.DataFrame(data)

    #print('Saving micrograph times to csv')
    #csvfile = main.report_dir+'/'+xml_session.sessionName+'_times.csv'
    #df.to_csv(csvfile, index = False, header=True)
    #print('Saved micrograph times to csv')

    print('Plotting exposure time data to graph')

    # Can write out exposure times but IO is slow
    #plotExposures(csvfile)

    # Or plot from dataframe
    plotExposures(df)

def verifyDataStructure(df):
    rows = len(df)
    print('Data structure has this many rows: '+str(rows))

def plotExposures(csv):

    # Read in csv (I/O) into dataframe and then plot from dataframe into time plot
    #df = pd.read_csv(csv)

    # Alternatively plot directly from data frame, much faster than I/O
    df = csv

    # Spit into global dataframe for use elsewhere... is this a good idea?
    plotExposures.micTimes = df

    # Plot exposures
    # https://stackoverflow.com/questions/4090383/plotting-unix-timestamps-in-matplotlib
    fig1 = plt.figure(1)

    ax0 = df.plot.scatter(x="Time", y="Count", s=10, figsize=(6,3), color='white', edgecolors='blue', marker='o')
    df['Time'] = pd.to_datetime(df['Time'])
    df['Time'] = df['Time'].dt.strftime("%Y-%m-%d %H:%M:%S")

    plt.title("Acquisition time")
    plt.xlabel("Time")
    plt.ylabel("Cumulative exposures")
    ax0.xaxis.set_major_locator(plt.MaxNLocator(10))
    plt.xticks( rotation=25 )
    #ax0.set_xticklabels(labels=df['Time'], rotation=25, rotation_mode="anchor", ha="right", fontsize=6)

    plt.tight_layout()
    fig1 = ax0.get_figure()

    #plt.show()
    if not args.report:
        print('No report, not saving plot')
    else:
        fig1.savefig(main.plot_dir+'/'+xml_session.sessionName+'_exposures.png', dpi=300)
    plt.figure(1).clear()
    plt.close(1)
    plt.cla()
    plt.clf()

    print('Plotted micrograph acquisition times')

def get_mic_date(micpath):
    # This will fetch the micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    # This will get the micrograph time stamp as written in the data acqusition xml file
    data = data["microscopeData"]["acquisition"]["acquisitionDateTime"]
    # Format date
    micDate = dateutil.parser.parse(data)
    micDate = micDate.strftime("%Y-%m-%d %H:%M:%S")

    return micDate

def get_mic_date_xmltree(micpath):
    ## Find data using direct xml query
    tree = ET.parse(micpath)
    root = tree.getroot()

    # FEI EPU xml stuff
    ns = {'p': 'http://schemas.datacontract.org/2004/07/Applications.Epu.Persistence'}
    ns['system'] = 'http://schemas.datacontract.org/2004/07/System'
    ns['so'] = 'http://schemas.datacontract.org/2004/07/Fei.SharedObjects'
    ns['g'] = 'http://schemas.datacontract.org/2004/07/System.Collections.Generic'
    ns['s'] = 'http://schemas.datacontract.org/2004/07/Fei.Applications.Common.Services'
    ns['a'] = 'http://schemas.datacontract.org/2004/07/System.Drawing'
    ns['t'] = 'http://schemas.datacontract.org/2004/07/Fei.Types'

    # Timestamp
    timestamp = root.find('so:microscopeData/so:acquisition/so:acquisitionDateTime', ns).text

    micDate = dateutil.parser.parse(timestamp)
    micDate = micDate.strftime("%Y-%m-%d %H:%M:%S")
    fileTimeString = micDate
    #fileTime = datetime.datetime.strptime(fileTimeString, '%Y-%m-%d %H:%M:%S')

    return fileTimeString

def get_img_sample():
    print(search_mics.randomName0Path)

def count_mics(filelist):
    micNumber = len(filelist)
    print('Found: '+str(micNumber))
    print('')
    print('For the detector: '+str(getScopeNameMag.camera)+' ('+getScopeNameMag.cameraW+' x '+getScopeNameMag.cameraH+'):')
    print('Detector area per shot (nm^2): '+str(getScopeNameMag.area))
    print('Total imaged area (nm^2): '+str(float(getScopeNameMag.area)*float(micNumber)))
    print('')

    if args.timemethod == 'F':
        # Get data acquisition time by file stamp - fast but system timestamp dependent
        # mtime is last modification of file contents - which should be the write time
        # TFS do not write out EXIF data to jpg so can't get creation date there
        # with mtime this comes out with format %Y-%m-%d %H:%M:%S
        # ctime is last change of file inode, permissions changed, file renamed etc
        # with ctime this comes out with format %Y-%m-%d %H:%M:%S.%f - formatFileDate will reformat
        # Most appropriate to use mtime to hope to get the file creation date
        timeFirst = datetime.datetime.fromtimestamp(os.path.getmtime(search_mics.firstXml))
        timeLast = datetime.datetime.fromtimestamp(os.path.getmtime(search_mics.lastXml))
        # Reporting
        #count_mics.first = formatFileDate(str(timeFirst))
        #count_mics.last = formatFileDate(str(timeLast))
        count_mics.first = timeFirst
        count_mics.last = timeLast

    elif args.timemethod == 'X':
        # Get data acquisition time by time written in xml file - more likely accurate
        timeFirst = datetime.datetime.strptime(get_mic_date(search_mics.firstXml),"%Y-%m-%d %H:%M:%S")
        timeLast = datetime.datetime.strptime(get_mic_date(search_mics.lastXml),"%Y-%m-%d %H:%M:%S")
        # Reporting
        count_mics.first = timeFirst
        count_mics.last = timeLast

    # Reporting
    count_mics.number = micNumber
    print('')

# DEV DEV Why do we have two functions called this, sort this out
def count_movies(filelist, dir):
    realPath = dir
    print('Counting movies in raw* directory')
    searchedFiles = filelist
    #searchedFiles = glob.glob(dir+"/**/Data/*xml", recursive=True)

    oneFile = searchedFiles[0]
    name = str(os.path.splitext(os.path.basename(oneFile))[0])
    result = findRaw(name, dir)
    raw = result[0]
    rawExt = result[1]
    # Count movies in appropriate raw folder
    rawNo = countMovies(realPath, raw)

    count_movies.number = rawNo
    count_movies.format = rawExt

    print(str(findRaw.dir)+' size (TB): '+ str(findRaw.size))

def countMovies(path, raw):
    d = path+'/'+raw
    rawMovies = glob.glob(d+'/**/FoilHole_*_Data_*', recursive = True)
    rawNo = len(rawMovies)
    rawName = os.path.basename(os.path.normpath(d))
    #print('Counting movies in '+sessionName+': '+rawNo)
    return rawNo

def findRaw(name, path):
    #print(name)
    #print(path)
    # Find movie and extension
    print('Searching for raw movie directory...')
    # https://stackoverflow.com/questions/1724693/find-a-file-in-python
    try:
        #movie = glob.glob(path+'/raw*/GridSquare*/Data/'+name+'*', recursive = True)[0]
        movie = glob.glob(path+'/raw*/**/Data/'+name+'*', recursive = True)[0]
    except:
        rawDir = 'None'
        findRaw.dir = 'None'
        findRaw.size = 0
        return rawDir
    else:
        extension = str(os.path.splitext(os.path.basename(movie))[1])
        #print(movie)
        #print(extension)
        # Find what raw directory this is in
        full_path = str(movie)
        parent_path = str(path)
        relative_path = os.path.relpath(full_path, parent_path)

        rawDir = Path(relative_path).parts[0]
        findRaw.dir = path+'/'+rawDir

        # raw size
        if args.sizecalc == 'Y':
            print('Calculating raw directory size (TB)')
            size = get_dir_size(findRaw.dir, 'TB')
            findRaw.size = size
            print(str(size)+' TB')
            print()
        else:
            findRaw.size = np.nan

        return rawDir, extension

def calcTimings():
    # The screening time won't make sense for multi grid collections
    # Calculate screening time
    c = xml_session.sessionDate-screeningSession_xml.date
    #print('Screening time: ', c)
    minutes = c.total_seconds() / 60
    hours = minutes / 60

    # Reporting
    calcTimings.screenTimeMin = minutes
    calcTimings.screenTimeHr = hours
    # The screening time won't make sense for multi grid collections

    # Calculate session setup time
    c = count_mics.first-xml_session.sessionDate
    print('Setup time: ', c)
    minutes = c.total_seconds() / 60
    hours = minutes / 60

    # Reporting
    calcTimings.setupTimeMin = minutes
    calcTimings.setupTimeHr = hours

    # Calculate micrograph collection times
    # https://www.geeksforgeeks.org/python-difference-between-two-dates-in-minutes-using-datetime-timedelta-method/
    c = count_mics.last-count_mics.first
    print('Collection time: ', c)
    minutes = c.total_seconds() / 60
    hours = minutes / 60
    print('Delta first to last micrograph (minutes): '+str(round(minutes)))
    print('')
    print('Average data rate (micrographs per hour):')
    collectionRate = round(count_mics.number/(minutes/60))
    print(collectionRate)
    print('')
    print('Average data rate (um^2 per hour):')
    collectionRateArea = str(round(float(collectionRate)*float(getScopeNameMag.area)))
    print(collectionRateArea)
    print('')

    # Reporting
    calcTimings.timeMin = minutes
    calcTimings.timeHr = hours
    calcTimings.rate = collectionRate
    calcTimings.rateArea = collectionRateArea

def plotTimings(f):
    #print('plot timings')
    df = pd.read_csv(f, usecols=['Session date', 'EPU session date', 'Setup time (hrs)', 'Collection time (hrs)', 'Total_EPU_mics'], parse_dates=['EPU session date','Session date'])
    fig0 = plt.figure(0)
    # https://medium.com/@jb.ranchana/easy-way-to-create-stacked-bar-graphs-from-dataframe-19cc97c86fe3

    figHeight = 4
    figWidth = 7

    ax = df.plot.bar(y=["Setup time (hrs)","Collection time (hrs)"], color=['lightblue','orange'], title='Session Timings', figsize=(figWidth,figHeight))

    plt.title("Session timings")
    plt.xlabel("Operation")
    plt.ylabel("Hours")
    #ax.set_xticklabels(labels=df['EPU session date'], rotation=70, rotation_mode="anchor", ha="right")

    plt.tight_layout()
    fig0 = ax.get_figure()
    fig0.savefig(main.plot_dir+'/'+xml_session.sessionName+'_session_timings.png', dpi=300)
    plt.figure(0).clear()
    plt.close()
    plt.cla()
    plt.clf()
    #plt.show()

def formatFileDate(d):
    # https://stackoverflow.com/questions/17594298/date-time-formats-in-python
    # https://www.w3schools.com/python/python_datetime.asp
    # Read in EPU formatted date and time
    new_date = datetime.datetime.strptime(d,"%Y-%m-%d %H:%M:%S.%f")
    # Reformat date into YY-MM-DD HH-MM-SS
    return new_date.strftime("%Y-%m-%d %H:%M:%S")

def getScopeName(micpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    # Diagnostic print data xml file
    #data = json.loads(json.dumps(data))
    #pprint(data)

    # Using FEI xml - Look up instrument, magnification and cross reference against mag table to get real pixel size
    temModel = data["microscopeData"]["instrument"]["InstrumentModel"]
    #print('TEM instrument: '+str(temModel))
    temMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
    #print('TEM magnification: '+str(temMag))

    # Using eBIC dictionary (above) - Look up instrument name and apix
    temName = magTable[temModel]['#name']
    temDlsName = magTable[temModel]['mXX']
    print('eBIC instrument: '+str(temName))
    temApix = magTable[temModel][temMag]
    print('Magnified pixel size: '+str(temApix))
    print('')
    camera = magTable[temModel]["camera"]
    cameraOptimalDose = magTable[temModel]["optimalDose"]
    cameraWpx = magTable[temModel]["cameraWpx"]
    cameraHpx = magTable[temModel]["cameraHpx"]
    cameraWnano = float(temApix) * float(cameraWpx) / 10
    cameraHnano = float(temApix) * float(cameraHpx) / 10
    cameraWmicron = float(temApix) * float(cameraWpx) / 10000
    cameraHmicron = float(temApix) * float(cameraHpx) / 10000
    cameraAreaNano = cameraWnano * cameraHnano
    cameraAreaMicron = cameraWmicron * cameraHmicron
    cameraRatio = float(cameraWpx)/float(cameraHpx)

    # DEV need some error checking in here if it can't find the xml, or the data in the xml

    # Report
    getScopeNameMag.scope = temName
    getScopeNameMag.dlsName = temDlsName
    getScopeNameMag.mag = temMag
    getScopeNameMag.apix = temApix
    getScopeNameMag.area = cameraAreaMicron
    getScopeNameMag.camera = camera
    getScopeNameMag.optimalDose = cameraOptimalDose
    getScopeNameMag.cameraW = cameraWpx
    getScopeNameMag.cameraH = cameraHpx
    getScopeNameMag.cameraRatio = cameraRatio

def getScopeNameMag(micpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    # Diagnostic print data xml file
    #data = json.loads(json.dumps(data))
    #pprint(data)

    # Using FEI xml - Look up instrument, magnification and cross reference against mag table to get real pixel size
    temModel = data["microscopeData"]["instrument"]["InstrumentModel"]
    #print('TEM instrument: '+str(temModel))
    temMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
    #print('TEM magnification: '+str(temMag))

    # Using eBIC dictionary (above) - Look up instrument name and apix
    temName = magTable[temModel]['#name']
    temDlsName = magTable[temModel]['mXX']
    print('eBIC instrument: '+str(temName))
    temApix = magTable[temModel][temMag]
    print('Magnified pixel size: '+str(temApix))
    print('')
    camera = magTable[temModel]["camera"]
    cameraOptimalDose = magTable[temModel]["optimalDose"]
    cameraWpx = magTable[temModel]["cameraWpx"]
    cameraHpx = magTable[temModel]["cameraHpx"]
    cameraWnano = float(temApix) * float(cameraWpx) / 10
    cameraHnano = float(temApix) * float(cameraHpx) / 10
    cameraWmicron = float(temApix) * float(cameraWpx) / 10000
    cameraHmicron = float(temApix) * float(cameraHpx) / 10000
    cameraAreaNano = cameraWnano * cameraHnano
    cameraAreaMicron = cameraWmicron * cameraHmicron
    cameraRatio = float(cameraWpx)/float(cameraHpx)

    # DEV need some error checking in here if it can't find the xml, or the data in the xml

    # Report
    getScopeNameMag.scope = temName
    getScopeNameMag.dlsName = temDlsName
    getScopeNameMag.mag = temMag
    getScopeNameMag.apix = temApix
    getScopeNameMag.area = cameraAreaMicron
    getScopeNameMag.camera = camera
    getScopeNameMag.optimalDose = cameraOptimalDose
    getScopeNameMag.cameraW = cameraWpx
    getScopeNameMag.cameraH = cameraHpx
    getScopeNameMag.cameraRatio = cameraRatio

def getFractionNo(micpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    # The FractionationSettings are not always in the same list position in a:KeyValueOfstringanyType
    # DEV DEV The SuperResolutionFactor is also in this list
    keyValueList = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"]

    # Loop through the list to find the FractionationSettings list position
    i = 0
    for value in keyValueList:
       key = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][i]["a:Key"]
       if key == "FractionationSettings":
          FractionationSettings = i
       i = i+1

    # Data structure depends on camera type K3 versus F4
    # Work based on camera type in getScopeNameMag.camera
    if getScopeNameMag.camera == 'K3':

       try:
          fractionNo = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][FractionationSettings]["a:Value"]["b:NumberOffractions"]
       except:
          fractionNo = 'Unknown'

       fractionationType = 'K3-fractionation'

    elif getScopeNameMag.camera == 'F4' or getScopeNameMag.camera == 'F4i':

       # EER has no fraction dose, as is fraction-less, otherwise F4 stores fractions in list, count to get number of fractions
       fractionationType = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][FractionationSettings]["a:Value"]["@i:type"]

       if fractionationType == 'b:VaridistantFractionation':
          try:
             fractionList = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][FractionationSettings]["a:Value"]["b:DoseFractions"]["a:anyType"]
             fractionNo = len(fractionList)
          except:
             fractionNo = 'Unknown'
       elif fractionationType == 'b:EerFractionation':
          fractionNo = 'Unknown'

    print()
    print('Fraction number determined as: '+str(fractionNo)+' (in fractionation mode: '+fractionationType+')')

    return fractionNo

def getCameraDim(micpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    print('Sensor dimensions at DataAcquisition:')

    # Find the reported camera readout size and Binning
    cameraXbin = data["microscopeData"]["acquisition"]["camera"]["Binning"]["a:x"]
    cameraYbin = data["microscopeData"]["acquisition"]["camera"]["Binning"]["a:y"]
    cameraXpx = data["microscopeData"]["acquisition"]["camera"]["ReadoutArea"]["a:width"]
    cameraYpx = data["microscopeData"]["acquisition"]["camera"]["ReadoutArea"]["a:height"]

    cameraXpxscaled = float(cameraXpx) * float(cameraXbin)
    cameraYpxscaled = float(cameraYpx) * float(cameraYbin)
    print('Sensor x (px): '+str(cameraXpxscaled))
    print('Sensor y (px): '+str(cameraYpxscaled))

    # Repeat code to get pixel size at data acquisition
    temModel = data["microscopeData"]["instrument"]["InstrumentModel"]
    temMag = data["microscopeData"]["optics"]["TemMagnification"]["NominalMagnification"]
    temApix = magTable[temModel][temMag]

    # Calculate sensor dimensions in Angstroms
    cameraXang = float(temApix) * cameraXpxscaled
    cameraYang = float(temApix) * cameraYpxscaled
    cameraXmicron = float(temApix) * cameraXpxscaled / 1e4
    cameraYmicron = float(temApix) * cameraYpxscaled / 1e4
    print('Sensor x (um): '+str(cameraXmicron))
    print('Sensor y (um): '+str(cameraYmicron))

    getCameraDim.cameraXpx = cameraXpxscaled
    getCameraDim.cameraYpx = cameraYpxscaled

    return [cameraXmicron, cameraYmicron]

def getStageTilt(micpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    # Find the stage Alpha (DEV DEV think the units might be 1/100th)
    stageAlpha = data["microscopeData"]["stage"]["Position"]["A"]
    stageBeta = data["microscopeData"]["stage"]["Position"]["B"]

    return [stageAlpha, stageBeta]

def roundup(n, decimals=0):
    # https://realpython.com/python-rounding/
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

def screeningSession_xml(path):
    # Get
    screeningDataModel = (path+'/ScreeningSession.dm')
    # This will fetch the atlas xml data
    with open(screeningDataModel, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["ScreeningSessionXml"]

    try:
        screenTime = data["StartDateTime"]
    except:
        screenTime = 'Unknown'

    #print('Screening session start time: '+str(formatEPUDate(screenTime)))

    try:
        sessionName = data["Name"]['#text']
    except:
        sessionName = 'Unknown'

    # Count the number of samples, assumed by number of atlases collected
    searchedFiles = glob.glob(path+"/Sample*/Atlas/Atlas_*.jpg")

    #print('Atlas screening session name: '+str(sessionName))

    # Report
    screeningSession_xml.date = formatEPUDate(screenTime)
    screeningSession_xml.name = str(sessionName)
    screeningSession_xml.atlasNo = len(searchedFiles)

def atlas_xml(atlas_dir: Path) -> Dict[str, Any]:
    # Find atlas first
    find_atlas(atlas_dir)
    # This will fetch the atlas xml data
    with open(main.atlasxml, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["SampleXml"]

    # Get autoloader slot name
    try:
        autoSlotName = data["Name"]["#text"]
    except:
        autoSlotName = 'None' # Probably nothing entered into TUI

    #print('Autoloader position name: '+str(autoSlotName))

    #Report
    atlas_xml.atlasName = autoSlotName

# This function is general use, provide the atlas path up to Sample number
# find image, autoslotname and atlas xml path
# Could use this function for all atlas operations?
# DEV DEV not finished, finish it and then use in PDF reporting
def atlas_image(sample_path: Path) -> Dict[str, Any]:
    atlas_xml = sample_path+'/Sample.dm'
    # Glob can return empty list as not all slots will have an atlas
    try:
        atlas_img = glob.glob(sample_path+'/Atlas/Atlas_*.jpg')[-1]
    except:
        atlas_img = ''

    try:
        with open(atlas_img) as f:
            # This will fetch the atlas xml data
            with open(atlas_xml, "r") as xml:
                for_parsing = xml.read()
                data = xmltodict.parse(for_parsing)
            data = data["SampleXml"]

            # Get autoloader slot name
            try:
                autoSlotName = data["Name"]["#text"]
            except:
                autoSlotName = 'None' # Probably nothing entered into TUI

            atlas_autoslot_name = autoSlotName

    except IOError:
            #print("No atlas found in "+sample_path+" using atlas_image function")
            atlas_img = 'None'
            atlas_xml = 'None'
            atlas_autoslot_name = 'None'

    # Return as list for use elsewhere
    return [atlas_img, atlas_xml, atlas_autoslot_name]

def find_atlas(atlas_dir: Path) -> Dict[str, Any]:
    with open(main.xml, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    # Autoloader slot (starts at 0)
    autoSlot = data["AutoloaderSlot"]
    autoSlotReal = float(autoSlot)+1

    # So the SampleID is
    sampleNo = str('Sample'+str(round(autoSlotReal)))

    # Look for EpuSession.dm in defined top level EPU session directory
    path_to_file= str(main.atlas+'/'+sampleNo+'/Sample.dm')
    path = Path(path_to_file)

    # Look for atlas
    if path.is_file():
        #print(f'The file {path_to_file} exists')
        # This is the path to the epusession xml, if present, populate atlas variables
        main.atlasxml = path_to_file
        main.atlasxmlPath = path

        main.sampledir = str(main.atlas+'/'+sampleNo)
        main.atlasdir = str(main.sampledir+'/Atlas')

        main.atlasimg = str(glob.glob(main.sampledir+'/Atlas/Atlas_*.jpg')[-1]) # glob.glob returns full path
        main.atlasimgPath = str(main.atlasimg)
    else:
        print(f'The file {path_to_file} does not exist, exiting with error')
        exit(1)

def print_atlas_xml(atlas_dir: Path) -> Dict[str, Any]:
    # Find atlas first
    find_atlas(atlas_dir)
    # This will print the atlas xml data
    with open(main.atlasxml, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["SampleXml"]

    data = json.loads(json.dumps(data))
    pprint(data)

def print_epu_xml(xml_path: Path) -> Dict[str, Any]:
    # Use this function for troubleshooting/viewing the raw xml to find data structure
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    data = json.loads(json.dumps(data))
    pprint(data)

def xml_session(xml_path: Path) -> Dict[str, Any]:
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    # Location of EPU session directory on which this script was ran
    realPath = os.path.realpath(xml_path)

    # EPU version
    epuId = data["Version"]["@z:Id"]
    epuBuild = data["Version"]["a:_Build"]
    epuMajor = data["Version"]["a:_Major"]
    epuMinor = data["Version"]["a:_Minor"]
    epuRevision = data["Version"]["a:_Revision"]
    epuVersion = str(epuMajor)+'.'+str(epuMinor)+'.'+str(epuRevision)+'-'+str(epuId)+'.'+str(epuBuild)

    # Output format
    doseFractionOutputFormat = data["DoseFractionsOutputFormat"]["#text"]

    # Autoloader slot (starts at 0)
    autoSlot = data["AutoloaderSlot"]
    autoSlotReal = float(autoSlot)+1
    #print('Autoloader position: '+str(autoSlotReal))

    # SampleXml is a list denoed inside [], pull 0 out
    # Then the text of AtlasId is in a key pair in the dictionary denoted by {}, so call that key pair
    # Use the print_epu_xml def to see full [] and {} formatted data structure
    atlasDir = data["Samples"]["_items"]["SampleXml"][0]["AtlasId"]["#text"]
    #print('Atlas location: '+atlasDir)
    #print('')

    # Session name and creation time
    sessionName = xml_sessionName(xml_path)
    #sessionName = data["Name"]["#text"]
    #sessionDate = data["StartDateTime"]
    sessionDate = data["Samples"]["_items"]["SampleXml"][0]["StartDateTime"]
    sessionDateFormat = formatEPUDate(sessionDate)
    #print('EPU session name: '+sessionName)
    #print('EPU session created: '+str(sessionDateFormat))
    #print('EPU session created: '+str(sessionDate))
    #print('')

    # Grid type - lacey or holeycarbon or holeygold
    try:
       gridType = data["Samples"]["_items"]["SampleXml"][0]["GridType"]
    except:
       gridType = 'Unknown'

    # The I0 filter settings may hint at what grid type is being used
    try:
       I0set = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["IsCalibrated"]
    except:
       I0set = 'Unknown'

    try:
       I0MaxInt = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["MaximumIntensity"]
    except:
       I0MaxInt = 'Unknown'

    try:
       I0MinInt = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["MinimumIntensity"]
    except:
       I0MinInt = 'Unknown'

    # Clustering method
    clustering = data["ClusteringMode"]
    #print('Clustering mode: '+clustering)
    clusteringRadius = data["ClusteringRadius"]
    clusteringRadiusMicron = float(clusteringRadius)*1e6
    #print('Clustering radius: '+str(clusteringRadiusMicron)+' microns')
    #print('')
    if clustering == 'ClusteringWithImageBeamShift':
        afisMode = 'AFIS'
        afisRadius = str(clusteringRadiusMicron)
    else:
        afisMode = 'Accrt'
        afisRadius = np.nan

    focusWith = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["AutoFocusArea"]["FocusWith"]["#text"]
    focusRecurrence = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["AutoFocusArea"]["Recurrence"]["#text"]
    #print('Focus with: '+focusWith)
    #print('Focus recurrence: '+focusRecurrence)
    #print('')

    delayAfterImageShift = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DelayAfterImageShift"]
    delayAfterStageShift = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DelayAfterStageShift"]
    #print('Delay after Image Shift: '+delayAfterImageShift+' seconds')
    #print('Delay after Stage Shift: '+delayAfterStageShift+' seconds')
    #print('')

    # Send xml dict over to function to get defocus range
    defocusRange=getDefocusRange(data)

    # In some cases the defocus list is a single value, test to deal with
    # Also find max and min defocus values
    if isinstance(defocusRange, str):
        if defocusRange == "xml read error":
           defocusRangeMicron = "xml read error"
           defocusRangeRound = "xml read error"
           defocusMax = "xml read error"
           defocusMin = "xml read error"
        else:
           #defocusRangeMicron = float(defocusRange) * 1e6
           defocusRangeMicron = float(defocusRange) * 1e6
           defocusRangeRound = roundup(defocusRangeMicron, 2)
           defocusMax = min(defocusRangeRound)
           defocusMin = max(defocusRangeRound)
    elif isinstance(defocusRange, list):
        if defocusRange[0] == "xml read error":
           defocusRangeMicron = "xml read error"
           defocusRangeRound = "xml read error"
           defocusMax = "xml read error"
           defocusMin = "xml read error"
        else:
           #defocusRangeMicron = [float(element) * 1e6 for element in defocusRange]
           defocusRangeMicron = [float(element) * 1 for element in defocusRange]
           defocusRangeRound = [roundup(num, 2) for num in defocusRangeMicron]
           defocusMax = min(defocusRangeRound)
           defocusMin = max(defocusRangeRound)
    #print('Defocus range (um): ')
    #print(defocusRangeRound)
    #print('Defocus max: '+str(defocusMax))
    #print('Defocus min: '+str(defocusMin))

    # Report
    xml_session.epuVersion = epuVersion
    xml_session.doseFractionOutputFormat = doseFractionOutputFormat
    xml_session.sessionName = sessionName
    xml_session.sessionDate = sessionDateFormat
    #xml_session.sessionDate = sessionDate
    xml_session.clustering = clustering
    xml_session.clusteringRadius = clusteringRadiusMicron
    xml_session.focusWith = focusWith
    xml_session.focusRecurrence = focusRecurrence
    xml_session.xmlPath = realPath
    # This should get the directory name two levels up from EpuSession.dm and on diamond systems be the DLS Visit
    xml_session.autoSlot = autoSlotReal
    xml_session.atlasDirOrig = atlasDir

    xml_session.defocus = defocusRangeRound
    xml_session.defocusMax = defocusMax
    xml_session.defocusMin = defocusMin

    xml_session.delayImageShift = delayAfterImageShift
    xml_session.delayStageShift = delayAfterStageShift
    xml_session.afisMode = afisMode
    xml_session.afisRadius = afisRadius
    xml_session.gridType = gridType
    xml_session.I0set = I0set
    xml_session.I0MaxInt = round(float(I0MaxInt))
    xml_session.I0MinInt = round(float(I0MinInt))

    print('Finished gathering metadata from main EPU session file')
    print('')

def xml_sessionName(xml_path):
    # It is necessary to have a function for getting xml session name elsewhere in script
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    sessionName = data["Name"]["#text"]

    return sessionName

def getDefocusRange(data):
    # Need to deal with cases where it's multi shot or single shot
    templates = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]

    if isinstance(templates, list):
        shotType = 'Multishot'
        # DEV DEV This defocus code might need a revisit if you get lots of errors
        d = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][0]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]
        if d.get("a:double"):
           try:
              df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][0]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]["a:double"]
              #Sometimes the values contain unicode en-dash and not ASCII hyphen
              #df.replace('\U00002013', '-')
           except:
              print('Warning, could not find defocus range in xml file')
              df = ['xml read error']
        else:
            try:
               df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]["a:_items"]["a:double"]
               #Sometimes the values contain unicode en-dash and not ASCII hyphen
            #df.replace('\U00002013', '-')
            except:
               print('Warning, could not find defocus range in xml file')
               df = ['xml read error']
    else:
        shotType = 'Single'
        # There is sometimes a data structure change I think in single shot acqusition, cause currently unknown, check for it
        d = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]
        if d.get("a:double"):
            try:
               df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]["a:double"]
               #Sometimes the values contain unicode en-dash and not ASCII hyphen
               #df.replace('\U00002013', '-')
            except:
               print('Warning, could not find defocus range in xml file')
               df = ['xml read error']
        else:
            try:
               df = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ImageAcquisitionSettingXml"]["Defocus"]["a:_items"]["a:double"]
               #Sometimes the values contain unicode en-dash and not ASCII hyphen
               #df.replace('\U00002013', '-')
            except:
               print('Warning, could not find defocus range in xml file')
               df = ['xml read error']

    getDefocusRange.shotType = shotType

    # Remember df in this case means defocus, not dataframe!!
    # Sometimes there is a single value in the defocus list and then this gets stored as a single string
    if isinstance(df, str):
       # Convert to list is str and thus single value
       df = df.split(sep=None, maxsplit=-1)

    # Check for error, which is stored as single item list
    read = df[0]
    if df[0] == "xml read error":
       return df
    # Otherwise convert metres into microns
    else:
       dfMicron = [float(item) * 1e6 for item in df]
       return dfMicron

def formatEPUDate(d):
    # Returns formatted in datetime, needs to be string for printing

    # https://stackoverflow.com/questions/17594298/date-time-formats-in-python
    # https://www.w3schools.com/python/python_datetime.asp
    # https://www.tutorialexample.com/python-detect-datetime-string-format-and-convert-to-different-string-format-python-datetime-tutorial/amp/
    # Read in EPU formatted date and time - remember input is a string
    epuDate = dateutil.parser.parse(d)
    epuDate = epuDate.strftime("%Y-%m-%d %H:%M:%S")
    #new_date = datetime.datetime.strptime(d,"%Y-%m-%dT%H:%M:%S.%fZ")
    # Reformat date into YY-MM-DD HH-MM-SS
    #return new_date.strftime("%Y-%m-%d %H:%M:%S")
    return datetime.datetime.strptime(epuDate,"%Y-%m-%d %H:%M:%S")

def xml_templateShotsPerHole(xml_path: Path) -> Dict[str, Any]:
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    templates = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]
    if isinstance(templates, list):
        print('Multishot detected')
        templates = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]
        length=len(templates)
        shotsPerHole = length
    else:
        print('Single shot detected...')
        length=1
        shotsPerHole = 1

    return shotsPerHole

def xml_targets(xml_path: Path) -> Dict[str, Any]:
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    # get epu directory path
    epudir = os.path.dirname(xml_path)

    ## Look through the square metadata to find all squares
    # You can use this to see which foilholes were selected, which have been corrected
    print('Analysing GridSquares, FoilHoles & DataAcquisitions selected by user')
    squareData = glob.glob(epudir+'/Metadata/*.dm', recursive = True)

    # Get shots per hole
    shotsPerHole = xml_templateShotsPerHole(xml_path)

    # Empty dataframe for collecting square meta data
    sqMetadataDf = df = pd.DataFrame({"SquareName":[],
                             "DataAqNo":[]})
    # Loop through each square metadata file
    for dm in tqdm(squareData):
       # This function will return the square identifier and number of exposures on it
       metadata = squareMetadata(dm)
       FoilAqNumber = metadata[1]
       squareName = metadata[0]
       df2 = pd.DataFrame({"SquareName":[squareName],
                           "FoilAqNo":[FoilAqNumber]})
       sqMetadataDf = pd.concat([sqMetadataDf, df2], ignore_index=True, sort=False)
    # Calculate how many acquitions per square
    sqMetadataDf['DataAqNo'] = sqMetadataDf['FoilAqNo'] * shotsPerHole

    #print(sqMetadataDf.dropna())
    #print(sqMetadataDf['DataAqNo'].dropna().sum())
    totalFoilHoleNo = sqMetadataDf['FoilAqNo'].dropna().sum()

    # Do some statistics on target selections
    # DEV DEV DEV errors in bi28713-27, there is a different key pair path and then a list inside, see glom use attempt in squareMetadata function
    # This needs fixing
    # The try statements below will at least handle the exception so analysis doesn't bail out
    # Also finding some sessions have very low foilhole numbers per square, is this right? see bi23268-92

    import statistics

    sqDataAqNo = sqMetadataDf['DataAqNo'].dropna().tolist()
    try:
       avSqDataAqNo = round(statistics.mean(sqDataAqNo))
    except:
       avSqDataAqNo = np.nan
    try:
       stdevSqDataAqNo = round(statistics.stdev(sqDataAqNo))
    except:
       stdevSqDataAqNo = np.nan

    sqFoilHoleNo = sqMetadataDf['FoilAqNo'].dropna().tolist()
    try:
       avSqFoilHoleNo = round(statistics.mean(sqFoilHoleNo))
    except:
       avSqFoilHoleNo = np.nan
    try:
       stdevSqFoilHoleNo = round(statistics.stdev(sqFoilHoleNo))
    except:
       stdevSqFoilHoleNo = np.nan

    print()
    print('Average Data Acqutisitions targetted per Square: '+str(avSqDataAqNo)+' (+/- '+str(stdevSqDataAqNo)+')')
    print('Average Foil Holes selected per Square: '+str(avSqFoilHoleNo)+' (+/- '+str(stdevSqFoilHoleNo)+')')
    print('Shots per Foil Hole: '+str(shotsPerHole))
    print()

    #input("Press Enter to continue...")

    #squareImages = epudir+'/Images-Disc1'
    squareImages = main.epuSessionPath+'/Images-Disc1'
    print('Square Image path: '+squareImages)
    print()

    # Get count of all squares on grid number from xml
    squares = data["Samples"]["_items"]["SampleXml"][0]["GridSquares"]["a:m_serializationArray"]["b:KeyValuePairOfintSerializedReferenceOfGridSquareXml_PsqDC3X0m6koiAE_P"]
    totalSquares=len(squares)
    # Get count of targetted squares on grid number square images
    targetSquares=len([name for name in os.listdir(squareImages) if os.path.isdir(os.path.join(squareImages, name))])
    # Get count of number of squares that were actually collected on
    squares = glob.glob(epudir+'/Images-Disc1/GridSquare*', recursive = True)
    collectedSquares = len(squares)

    # Get hole size, spacing and rotation
    holeSize = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["HoleSize"]
    holeSizeMicron = float(holeSize)*1e6
    #print('Hole size: '+str(holeSizeMicron)+' microns (measured in EPU)')

    holeSpace = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["HoleSpacing"]
    holeSpaceMicron = float(holeSpace)*1e6/2
    #print('Hole spacing: '+str(holeSpaceMicron)+' microns (measured in EPU)')
    #print('')

    holeRotation = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["Rotation"]

    # Match measured hole size to regularised hole size
    holeSizeRegular, holeSpaceRegular = holeRegularisation(holeSizeMicron, holeSpaceMicron)

    # Report
    xml_targets.totalSquares = totalSquares
    xml_targets.targetSquares = targetSquares
    xml_targets.collectedSquares = collectedSquares
    # Note these are only squares that actually have data acquisition images collected on them
    xml_targets.listSquares = squares

    xml_targets.totalFoilHoles = totalFoilHoleNo
    xml_targets.avSqFoilHoleNo = avSqFoilHoleNo
    xml_targets.stdevSqFoilHoleNo = stdevSqFoilHoleNo

    xml_targets.holeSize = roundup(holeSizeMicron,1)
    xml_targets.holeSpace = roundup(holeSpaceMicron,1)

    xml_targets.holeSizeRegular = holeSizeRegular # This comes as a list, you need to sort it out so it's a single string
    xml_targets.holeSpaceRegular = holeSpaceRegular

    xml_targets.holeRotation = holeRotation

    xml_targets.dataAqList = sqMetadataDf

def holeRegularisation(size, space):
    print('Regularising measured hole size and spacing')
    # https://www.quantifoil.com/products/quantifoil/quantifoil-circular-holes
    # R N/N listed as of 10/02/2023
    # 0.6/1, 1/1, 1/2, 1/4, 1.2/1.3, 1.2/20, 2/1, 2/2, 2/4, 3/3, 3/5, 3.5/1, 5/2, 5/20, 6/6.5, 8/3, 10/5, 10/10, 10/20, 17/5

    #create DataFrame
    # Currently seperate data frames but should you use them in a single dataframe?
    #dfSize = pd.DataFrame({'size': [0.6,1.0, 1.2, 2.0,3.0,3.5,5.0,6.0,8.0,10.0,17.0]})
    #dfSpace = pd.DataFrame({'space': [1.0, 1.3, 2.0, 3.0, 4.0, 5.0, 6.5, 10, 20]})

    # Difficult to tell between 1 and 1.2 micron hole so remove 1 micron from regularisation list
    dfSize = pd.DataFrame({'size': [0.6,1.2,2.0,3.0,3.5,5.0,6.0,8.0,10.0,17.0]})
    dfSpace = pd.DataFrame({'space': [1.0,1.3,2.0,3.0,4.0,5.0,6.5,10,20]})

    dfSize_closest = dfSize.iloc[(dfSize['size']-float(size)).abs().argsort()[:1]]
    size = dfSize_closest['size'].tolist()[0]

    dfSpace_closest = dfSpace.iloc[(dfSpace['space']-float(space)).abs().argsort()[:1]]
    space = dfSpace_closest['space'].tolist()[0]

    # Overrule some spacings spacing
    # If hole size found to be 1.2, then force spacing to be 1.3
    if size == 1.2:
        space = 1.3

    # If hole size is 2 and space incorrectly 1.3, force to 1.0
    if size == 2.0 and space == 1.3:
        space = 1.0

    return str(size), str(space)

    # Now write a function that calculates upper and lower bounds for each hole size
    # Round measured hole size to nearest value within calculated bounds

def gatherScreening(squareNames, report):
    print('Gathering gridsquare data to create GridSquare report...')
    # THIS SECTION IS THE NEXT TO WORK ON, FOR SPT SCREENING PROJECT
    # Now you can reference search_mics.sessionStructure dataframe and pull out Atlas, GridSquare and some randomly selected Micrographs

    # Create FPDF page, insert atlas, grid square and 6 random micrographs
    print('Generating session PDF screening report...')
    # You need to look into another PDF tool, this is getting silly
    # https://stackoverflow.com/questions/47263621/page-breaks-for-images-in-python-fpdf
    colwidth = 40
    col0width = 65
    col1width = 50
    col2width = 26 # 420 pixels
    colH = 4
    # https://towardsdatascience.com/how-to-create-pdf-reports-with-python-the-essential-guide-c08dd3ebf2ee
    # https://stackoverflow.com/questions/51864730/python-what-is-the-process-to-create-pdf-reports-with-charts-from-a-db
    # Initialise new PDF
    pdf = FPDF()

    # Loop through square list
    for path in tqdm(squareNames):
        #print(path)
        # Convert square path to square name
        gridsqpath = os.path.normpath(path)
        #print(gridsqpath)
        pathList = gridsqpath.split(os.sep)
        GridSquare = [match for match in pathList if "GridSquare" in match]
        #print('Working on: '+str(GridSquare[0]))

        # Look up square name in search_mics.sessionStructure dataframe and pull out relevant rows
        #print(search_mics.sessionStructure)
        df = search_mics.sessionStructure[(search_mics.sessionStructure == GridSquare[0]).any(axis=1)]
        #print(df)
        #print(df['GridSquare_path'])
        #print(df['DataAcquisition_name'])

        # Atlas
        atlasimgpath = main.atlasimgPath

        # Gridsquare
        try:
            gridsqpath = df['GridSquare_path'].iloc[0]
        except:
           continue # skip this grid square
        gridsqimg = df['GridSquare_img'].iloc[0]
        gridsqimgpath = str(gridsqpath+'/'+gridsqimg)

        # Count how many micrographs on this grid square
        miccount = len(df['DataAcquisition_path'])
        #print('Mic count on square: '+str(miccount))
        # Select up to 4 random micrographs from data structure dataframe
        if miccount > 4:
            selectsize = 4
        else:
            selectsize = miccount
        #print('Selected '+str(selectsize)+' micrographs')
        micpath = df['DataAcquisition_path'].iloc[0]
        dfRandomMics = df.sample(n = selectsize)['DataAcquisition_name']
        #print(dfRandomMics)

        # New page for template
        pdf.add_page()
        pdf.set_xy(10, 0)
        pdf.set_font('arial','B', 14)
        pdf.cell(190, 5, "", 0, 2)
        pdf.cell(190, 8, "eBIC Screening Report", 0, 2, 'R')
        pdf.set_font('arial','I', 8)
        pdf.cell(190, 6, str(GridSquare[0]), 0, 2, 'R')
        # Establish initial position for text
        pdf.set_xy(0, 22)
        pdf.set_font('arial','', 8)
        pdf.cell(10,0)

        # Establish initial position for text
        pdf.set_xy(0, 22)
        pdf.set_font('arial','', 8)
        pdf.cell(10,0)

        # Insert some random micrographs
        micSizeWidth = 70
        micSizeHeight = micSizeWidth / float(getScopeNameMag.cameraRatio)
        xspacing = micSizeWidth+10
        yspacing = micSizeHeight

        # Need a function to look at ratio of image to know how to insert image, code reuse from session report image insertion
        # Insert images
        startx = pdf.get_x()+5
        starty = pdf.get_y()
        pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
        pdf.cell(-20)
        ##Atlas and GridSquare
        pdf.image(atlasimgpath, x = startx+(xspacing*0), y = starty+(100*0), w = micSizeWidth, h = micSizeWidth, type = '', link = '')
        pdf.image(gridsqimgpath, x = startx+(xspacing*1), y = starty+(100*0), w = micSizeWidth, h = micSizeHeight, type = '', link = '')

        ##Micrographs - may be different numbers of mics, add up to n in a loop
        i = 0
        for mic in dfRandomMics:
            if (i % 2) == 0:
                multiplier = 1
            else:
                multiplier = 0

            #print('Inserting '+str(mic))
            pdf.image(str(micpath)+'/'+mic+'.jpg', x = startx+(xspacing*multiplier), y = starty+(yspacing*1), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
            i +=1
        #pdf.image(str(micpath)+'/'+dfRandomMics.iloc[0]+'.jpg', x = startx+(xspacing*0), y = starty+(yspacing*1), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
        #pdf.image(str(micpath)+'/'+dfRandomMics.iloc[1]+'.jpg', x = startx+(xspacing*1), y = starty+(yspacing*1), w = micSizeWidth, h = micSizeHeight, type = '', link = '')

        ##Motion and CTF analyses
        if analyseProcessed.motionCor == 'Y':
            #plotCorPerMic(df, 'rlnAccumMotionTotal', 'Total motion per mic ('+str(GridSquare[0])+')', '_motion_total_per_mic_'+str(GridSquare[0]))
            # plotPerMic(df, column, bins, fill, edge, title, yaxis, xaxis, filename)
            plotPerMic(df, 'rlnAccumMotionTotal', 'filter', 40, 'lightgreen', 'green', 'Total motion per mic ('+str(GridSquare[0])+')', 'Micrograph count', 'Total motion per micrograph', '_motion_total_per_mic_'+str(GridSquare[0]))
            pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_motion_total_per_mic_'+str(GridSquare[0])+'.png', x = startx+(xspacing*0), y = starty+(yspacing*2), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
        if analyseProcessed.ctfFind == 'Y':
            #plotCtfPerMic(df, 'rlnCtfMaxResolution', 'CTF max per mic ('+str(GridSquare[0])+')', '_ctf_max_per_mic_'+str(GridSquare[0]))
            # plotPerMic(df, column, bins, fill, edge, title, yaxis, xaxis, filename)
            plotPerMic(df, 'rlnCtfMaxResolution', 'filter', 40, 'gold', 'darkorange', 'CTF max per mic ('+str(GridSquare[0])+')', 'Micrograph count', 'CTF max per micrograph', '_ctf_max_per_mic_'+str(GridSquare[0]))
            pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_ctf_max_per_mic_'+str(GridSquare[0])+'.png', x = startx+(xspacing*1), y = starty+(yspacing*2), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
        ##Picking analyses
        if analyseProcessed.extract == 'Y':
            print('extract found')
            # Plot ptcls per mic for this gridsquare
            #plotPerMic(df, 'ptcls_per_mic', 40, 'deeppink', 'darkviolet', 'Particles per mic ('+str(GridSquare[0])+')', 'Micrograph count', 'Particle per micrograph', '_ptcls_per_mic_'+str(GridSquare[0]))
            plotPerMic(df, 'ideal_picking', 'filter', 40, 'violet', 'indigo', 'Particles per mic ('+str(GridSquare[0])+')', 'Micrograph count', 'Particle per micrograph', '_ptcls_per_mic_'+str(GridSquare[0])+'_norm')
            plotPerMic(df, 'pct_clustered', 'nofilter', 40, 'tomato', 'k', 'Particles clustering ('+str(GridSquare[0])+')', 'Micrograph count', 'Percentage of clustered particles', '_ptcls_per_mic_clustering_'+str(GridSquare[0]))
            plotPerMic(df, 'n_clusters_', 'nofilter', 40, 'tomato', 'k', 'Particles cluster number ('+str(GridSquare[0])+')', 'Micrograph count', 'Number of particle clusters', '_ptcls_per_mic_clusters_'+str(GridSquare[0]))
            # Insert images
            #pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_ptcls_per_mic_'+str(GridSquare[0])+'.png', x = startx+(xspacing*0), y = starty+(yspacing*3), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
            pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_ptcls_per_mic_'+str(GridSquare[0])+'_norm.png', x = startx+(xspacing*0), y = starty+(yspacing*3), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
            pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_ptcls_per_mic_clustering_'+str(GridSquare[0])+'.png', x = startx+(xspacing*1), y = starty+(yspacing*3), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
            pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_ptcls_per_mic_clusters_'+str(GridSquare[0])+'.png', x = startx+(xspacing*1), y = starty+(yspacing*4), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
        ##2D analyses

    # When out of loop, save GridSquare screening report
    # Write out report
    pdf.output(report, 'F')

    # Report
    print('Generated PDF screening report in '+report)
    print()

def squareMetadata(xml_path: Path) -> Dict[str, Any]:
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["GridSquareXml"]

    # This will allow you to build nested dictionary paths based on conditional statements
    # DEV DEV DEV Except haven't figured out how to get it to work with nested lists
    # Hence still failing on some sessions - see bi28713-27
    # https://bond-kirill-alexandrovich.medium.com/nested-dictionaries-in-python-5a362e03ec6c
    import glom

    # This function looks through the Metadata dm's to find out which squares had foil holes selected on this during the EPU session
    squareName = data["TargetLocations"]["GridSquareId"]

    # Foil hole list won't exist if no Foil Holes were defined on the square
    try:
       # This returns a list of every exposure planned on a Square
       dataAqList = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"]
       #base = "TargetLocations.TargetLocationsEfficient.a:m_serializationArray.b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"
       #dataAqList = glom.glom(data, base)
       empty = 'N'
    except:
       try:
          # Sometimes it doesn't have "TargetLocationsEfficient" - as in bi28713-27
          dataAqList = data["TargetLocations"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"]
          #base = "TargetLocations.a:m_serializationArray.b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"
          #dataAqList = glom.glom(data, base)
          empty = 'N'
       except:
          dataAqList = []
          empty = 'Y'

    if empty == 'N':
       # Can we loop through the dataacquisitions to create a dataframe of which are on/off, collected/not
       # Loop through the list of data acquisitions

       # Empty dataframe for collecting Foil Hole data
       # Keep in mind that this dataframe contains lots of useful information, including what was collected or not
       # You may want to use it for something another time
       foilHoleMetadataDf = df = pd.DataFrame({"GridSquareId":[],
                                #"TargetLocationId":[],
                                #"DataAcquisitionAreaId":[],
                                #"ImageAcquisitionSettingId":[],
                                "FoilHolePx_x":[],
                                "FoilHolePx_y":[],
                                "Selected":[],
                                "State":[]})

       i = 0
       aqCount = 0
       for value in dataAqList:
          # Retrieve info on each exposure
          GridSquareId = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][i]["b:value"]["GridSquareId"]
          # This is roughly how you would then build a new path to pull out data, list syntax is not working though, see bi28713-27
          #path = base+"["+str(i)+"]."+"b:value.GridSquareId"
          #print(xml_path)
          #print(path)
          #test = glom.glom(data, path)

          path = "b:value.GridSquareId"
          PixelCenter_x = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][i]["b:value"]["PixelCenter"]["c:x"]
          PixelCenter_y = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][i]["b:value"]["PixelCenter"]["c:y"]

          selected = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][i]["b:value"]["Selected"]
          state = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][i]["b:value"]["State"]

          # For the following, not every FoilHole has these defined, come back to this DEV DEV DEV
          #TargetLocationId = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][i]["b:value"]["ParticleAcquisitions"]["KeyValuePairs"]["b:KeyValuePairOfintParticleAcquisitionDataXmlBpEWF4JT"]["b:value"]["c:TargetLocationId"]

          #DataAcquisitionAreaId = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][i]["b:value"]["ParticleAcquisitions"]["KeyValuePairs"]["b:KeyValuePairOfintParticleAcquisitionDataXmlBpEWF4JT"]["b:value"]["c:DataAcquisitionAreaId"]

          #ImageAcquisitionSettingId = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][i]["b:value"]["ParticleAcquisitions"]["KeyValuePairs"]["b:KeyValuePairOfintParticleAcquisitionDataXmlBpEWF4JT"]["b:value"]["c:ImageAcquisitionSettingId"]

          # Create new dataframe with info on acquisitions
          df2 = pd.DataFrame({"GridSquareId":[GridSquareId],
                                #"TargetLocationId":[TargetLocationId],
                                #"DataAcquisitionAreaId":[DataAcquisitionAreaId],
                                #"ImageAcquisitionSettingId":[ImageAcquisitionSettingId],
                                "FoilHolePx_x":[PixelCenter_x],
                                "FoilHolePx_y":[PixelCenter_y],
                                "Selected":[selected],
                                "State":[state]})
          # Append to dataframe
          foilHoleMetadataDf = pd.concat([foilHoleMetadataDf, df2], ignore_index=True, sort=False)

          #if selected == "true":
             #aqCount = aqCount + 1

          i = i+1

       # Sanity check foilhose acquisitions
       #print(foilHoleMetadataDf)

       # Count how many aquisitions on the foil holes were selected for collection
       #aqCount = foilHoleMetadataDf['Selected'].str.contains('true').value_counts()[True]
       #aqCount = len(df[foilHoleMetadataDf['Selected'].str.contains('true')])

       # This is actually how many foil hoels were selected
       aqCount = sum(foilHoleMetadataDf['Selected'] == 'true')

    # if no foil holes selected, return np.nan
    if empty == 'Y':
       foilAqNumber = np.nan
    else:
       #foilAqNumber = len(dataAqList)
       foilAqNumber = aqCount

    #foilHoleName = data["TargetLocations"]["TargetLocationsEfficient"]["a:m_serializationArray"]["b:KeyValuePairOfintTargetLocationXmlBpEWF4JT"][0]["b:key"]

    return [squareName, foilAqNumber]

def xml_template_dataframe(xml_path: Path) -> Dict[str, Any]:
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    # Get hole size (do so independently in case running this function alone)
    holeSize = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["HoleSize"]
    holeSizeMicron = float(holeSize)*1e6

    holeSpace = data["Samples"]["_items"]["SampleXml"][0]["FilterHolesSettings"]["HoleSpacing"]
    holeSpaceMicron = float(holeSpace)*1e6/2

    # Match measured hole size to regularised hole size
    holeSizeRegular, holeSpaceRegular = holeRegularisation(holeSizeMicron, holeSpaceMicron)

    # Hole template image pixel size, converted into microns per pixel
    templateScaleH = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["TemplateImagePixelSize"]["a:height"]
    templateScaleW = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["TemplateImagePixelSize"]["a:width"]
    #templateApxH = float(templateScaleH)*1e9
    #templateApxW = float(templateScaleW)*1e9
    templateMpxH = float(templateScaleH)*1e6
    templateMpxW = float(templateScaleW)*1e6

    # Get shots per hole
    shotsPerHole = xml_templateShotsPerHole(xml_path)
    xml_template_dataframe.shotsPerHole = shotsPerHole

    print('Number of shots per hole: '+str(shotsPerHole))
    print('')

    # Calculate the data acquitision exposure dimensions in microns
    cameraXY = getCameraDim(search_mics.firstXml)
    cameraXmicron = cameraXY[0]
    cameraYmicron = cameraXY[1]

    # Create list for gathering template coordinates for plotting
    proposalList = []
    visitList = []
    sessionList = []
    templateNoList = []
    holeSizeList = []
    holeSizeRegList = []
    holeSpaceRegList = []
    beamDiameterList = []
    cameraYMicronList = []
    cameraXMicronList = []
    templateNameList = []
    templateXMicronList = []
    templateYMicronList = []

    if shotsPerHole == 1: # Is single shot
        # Extract the shot name, shot coordinate in pixels on template image
        templateName = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:key"]
        templateHpx = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ShiftInPixels"]["c:height"]
        templateWpx = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"]["b:value"]["ShiftInPixels"]["c:width"]
        # Convert to micron coordinates
        templateHmicron = float(templateHpx)*float(templateMpxH)
        templateWmicron = float(templateWpx)*float(templateMpxW)
        # Report to terminal shot name and coordinates
        print(str('1: Template '+str(templateName)+' (x, y µm): '+str(templateWmicron)+', '+str(templateHmicron)))

        # Quirk, you need to define x so correctly it is used later for numbering autofocus shot in dataframe
        x = 0

        # Append to list of shotID and coordinates
        proposalList.append(lookupBAG.proposal)
        visitList.append(main.visitNameDLS)
        sessionList.append(xml_session.sessionName)
        templateNoList.append(str(round(float(x)+1)))
        holeSizeList.append(holeSizeMicron)
        holeSizeRegList.append(holeSizeRegular)
        holeSpaceRegList.append(holeSpaceRegular)
        beamDiameterList.append(xml_presets.beamD)
        cameraYMicronList.append(cameraYmicron)
        cameraXMicronList.append(cameraXmicron)
        templateNameList.append(templateName)
        templateYMicronList.append(templateHmicron)
        templateXMicronList.append(templateWmicron)

    else: # Is multishot
        # Loop through multishot templates to calculate real coordinates in microns and add to template plot
        for x in range(0, shotsPerHole):
            # Extract the shot name, shot coordinate in pixels on template image
            templateName = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][x]["b:key"]
            templateHpx = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][x]["b:value"]["ShiftInPixels"]["c:height"]
            templateWpx = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["DataAcquisitionAreas"]["a:m_serializationArray"]["b:KeyValuePairOfintDataAcquisitionAreaXmlBpEWF4JT"][x]["b:value"]["ShiftInPixels"]["c:width"]
            # Convert to micron coordinates
            templateHmicron = float(templateHpx)*float(templateMpxH)
            templateWmicron = float(templateWpx)*float(templateMpxW)
            # Report to terminal shot name and coordinates
            print(str(round(float(x)+1))+': Template '+str(templateName)+' (x, y µm): '+str(templateWmicron)+', '+str(templateHmicron))

            # Append to list of shotID and coordinates
            proposalList.append(lookupBAG.proposal)
            visitList.append(main.visitNameDLS)
            sessionList.append(xml_session.sessionName)
            templateNoList.append(str(round(float(x)+1)))
            holeSizeList.append(holeSizeMicron)
            holeSizeRegList.append(holeSizeRegular)
            holeSpaceRegList.append(holeSpaceRegular)
            beamDiameterList.append(xml_presets.beamD)
            cameraYMicronList.append(cameraYmicron)
            cameraXMicronList.append(cameraXmicron)
            templateNameList.append(templateName)
            templateYMicronList.append(templateHmicron)
            templateXMicronList.append(templateWmicron)

    # Extract the autofocus, shot coordinate in pixels on template image
    autoFocusName = 'autofocus'
    autoFocusHpx = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["AutoFocusArea"]["ShiftInPixels"]["a:height"]
    autoFocusWpx = data["Samples"]["_items"]["SampleXml"][0]["TargetAreaTemplate"]["AutoFocusArea"]["ShiftInPixels"]["a:width"]
    # Convert to micron coordinates
    autoFocusHmicron = float(autoFocusHpx)*float(templateMpxH)
    autoFocusWmicron = float(autoFocusWpx)*float(templateMpxW)
    # Report to terminal shot name and coordinates
    print(str(round(float(shotsPerHole+1)))+': Auto Focus (x, y µm): '+str(autoFocusWmicron)+', '+str(autoFocusHmicron))

    # Append to list of shotID and coordinates
    proposalList.append(lookupBAG.proposal)
    visitList.append(main.visitNameDLS)
    sessionList.append(xml_session.sessionName)
    templateNoList.append(str(round(float(x)+2)))
    holeSizeList.append(holeSizeMicron)
    holeSizeRegList.append(holeSizeRegular)
    holeSpaceRegList.append(holeSpaceRegular)
    beamDiameterList.append(xml_presets.beamD)
    cameraYMicronList.append(cameraYmicron)
    cameraXMicronList.append(cameraXmicron)
    templateNameList.append(autoFocusName)
    templateYMicronList.append(autoFocusHmicron)
    templateXMicronList.append(autoFocusWmicron)

    # Initialise shot lists into dataframe
    xml_template_dataframe.dfShot = pd.DataFrame(
        {
            'proposal': proposalList,
            'visit': visitList,
            'EPU session date': xml_session.sessionDate,
            'session': sessionList,
            'shotNo': templateNoList,
            'shotID': templateNameList,
            'shotCoordXmicron': templateXMicronList,
            'shotCoordYmicron': templateYMicronList,
            'shotDimXmicron': cameraXMicronList,
            'shotDimYmicron': cameraYMicronList,
            'beamDiameterMicron': beamDiameterList,
            'holeSizeMicron': holeSizeList,
            'holeSizeRegularised': holeSizeRegList,
            'holeSpaceRegularised': holeSpaceRegList
            }
    )

    # Write out csv with shot information
    # Look in analyseProcessedPerShot function for how analysis statistics per shot are addded to the shots.csv
    print('Writing shot per hole information to csv')
    xml_template_dataframe.dfShot.to_csv(main.csv_dir+'/'+xml_session.sessionName+'_shots.csv', index = False, header=True)

    # Global lists
    xml_template_dataframe.templateNameList = templateNameList

def xml_template_dataframe_plot(df):
    print('Plotting shot templates to graphical representation')
    print()

    # Make a plot of the hole, using the measured hole diameter to input radius of circle
    # https://www.codespeedy.com/how-to-draw-shapes-in-matplotlib-with-python/
    fig2 = plt.figure(2)
    ax1 = fig2.add_subplot()
    plt.gcf().set_size_inches(4, 4)

    # Plot the hole circle based on the measured hole size
    #holeSizeMicron = df['holeSizeMicron'].iloc[0]
    holeSizeMicron = df['holeSizeRegularised'].iloc[0]
    holeSpaceMicron = df['holeSpaceRegularised'].iloc[0]
    circle = plt.Circle((0,0),float(holeSizeMicron)/2, fc='lightgrey',ec="black")
    plt.gca().add_patch(circle)
    ax1.axis('scaled')

    # Loop through the shot information dataframe to plot shot templates and autofocus template
    for index, row in df.iterrows():
        shotID = row['shotID']
        shotNo = row['shotNo']
        templateWmicron = row['shotCoordXmicron']
        templateHmicron = row['shotCoordYmicron']
        cameraXmicron = row['shotDimXmicron']
        cameraYmicron = row['shotDimYmicron']
        templateBeamDiameter = row['beamDiameterMicron']

        if row['shotID'] == 'autofocus':
            color = 'b'
        else:
            color ='g'

        ## Add the templates to the plot using the dataframe as input
        # Beam
        circle = plt.Circle((float(templateWmicron), float(templateHmicron)), float(templateBeamDiameter)/2, color=color, alpha = 0.3)
        # Sensor
        sensor = plt.Rectangle((float(templateWmicron)-(cameraXmicron/2), float(templateHmicron)-(cameraYmicron/2)), cameraXmicron, cameraYmicron, fc=color,ec=color, alpha=0.2)
        # Label
        #label = plt.text(float(templateWmicron)-0.04, float(templateHmicron)-0.04, str(shotNo), style='oblique', color=color, fontsize=6)
        label = plt.text(float(templateWmicron)-0.25, float(templateHmicron)-0.04, str(shotID), style='oblique', color=color, fontsize=6)
        # Add to plot
        ax1 = plt.gca()
        ax1.add_patch(circle)
        ax1.add_patch(sensor)

    # Show and save the plot
    ax1.set_xlim((-2, 2))
    ax1.set_ylim((-2, 2))
    ax1.set_facecolor('silver')
    ax1.set_title('Template definition (R'+str(holeSizeMicron)+'/'+str(holeSpaceMicron)+")\n"+xml_session.sessionName,fontweight = "bold", fontsize=8)
    #ax1.set_title('Template definition (R'+str(xml_targets.holeSize)+'/'+str(xml_targets.holeSpace)+")\n"+xml_session.sessionName,fontweight = "bold", fontsize=8)
    ax1.set_ylabel('Hole image Y coordinate (microns)', fontsize=5)
    ax1.set_xlabel('Hole image X coordinate (microns)', fontsize=5)
    ax1.tick_params(axis='both', which='major', labelsize=5)
    ax1.tick_params(axis='both', which='minor', labelsize=5)
    plt.tight_layout()
    if not args.report:
        print('No report, not saving plot')
    else:
        fig2.savefig(main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate.png',dpi=300)
    #plt.show()
    plt.figure(2).clear()
    plt.close(2)

def xml_template_dataframe_plot_stats(df, column, cmapname, unit):
    print('Plotting shot templates to graphical representation with data from: '+str(column))
    print()

    suffix = column.strip("_") # Strip out '_' prefix for relion parameters

    # Make a plot of the hole, using the measured hole diameter to input radius of circle
    # https://www.codespeedy.com/how-to-draw-shapes-in-matplotlib-with-python/
    fig2 = plt.figure(2)
    ax1 = fig2.add_subplot()
    plt.gcf().set_size_inches(4.5, 5)

    # Plot the hole circle based on the measured hole size
    holeSizeMicron = df['holeSizeMicron'].iloc[0]
    circle = plt.Circle((0,0),float(holeSizeMicron)/2, fc='lightgrey',ec="black")
    plt.gca().add_patch(circle)
    ax1.axis('scaled')

    # Apply color map (cmap) to values in 'column' and return color values in df['cmap_color']
    df = applyCmapDf(df, column, cmapname)

    # Loop through the shot information dataframe to plot shot templates and autofocus template
    for index, row in df.iterrows():
        shotID = row['shotID']
        shotNo = row['shotNo']
        label = row[column]

        templateWmicron = row['shotCoordXmicron']
        templateHmicron = row['shotCoordYmicron']
        cameraXmicron = row['shotDimXmicron']
        cameraYmicron = row['shotDimYmicron']
        templateBeamDiameter = row['beamDiameterMicron']

        if row['shotID'] == 'autofocus':
            color = 'grey'
            label = 'autofocus'
            label_color = 'grey'
            adjust = 0.25
        else:
            color = row['cmap_color']
            #label = str(label)+' '+str(unit)
            label = str(label)
            label_color = 'black'
            adjust = 0.1

        ## Add the templates to the plot using the dataframe as input
        # Beam
        circle = plt.Circle((float(templateWmicron), float(templateHmicron)), float(templateBeamDiameter)/2, color=color, alpha = 0.3)
        # Sensor
        sensor = plt.Rectangle((float(templateWmicron)-(cameraXmicron/2), float(templateHmicron)-(cameraYmicron/2)), cameraXmicron, cameraYmicron, fc=color,ec=color, alpha=0.2)
        # Label
        #label = plt.text(float(templateWmicron)-0.04, float(templateHmicron)-0.04, str(shotNo), style='oblique', color=color, fontsize=6)
        label = plt.text(float(templateWmicron)-adjust, float(templateHmicron)-0.04, str(label), color=label_color, style='oblique', fontsize=6)
        # Add to plot
        ax1 = plt.gca()
        ax1.add_patch(circle)
        ax1.add_patch(sensor)

    ## Create a color map legend
    # Find min and max to assign these to min and max of cmap
    min_val = df[column].min()
    max_val = df[column].max()
    sm = plt.cm.ScalarMappable(cmap=cmapname, norm=plt.Normalize(vmin=min_val, vmax=max_val))
    sm._A = []  # Empty array to fool the colorbar into thinking there's data
    # Create a divider for existing axes instance
    divider = make_axes_locatable(ax1)
    # Define the size and position of the colorbar legend
    cax = divider.append_axes("bottom", size="5%", pad=0.5)  # Adjust 'size' and 'pad' as needed

    # Add colorbar legend
    cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label(str(unit), size=6)  # Set the label for the colorbar

    # Set the font size of the colorbar tick labels
    cbar.ax.tick_params(labelsize=6)

    # Show and save the plot
    ax1.set_xlim((-2, 2))
    ax1.set_ylim((-2, 2))
    ax1.set_facecolor('silver')
    #ax1.set_title('Shot per hole average values for: '+str(column)+"\n"+xml_session.sessionName,fontweight = "bold", fontsize=8)
    ax1.set_title(str(column)+"\n"+xml_session.sessionName,fontweight = "bold", fontsize=8)
    ax1.set_ylabel('Hole image Y coordinate (microns)', fontsize=5)
    ax1.set_xlabel('Hole image X coordinate (microns)', fontsize=5)
    ax1.tick_params(axis='both', which='major', labelsize=5)
    ax1.tick_params(axis='both', which='minor', labelsize=5)
    plt.tight_layout()
    if not args.report:
        print('No report, not saving plot')
    else:
        fig2.savefig(main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate_'+str(suffix)+'.png',dpi=300)
    #plt.show()
    plt.figure(2).clear()
    plt.close(2)

def applyCmapDf(df, column, cmapname):
    # Assign values in dataframe [column] to color map for plotting
    cmap_name = cmapname  # You can use other colormaps like 'plasma', 'inferno', etc.
    cmap = get_cmap(cmap_name)
    # Find min and max to assign these to min and max of cmap
    min_val = df[column].min()
    max_val = df[column].max()
    print('Normalising color in column '+str(column)+' to range '+str(min_val)+'-'+str(max_val))
    print()
    # Normalise values to apply across scale of cmap
    df['normalised'] = (df[column] - min_val) / (max_val - min_val)
    # Apply cmap to normalised values
    df['cmap_color'] = df['normalised'].apply(cmap)

    # Some data is colored as an outlier, brutish filter out here
    black_threshold=0.1
    # Define a function to check if a color is close to black
    def is_close_to_black(rgba_tuple):
        r, g, b, a = rgba_tuple
        return r < black_threshold and g < black_threshold and b < black_threshold and a == 0.0

    # Filter out rows where cmap_color is close to black
    #df = df[~df['cmap_color'].apply(is_close_to_black)]

    return df

# For some reason the applyCmapDf function leaves the normalised and cmap_color columns in the xml_template_dataframe.dfShot dataframe... have this function in case
def dropCmapDf(df):
    df = df.drop(columns=['normalised', 'cmap_color'])

    return df

def xml_presets(xml_path: Path) -> Dict[str, Any]:
    with open(xml_path, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["EpuSessionXml"]

    ## Presets
    # Loop through the presets in the Microscope Settings list
    presets = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"]
    length=len(presets)
    camera = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][0]["value"]["b:Acquisition"]["c:camera"]["c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"]
    lengthCam=len(camera)

    # Create list for gathering preset conditions for reporting
    # DEV DEV DEV might be more flexible to have these going into dataframe
    namePresetList = []
    probePresetList = []
    magPresetList = []
    apixPresetList = []
    spotPresetList = []
    c2PresetList = []
    beamDPresetList = []
    defocusPresetList = []
    timePresetList = []
    binPresetList = []

    # Loop to gather all microscope presets used for session
    for x in range(0, length):
        name = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["key"]

        # Get magnifications from image xml, they are not stored in the epu session file
        if  name == 'Atlas':
            mag = getXmlMag(searchSupervisorAtlas.xmlAtlas)[0]
            apix = getXmlMag(searchSupervisorAtlas.xmlAtlas)[1]
        elif name == 'GridSquare':
            mag = getXmlMag(searchSupervisorData.xmlSquare)[0]
            apix = getXmlMag(searchSupervisorData.xmlSquare)[1]
        elif name == 'Hole':
            mag = getXmlMag(searchSupervisorData.xmlHole)[0]
            apix = getXmlMag(searchSupervisorData.xmlHole)[1]
        elif name == 'Acquisition':
            mag = getXmlMag(searchSupervisorData.xmlData)[0]
            apix = getXmlMag(searchSupervisorData.xmlData)[1]
        else:
            mag = 0
            apix = 0

        probeMode = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:ProbeMode"]
        spot = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:SpotIndex"]
        c2 = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:Apertures"]["c:C2Aperture"]["c:Diameter"]
        beamD = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:BeamDiameter"]
        # Deal with two condensor lens systems that don't know beam diameter
        if isinstance(beamD, dict):
            beamD = 0
        else:
            beamD = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:BeamDiameter"]
        beamDmicron = float(beamD)*1e6
        DF = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Optics"]["c:Defocus"]
        DFmicron = float(DF)*1e6
        time = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:ExposureTime"]
        epuBin = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:Binning"]["d:x"]

        # Here we face data in lists and not always in the same list position so need to loop to find position
        #for y in range(0, lengthCam):
            #listKeyValue = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"][y]["key"]["#text"]
            #if listKeyValue == 'SuperResolutionFactor':
                #superRes = listKeyValue
                #superResOn = data["Samples"]["_items"]["SampleXml"][0]["MicroscopeSettings"]["KeyValuePairs"]["KeyValuePairOfExperimentSettingsIdMicroscopeSettingsCG2rZ1D8"][x]["value"]["b:Acquisition"]["c:camera"]["c:CameraSpecificInput"]["KeyValuePairs"]["KeyValuePairOfstringanyType"][y]["value"]["#text"]

        # Rounding
        spot = round(float(spot))
        c2 = round(float(c2))
        beamDmicron = roundup(float(beamDmicron),1)
        DFmicron = roundup(float(DFmicron),1)
        time = roundup(float(time),2)
        mag = round(float(mag))
        apix = roundup(float(apix),3)

        # Report presets or continue silentily
        if not args.silent or args.silent == 'N':
            print(name)
            print('Nominal magnification: '+str(mag)+' X')
            print('Pixel size: '+str(apix)+' apix')
            print('Probe mode: '+str(probeMode))
            print('Spot: '+str(spot))
            print('C2 apeture: '+str(c2)+' microns')
            print('Beam diameter: '+str(beamDmicron)+' microns')
            print('Defocus: '+str(DFmicron)+' microns')
            print('Exposure time: '+str(time)+' seconds')
            print('')

        # Append presets to the preset lists for reporting
        namePresetList.append(name)
        magPresetList.append(mag)
        apixPresetList.append(apix)
        probePresetList.append(probeMode)
        spotPresetList.append(spot)
        c2PresetList.append(c2)
        beamDPresetList.append(beamDmicron)
        defocusPresetList.append(DFmicron)
        timePresetList.append(time)
        binPresetList.append(epuBin)

        # Gather main params for reporting
        if  name == 'Acquisition':
            xml_presets.time = time
            xml_presets.beamD = beamDmicron
            xml_presets.probe = probeMode
            xml_presets.C2 = c2
            xml_presets.spot = spot
            xml_presets.epuBin = epuBin
            xml_presets.mag = mag
            xml_presets.apix = apix
            #xml_presets.superRes = superRes
            #xml_presets.superResOn = superResOn
        if  name == 'AutoFocus':
            xml_presets.beamDAutoFocus = beamDmicron

        # Gather all presets for mass reporting
        xml_presets.namePresetList = namePresetList
        xml_presets.magPresetList = magPresetList
        xml_presets.apixPresetList = apixPresetList
        xml_presets.probePresetList = probePresetList
        xml_presets.spotPresetList = spotPresetList
        xml_presets.c2PresetList = c2PresetList
        xml_presets.beamDPresetList = beamDPresetList
        xml_presets.defocusPresetList = defocusPresetList
        xml_presets.timePresetList = timePresetList
        xml_presets.binPresetList = binPresetList

    # report complete
    print('Finished gathering all microscope presets')

def xml_presets_data(micpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(micpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    ## SuperResolutionBinning Factor
    # The SuperResolutionFactor is not always in the same list position in a:KeyValueOfstringanyType
    keyValueList = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"]

    # Loop through the list to find the SuperResolutionFactor list position
    i = 0
    for value in keyValueList:
       key = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][i]["a:Key"]
       if key == "SuperResolutionFactor":
          j = i
       i = i+1

    try:
       superResBin = data["microscopeData"]["acquisition"]["camera"]["CameraSpecificInput"]["a:KeyValueOfstringanyType"][j]["a:Value"]["#text"]
    except:
       superResBin = 'Unknown'

    ## Energy filter
    # Known error in nt29493-49 - glacios
    try:
       filterSlit = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitInserted"]
    except:
       filterSlit = 'None'
    else:
       filterSlit = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitInserted"]
    try:
       filterSlitWidth = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitWidth"]
    except:
       filterSlitWidth = 'None'
    else:
       filterSlitWidth = data["microscopeData"]["optics"]["EnergyFilter"]["EnergySelectionSlitWidth"]

    # Aperture(s)
    # Loop through the list to find the Objective aperture list position
    i = 0
    for value in keyValueList:
       key = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Key"]
       if key == "Aperture[OBJ].Name":
          keyvalue = i
          objectiveAperture = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Value"]["#text"]
       i = i+1

    # Stage tilt
    stageAlphaRad = getStageTilt(micpath)[0]
    stageBetaRad = getStageTilt(micpath)[1]
    stageAlpha = roundup(math.degrees(float(stageAlphaRad)),2)
    stageBeta = roundup(math.degrees(float(stageBetaRad)),2)

    # Report
    xml_presets_data.superResBin = superResBin
    xml_presets_data.filterSlitWidth = filterSlitWidth
    xml_presets_data.filterSlit = filterSlit
    xml_presets_data.stageAlpha = roundup(float(stageAlpha),1)
    xml_presets_data.stageBeta = roundup(float(stageBeta),1)
    xml_presets_data.objective = objectiveAperture

def calc_params():
    #Code in here to deal with no dose input and defaulting to something reasonable based on detector
    print('\033[1m'+'Calculating important parameters from EPU session:'+'\033[0m')
    print('')
    if not args.apix:
        print('Pixel size: '+str(getScopeNameMag.apix)+' Å/px (read from xml)')
    else:
        print('Pixel size: '+str(getScopeNameMag.apix)+' Å/px (entered manually)')

    # Before trying to estimate dose try to read it from the gain preparation script README
    # Remember that args.dose is the input, main.dose is the variable that args.dose is copied into
    if not args.dose:
        # if K3 then try to get dose from gain script readme metadata, if Falcon, look in xml file
        if getScopeNameMag.camera == 'K3':
            # Does the README from the gain processing script exist?
            if os.path.exists(main.visitPath+'/processing/gain/README'):
               print('Found gain readme, will attempt to calibrate the dose rate')
               f = open(main.visitPath+'/processing/gain/README','r')
               search_phrase = "Dose rate over vacuum"# (e/px/s):"
               # Search for phrase and get line number
               line_num = 0
               for line in f.readlines():
                  if line.find(search_phrase) >= 0:
                     #print(line_num)
                     line_found = line_num
                  line_num += 1

               try:
                  line_found
               except:
                  print('Could not find calibrated dose rate in README')
                  print('Dose rate will be estimated based on camera optimal dose (e/px/s)')
                  ## If no dose entered then make a reaonsable estimation of for K3 or for F4i
                  dose = getScopeNameMag.optimalDose
                  doseMessage = 'estimated'
               else:
                  f = open(main.visitPath+'/processing/gain/README','r')
                  # Read the following line to get the dose rate
                  line_num = 0
                  for line in f.readlines():
                     if line_num == line_found+1:
                        dose = line
                        doseMessage = 'calibrated'
                     line_num += 1

            else:
               #print('No dose rate entered (e/px/s)')
               #doseMessage = 'uncalibrated'
               print('Dose rate will be estimated based on camera optimal dose (e/px/s)')
               ## If no dose entered then make a reaonsable estimation of for K3 or for F4i
               dose = getScopeNameMag.optimalDose
               doseMessage = 'estimated'
        elif getScopeNameMag.camera == 'F4' or getScopeNameMag.camera == 'F4i':
            # DEV DEV DEV discovered that does rate is measured for every single micrograph, so cannot pull from data aq
            # Is it stored somewhere else globally for falcon images?
            #dose = xmlDoseRate(search_mics.firstXml)
            #doseMessage = 'calibrated'

            print('Dose rate will be estimated based on camera optimal dose (e/px/s)')
            ## If no dose entered then make a reaonsable estimation of for K3 or for F4i
            dose = getScopeNameMag.optimalDose
            doseMessage = 'estimated'
    elif args.dose == 'e':
	# manually force estimation
        print('Dose rate will be estimated based on camera optimal dose (e/px/s)')
        ## If no dose entered then make a reaonsable estimation of for K3 or for F4i
        dose = getScopeNameMag.optimalDose
        doseMessage = 'estimated'
    elif args.dose == '1':
        dose = args.dose
        doseMessage = 'uncalibrated'
    else:
        print('Dose rate entered (e/px/s)')
        dose = args.dose
        doseMessage = 'manual'

    doseRate = roundup((float(dose)/(float(getScopeNameMag.apix)*float(getScopeNameMag.apix))),3)
    doseTotal = roundup((float(doseRate)*float(xml_presets.time)),1)
    print('Dose rate: '+str(dose)+' e/px/s '+str(doseMessage))
    print('Dose rate: '+str(doseRate)+' e/Å2/s '+str(doseMessage))

    # Frame dose
    fractionNo = getFractionNo(search_mics.firstXml)
    if fractionNo == 'Unknown':
       doseFrameTotal = 'N.D.'
    else:
       doseFrameTotal = roundup(float(doseTotal)/float(fractionNo),1)

    print('')
    print('Exposure time was: '+str(xml_presets.time)+' seconds')
    print('Total dose was: '+str(doseTotal)+' e/Å2')
    print('Total frame dose: '+str(doseFrameTotal)+' e/Å2')

    #Report
    calc_params.dose = dose
    calc_params.doseRate = doseRate
    calc_params.doseTotal = doseTotal
    calc_params.fractionNo = fractionNo
    calc_params.doseFrameTotal = doseFrameTotal
    calc_params.doseMessage = doseMessage

# Currently only designed to look through Falcon 4 xml for dose rate
def xmlDoseRate(xmlpath: Path) -> Dict[str, Any]:
    # This will fetch the first micrograph xml data
    with open(xmlpath, "r") as xml:
        for_parsing = xml.read()
        data = xmltodict.parse(for_parsing)
    data = data["MicroscopeImage"]

    # The values are not always in the same list position in a:KeyValueOfstringanyType
    keyValueList = data["CustomData"]["a:KeyValueOfstringanyType"]

    # Loop through the list to find the DoseRate list position
    i = 0
    for value in keyValueList:
       key = data["CustomData"]["a:KeyValueOfstringanyType"][i]["a:Key"]
       # Krios III has EF selectris falcon, Glacios II has standard Falcon no energy filter, BM is base model?
       if key == "Detectors[BM-Falcon].DoseRate" or key == "Detectors[EF-Falcon].DoseRate":
          keyvalue = i
       i = i+1

    xmlDoseRate = data["CustomData"]["a:KeyValueOfstringanyType"][keyvalue]["a:Value"]["#text"]

    return xmlDoseRate

def construct_session_report(report):
    print('Generating session PDF report...')
    # You need to look into another PDF tool, this is getting silly
    # https://stackoverflow.com/questions/47263621/page-breaks-for-images-in-python-fpdf
    colwidth = 40
    col0width = 65
    col1width = 50
    col2width = 26 # 420 pixels
    colH = 4
    # https://towardsdatascience.com/how-to-create-pdf-reports-with-python-the-essential-guide-c08dd3ebf2ee
    # https://stackoverflow.com/questions/51864730/python-what-is-the-process-to-create-pdf-reports-with-charts-from-a-db
    pdf = FPDF()

    # Add page
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(130, 18)
    pdf.cell(75, 18, "eBIC Session Report (EPU)", 0, 0, 'C')

    # Add headers
    pdf.set_xy(0, 13)
    pdf.set_font('Arial', 'I', 8)
    pdf.cell(27,5)
    pdf.cell(75,5,'Report creation is in beta testing, report errors to kyle.morris@diamond.ac.uk', 0, 0, 'C')

    pdf.set_xy(0, 18)
    pdf.set_font('Arial', 'I', 8)
    pdf.cell(34,5)
    pdf.cell(75,5,'Confidential document: only to be shared with the data owner and Diamond Light Source', 0, 0, 'C')

    # Establish initial position for text
    pdf.set_xy(0, 23)
    pdf.set_font('arial','B', 8)
    pdf.cell(15,0)

    # Begin making report
    pdf.cell(0, 4, " ", 0, 2, 'L')
    pdf.cell(colwidth, colH, 'DLS Visit:', 0, 0, 'L')
    pdf.cell(0, colH, str(main.visitNameDLS), 0, 2, 'L')
    pdf.cell(-colwidth)
    pdf.cell(colwidth, colH, 'BAG (institute):', 0, 0, 'L')
    pdf.cell(0, colH, str(lookupBAG.name)+' ('+str(lookupBAG.institute)+')', 0, 2, 'L')
    pdf.cell(-colwidth, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-colwidth) # Reset column position
    pdf.cell(colwidth, colH, 'EPU session name:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.sessionName), 0, 2, 'L')
    pdf.cell(-colwidth) # Reset column position
    pdf.cell(colwidth, colH, 'EPU session created:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.sessionDate), 0, 2, 'L')
    pdf.cell(-colwidth) # Reset column position
    pdf.cell(colwidth, colH, 'Atlas session name:', 0, 0, 'L')
    pdf.cell(0, colH, str(screeningSession_xml.name), 0, 2, 'L')
    pdf.cell(-colwidth) # Reset column position
    pdf.cell(colwidth, colH, 'Atlas session created:', 0, 0, 'L')
    pdf.cell(0, colH, str(screeningSession_xml.date), 0, 2, 'L')
    pdf.cell(-colwidth) # Reset column position
    pdf.cell(colwidth, colH, 'Session location:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.xmlPath), 0, 2, 'L')
    pdf.cell(-colwidth) # Reset column position
    pdf.cell(colwidth, colH, 'Movie location:', 0, 0, 'L')
    pdf.cell(0, colH, str(findRaw.dir), 0, 2, 'L')
    pdf.cell(-colwidth) # Reset column position
    pdf.cell(colwidth, colH, 'Report creation time (user):', 0, 0, 'L')
    pdf.cell(0, colH, str(main.runDate)+' ('+str(main.user)+')', 0, 2, 'L')
    pdf.cell(-colwidth, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-colwidth) # Reset column position
    pdf.set_font('arial','', 8)

    pdf.set_font('arial','B', 8)
    pdf.cell(col0width, colH, 'Instrument:', 0, 2, 'L')
    pdf.set_font('arial','', 8)
    # Microscope
    pdf.cell(col0width, colH, 'Instrument name:', 0, 0, 'L')
    pdf.cell(0, colH, str(getScopeNameMag.scope+' ('+getScopeNameMag.dlsName+')'), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position

    # Sample
    pdf.cell(col0width, colH, 'Autoloader slot:', 0, 0, 'L')
    pdf.cell(0, colH, str(round(xml_session.autoSlot)), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Autoloader slot name:', 0, 0, 'L')
    pdf.cell(0, colH, str(atlas_xml.atlasName), 0, 2, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width) # Reset column position

    pdf.set_font('arial','B', 8)
    pdf.cell(col0width, colH, 'Collection strategy:', 0, 2, 'L')
    pdf.set_font('arial','', 8)
    pdf.cell(col0width, colH, 'Grid type:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.gridType), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Measured hole size / spacing (µm):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_targets.holeSize)+' / '+str(xml_targets.holeSpace), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Shots per hole:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_template_dataframe.shotsPerHole), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'AFIS/Accurate:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.afisMode), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Stage shift wait (s):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.delayStageShift), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Beam Image shift wait (s):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.delayImageShift), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Defocus list (µm):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.defocus), 0, 2, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width) # Reset column position

    # Illumination
    pdf.set_font('arial','B', 8)
    pdf.cell(col0width, colH, 'Acquisition parameters:', 0, 2, 'L')
    pdf.set_font('arial','', 8)
    pdf.cell(col0width, colH, 'Tilt (degrees):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_presets_data.stageAlpha), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Acquisition spot size:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_presets.spot), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Acquisition C2 (µm):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_presets.C2), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Acquisition Objective (µm):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_presets_data.objective), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Acquisition beam diameter (µm):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_presets.beamD), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Acquisition magnification (X):', 0, 0, 'L')
    pdf.cell(0, colH, str(getScopeNameMag.mag+' '), 0, 2, 'L')
    pdf.cell(-col0width)
    #pdf.cell(-col0width)
    #pdf.cell(col0width, colH, str(xml_presets.superRes), 0, 0, 'L')
    #pdf.cell(0, colH, str(xml_presets.superResOn), 0, 2, 'L')
    #pdf.cell(col0width, colH, "Acquisition binning on camera:", 0, 0, 'L')
    #pdf.cell(0, colH, str(xml_presets_data.superResBin), 0, 2, 'L')
    #pdf.cell(-col0width)
    #pdf.cell(col0width, colH, 'Acquisition exposure time (sec):', 0, 0, 'L')
    #pdf.cell(0, colH, str(xml_presets.time), 0, 2, 'L')
    #pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Fraction no:', 0, 0, 'L')
    pdf.cell(0, colH, str(calc_params.fractionNo), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Filter slit width (eV) (inserted):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_presets_data.filterSlitWidth)+' ('+str(xml_presets_data.filterSlit)+')', 0, 2, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width) # Reset column position

    # Calculated dose conditions
    pdf.set_font('arial','B', 8)
    pdf.cell(col0width, colH, 'Dose calculations:', 0, 2, 'L')
    pdf.set_font('arial','', 8)

    pdf.cell(col0width, colH, 'Acquisition pixel size (Å/px):', 0, 0, 'L')
    if not args.apix:
        pdf.cell(0, colH, str(getScopeNameMag.apix)+' (calibrated)', 0, 2, 'L')
    else:
        pdf.cell(0, colH, str(getScopeNameMag.apix)+' (entered manually)', 0, 2, 'L')
    #pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Measured vacuum dose rate (e/px/s):', 0, 0, 'L')
    pdf.cell(0, colH, str(calc_params.dose)+' ('+str(calc_params.doseMessage)+')', 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Sample dose rate (e/Å2/s):', 0, 0, 'L')
    pdf.cell(0, colH, str(calc_params.doseRate)+' ('+str(calc_params.doseMessage)+')', 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Acquisition exposure time (sec):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_presets.time), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Total dose / Fraction dose (e/Å2):', 0, 0, 'L')
    pdf.cell(0, colH, str(calc_params.doseTotal)+' / '+str(calc_params.doseFrameTotal)+' ('+str(calc_params.doseMessage)+')', 0, 2, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width) # Reset column position

    # Number of mics
    pdf.set_font('arial','B', 8)
    pdf.cell(col0width, colH, 'Collection statistics (calculated at time of report creation):', 0, 2, 'L')
    pdf.set_font('arial','', 8)
    #pdf.cell(col0width, colH, 'Screening time (hrs):', 0, 0, 'L')
    #pdf.cell(0, colH, str(roundup(calcTimings.screenTimeHr,1)), 0, 2, 'L')
    #pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Setup time (hrs):', 0, 0, 'L')
    pdf.cell(0, colH, str(roundup(calcTimings.setupTimeHr,1)), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Collection time (hrs):', 0, 0, 'L')
    pdf.cell(0, colH, str(roundup(calcTimings.timeHr,1)), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Number of squares total / targeted / collected:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_targets.totalSquares)+' / '+str(xml_targets.targetSquares)+' / '+str(xml_targets.collectedSquares), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Average FoilHoles selected per Square (+/-):', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_targets.avSqFoilHoleNo)+' (+/- '+str(xml_targets.stdevSqFoilHoleNo)+')', 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Total FoilHoles selected:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_targets.totalFoilHoles), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Number of micrographs acquired:', 0, 0, 'L')
    pdf.cell(0, colH, str(count_mics.number), 0, 2, 'L')
    pdf.cell(-col0width)
    #pdf.cell(col0width, colH, 'First micrograph time:', 0, 0, 'L')
    #pdf.cell(0, colH, str(count_mics.first), 0, 2, 'L')
    #pdf.cell(-col0width)
    #pdf.cell(col0width, colH, 'Last micrograph time:', 0, 0, 'L')
    #pdf.cell(0, colH, str(count_mics.last), 0, 2, 'L')
    #pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Collection rate (exposure/hr):', 0, 0, 'L')
    pdf.cell(0, colH, str(calcTimings.rate), 0, 2, 'L')
    pdf.cell(-col0width)
    pdf.cell(col0width, colH, 'Collection rate (um^2/hr):', 0, 0, 'L')
    pdf.cell(0, colH, str(calcTimings.rateArea), 0, 2, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width)

    # Advisories
    pdf.set_font('arial','B', 8)
    pdf.cell(col0width, colH, 'Parameter advisories:', 0, 0, 'L')
    pdf.cell(0, colH, str(assessSession.advisoryCount), 0, 2, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width)

    # Presets
    pdf.set_font('arial','B', 8)
    pdf.cell(col0width, colH, 'All Presets:', 0, 2, 'L')
    #pdf.cell(-col0width)
    pdf.set_font('arial','', 8)

    # Populate preset list table
    pdf.set_font('arial','B', 6)
    pdf.cell(22,3, 'Preset', 1, 0, 'R')
    pdf.cell(12,3, 'Mag (X)', 1, 0, 'C')
    pdf.cell(10,3, 'apix', 1, 0, 'C')
    pdf.cell(15,3, 'Probe', 1, 0, 'C')
    pdf.cell(10,3, 'Spot', 1, 0, 'C')
    pdf.cell(10,3, 'C2 (µm)', 1, 0, 'C')
    pdf.cell(15,3, 'Beam (µm)', 1, 0, 'C')
    pdf.cell(18,3, 'Defocus (µm)', 1, 0, 'C')
    pdf.cell(8,3, 'Bin', 1, 0, 'C')
    pdf.cell(12,3, 'Time (s)', 1, 2, 'C')
    pdf.cell(-120)

    # Mass report all presets
    pdf.set_font('arial','', 6)
    j = 0
    for i in xml_presets.namePresetList:
        pdf.cell(22,3, str(xml_presets.namePresetList[j]), 1, 0, 'R')
        pdf.cell(12,3, str(xml_presets.magPresetList[j]), 1, 0, 'C')
        pdf.cell(10,3, str(xml_presets.apixPresetList[j]), 1, 0, 'C')
        pdf.cell(15,3, str(xml_presets.probePresetList[j]), 1, 0, 'C')
        pdf.cell(10,3, str(xml_presets.spotPresetList[j]), 1, 0, 'C')
        pdf.cell(10,3, str(xml_presets.c2PresetList[j]), 1, 0, 'C')
        pdf.cell(15,3, str(xml_presets.beamDPresetList[j]), 1, 0, 'C')
        pdf.cell(18,3, str(xml_presets.defocusPresetList[j]), 1, 0, 'C')
        pdf.cell(8,3, str(xml_presets.binPresetList[j]), 1, 0, 'C')
        pdf.cell(12,3, str(xml_presets.timePresetList[j]), 1, 2, 'C')
        pdf.cell(-120) # Always return to 0 which is pdf.cell sum(x) - final(x)
        j = j+1

    # New page for template
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(130, 18)
    pdf.cell(75, 18, "eBIC Session Report (EPU)", 0, 0, 'C')

    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(20,0)

    # Insert Atlas image
    pdf.cell(-5) # Reset column position
    pdf.cell(col0width, colH, str('Atlas image of grid used for data collection: '), 0, 2, 'L')
    pdf.cell(col0width, colH, str(main.atlasxml), 0, 2, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.image(main.atlasimgPath, x = 26, y = None, w = 55, h = 55, type = '', link = '')
    pdf.cell(None, colH, '', 0, 2, 'L') # Table space
    pdf.cell(2) # Reset column position
    #pdf.cell(-col0width) # Reset column position

    # Insert some random micrographs
    micSizeWidth = 42
    micSizeHeight = micSizeWidth / float(getScopeNameMag.cameraRatio)

    pdf.cell(col0width, colH, 'Randomly selected square and acqusition images:', 0, 2, 'L')
    startx = pdf.get_x()
    starty = pdf.get_y()+colH

    pdf.set_font('arial','', 6)
    pdf.cell(0, micSizeHeight+colH, '', 0, 2, 'L') # Table space
    pdf.cell(micSizeWidth*2+2, colH, str(search_mics.randomSquareName), 0, 0, 'L')
    pdf.cell(micSizeWidth*2+2, colH, str(search_mics.randomSquare1Name), 0, 2, 'L')
    pdf.cell(-micSizeWidth*2-2)
    pdf.cell(micSizeWidth*2+2, colH, str(search_mics.random0Name), 0, 0, 'L')
    pdf.cell(micSizeWidth*2+2, colH, str(search_mics.random1Name), 0, 2, 'L')
    pdf.cell(-micSizeWidth*2-2)
    pdf.cell(0, 45-colH, '', 0, 2, 'L') # Table space
    pdf.cell(micSizeWidth*2+2, colH, str(search_mics.randomSquare2Name), 0, 0, 'L')
    pdf.cell(micSizeWidth*2+2, colH, str(search_mics.randomSquare3Name), 0, 2, 'L')
    pdf.cell(-micSizeWidth*2-2)
    pdf.cell(micSizeWidth*2+2, colH, str(search_mics.random2Name), 0, 0, 'L')
    pdf.cell(micSizeWidth*2+2, colH, str(search_mics.random3Name), 0, 2, 'L')
    pdf.cell(-micSizeWidth*2-2)
    pdf.cell(10, colH, '', 0, 2, 'L') # Table space
    pdf.set_font('arial','', 8)

    # Random micrograph display
    pdf.image(search_mics.randomSquarePath, x = startx, y = starty, w = micSizeWidth, h = micSizeHeight, type = '', link = '')
    pdf.image(search_mics.randomMicPath, x = startx+micSizeWidth, y = starty, w = micSizeWidth, h = micSizeHeight, type = '', link = '')

    pdf.image(search_mics.randomSquare1Path, x = startx+(micSizeWidth*2)+2, y = starty, w = micSizeWidth, h = micSizeHeight, type = '', link = '')
    pdf.image(search_mics.randomMic1Path, x = startx+(micSizeWidth*2)+micSizeWidth+2, y = starty, w = micSizeWidth, h = micSizeHeight, type = '', link = '')

    pdf.image(search_mics.randomSquare2Path, x = startx, y = starty+45, w = micSizeWidth, h = micSizeHeight, type = '', link = '')
    pdf.image(search_mics.randomMic2Path, x = startx+micSizeWidth, y = starty+45, w = micSizeWidth, h = micSizeHeight, type = '', link = '')

    pdf.image(search_mics.randomSquare3Path, x = startx+(micSizeWidth*2)+2, y = starty+45, w = micSizeWidth, h = micSizeHeight, type = '', link = '')
    pdf.image(search_mics.randomMic3Path, x = startx+(micSizeWidth*2)+micSizeWidth+2, y = starty+45, w = micSizeWidth, h = micSizeHeight, type = '', link = '')

    # Create cell autoloader images are in right place
    #pdf.cell(50, colH*20, '', 0, 2, 'L') # Table space
    #pdf.cell(3)

    # Autoloader images
    pdf.cell(7) # Reset column position
    pdf.cell(col0width, colH, str('Autoloader casette atlases:'), 0, 2, 'L')
    #pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    #pdf.cell(-col0width) # Reset column position

    # Insert atlas thumbnails for 1 through 6
    #pdf.cell(col0width, colH, str('Atlas image of all grids from screening session: '+main.atlasxml), 0, 0, 'L')
    #pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    #pdf.cell(None, 5, '', 0, 2, 'L') # Table space

    # Slot names
    pdf.cell(col2width,3, 'Slot 1', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 2', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 3', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 4', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 5', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 6', 1, 2, 'C')
    pdf.cell(-col2width*5)
    empty = main.noAtlas

    startx = pdf.get_x()
    starty = pdf.get_y()

    # loop to get it looking for and populating atlas images
    for i in range(1, 7):
        sampleNo = str('Sample'+str(i))
        # Returns atlas image path
        atlasimg = atlas_image(main.atlas+'/'+sampleNo)[0]
        #atlasimg = str(main.atlas+'/'+sampleNo+'/Atlas/Atlas_1.jpg')
        if atlasimg == 'None':
            insert = empty
        else:
            insert = atlasimg

        pdf.image(insert, x = startx+((i-1)*26), y = starty, w = 26, h = 26, type = '', link = '')

        #try:
            #with open(atlasimg) as f:
                #print(f.readlines())
                # If atlas exists, insert it into PDF
                #pdf.image(atlasimg, x = startx*i+(i-1), y = starty, w = 26, h = 26, type = '', link = '')
        #except IOError:
                #print("No atlas in: "+sampleNo)
                #pdf.image(empty, x = startx*i+(i-1), y = starty, w = 26, h = 26, type = '', link = '')
    pdf.cell(None, 26, '', 0, 2, 'L') # Table space

    # Reserved for names entered into TUI
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 2, 'C')
    pdf.cell(-col2width*5)

    #pdf.cell(col2width,3, 'Name', 1, 0, 'C')
    #pdf.cell(col2width,3, 'Name', 1, 0, 'C')
    #pdf.cell(col2width,3, 'Name', 1, 0, 'C')
    #pdf.cell(col2width,3, 'Name', 1, 0, 'C')
    #pdf.cell(col2width,3, 'Name', 1, 0, 'C')
    #pdf.cell(col2width,3, 'Name', 1, 2, 'C')
    #pdf.cell(-col2width*5)
    pdf.cell(None, 5, '', 0, 2, 'L') # Table space

    # Insert atlas thumbnails for 7 through 12
    pdf.cell(col2width,3, 'Slot 7', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 8', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 9', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 10', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 11', 1, 0, 'C')
    pdf.cell(col2width,3, 'Slot 12', 1, 2, 'C')
    pdf.cell(-col2width*5)
    empty = main.noAtlas
    startx = pdf.get_x()
    starty = pdf.get_y()

    # loop to get it looking for and populating atlas images
    j = 1
    for i in range(7, 13):
        sampleNo = str('Sample'+str(i))
        # Returns atlas image path
        atlasimg = atlas_image(main.atlas+'/'+sampleNo)[0]
        #atlasimg = str(main.atlas+'/'+sampleNo+'/Atlas/Atlas_1.jpg')
        if atlasimg == 'None':
            insert = empty
        else:
            insert = atlasimg

        pdf.image(insert, x = startx+((j-1)*26), y = starty, w = 26, h = 26, type = '', link = '')

        j = j + 1
    pdf.cell(None, 26, '', 0, 2, 'L') # Table space

    # Reserved for names entered into TUI
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 0, 'C')
    pdf.cell(col2width,3, '', 1, 2, 'C')
    pdf.cell(-col2width*5)

    pdf.cell(None, 5, '', 0, 2, 'L') # Table space

    # New page for template
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(130, 18)
    pdf.cell(75, 18, "eBIC Session Report (EPU)", 0, 0, 'C')

    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(20,0)

    # Insert a nice bat chart for eBIC setup time, user setup time, collection time, pipeline time
    pdf.cell(col0width, colH, 'The following plot shows the sessions timings:', 0, 0, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width)
    # Insert template plot
    pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_session_timings.png', x = 40, y = None, w = 112, h = 64, type = '', link = '')
    pdf.cell(50, colH, '', 0, 2, 'L') # Table space

    # Insert template plot
    pdf.cell(col0width, colH, 'The following plot and table show the calibrated hole and template sizes and location:', 0, 0, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width)
    # Insert template plot
    pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate.png', x = 60, y = None, w = 80, h = 80, type = '', link = '')
    pdf.cell(50, colH, '', 0, 2, 'L') # Table space

    # Template information

    # Populate template list table
    #pdf.set_font('arial','B', 6)
    #pdf.cell(5,3, '#', 1, 0, 'C')
    #pdf.cell(col1width,3, 'Template ID', 1, 0, 'C')
    #pdf.cell(col1width,3, 'x (µm)', 1, 0, 'C')
    #pdf.cell(col1width,3, 'y (µm)', 1, 2, 'C')
    #pdf.cell(float(col1width*-2-5))

    #pdf.set_font('arial','', 6)
    #j = 0
    #for i in xml_template.nameTemplateList:
    #    pdf.cell(5,3, str(j+1), 1, 0, 'C')
    #    pdf.cell(col1width,3, str(xml_template.nameTemplateList[j]+' (green)'), 1, 0, 'C')
    #    pdf.cell(col1width,3, str(xml_template.xTemplateList[j]), 1, 0, 'C')
    #    pdf.cell(col1width,3, str(xml_template.yTemplateList[j]), 1, 2, 'C')
    #    pdf.cell(float(col1width*-2-5))
    #    j = j+1

    # Populate autofocus list
    #pdf.cell(5,3, '-', 1, 0, 'C')
    #pdf.cell(col1width,3, 'Autofocus (blue)', 1, 0, 'C')
    #pdf.cell(col1width,3, str(xml_template.xAutofocus), 1, 0, 'C')
    #pdf.cell(col1width,3, str(xml_template.yAutofocus), 1, 2, 'C')
    #pdf.cell(float(col1width*-2-5))
    #pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    #pdf.cell(-col0width/6)

    #
    pdf.set_font('arial','', 8)
    pdf.cell(col0width, colH, 'The following plot shows the exposures over time:', 0, 0, 'L')
    pdf.cell(-col0width)
    # Insert exposure plot
    pdf.cell(col1width*-2, colH, '', 0, 2, 'L') # Table space
    pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_exposures.png', x = 20, y = None, w = 140, h = 80, type = '', link = '')

    # New page for template
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(130, 18)
    pdf.cell(75, 18, "eBIC Session Report (EPU)", 0, 0, 'C')

    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(20,0)

    # Add advisories as table if any are present
    advisorycsv=(main.csv_dir+'/'+xml_session.sessionName+'_advisory.csv')

    if os.path.exists(advisorycsv):
       pdf.cell(col0width, colH, 'The following tabulates any identified session advisories:', 0, 0, 'L')
       pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
       pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
       pdf.cell(-col0width)
       # Read in csv, keep NaN so you can handle these as strings, force all floats to read as string
       df = pd.read_csv(advisorycsv, keep_default_na=False, na_values=[''])
       row = df.shape[0]  # Gives number of rows
       col = df.shape[1]  # Gives number of columns

       # This dataframe is used to create the global report table, you can add new lines to lookup data in the csv files
       reportTable = [{'name':'Proposal', 'unit':'', 'data':'Proposal', 'width':18, 'height':3, 'line':'1', 'return':'0', 'align':'R'},
                     {'name':'Visit', 'unit':'', 'data':'Visit', 'width':18, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                     {'name':'Advisory type', 'unit':'', 'data':'Advisory_type', 'width':18, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                     {'name':'Advisory', 'unit':'', 'data':'Advisory', 'width':100, 'height':3, 'line':'1', 'return':'2', 'align':'L'}]
       reportDf = pd.DataFrame(reportTable)

       # Reset needs to be total cell width, less one cell width
       widthSum = reportDf['width'].sum()
       widthLast = reportDf['width'].iloc[-1]
       reset = (widthSum-widthLast)*-1

       # Create header for advisory list table
       pdf.set_font('arial','B', 6)
       print('Gathering all advisory')

       # loop through dataframe to create all session column headers (units)
       j = 0
       for index, row in reportDf.iterrows():
          pdf.cell(int(reportDf.iloc[j]['width']),int(reportDf.iloc[j]['height']),str(reportDf.iloc[j]['unit']),int(reportDf.iloc[j]['line']),int(reportDf.iloc[j]['return']),str(reportDf.iloc[j]['align']))
          j = j+1
       pdf.cell(reset)

       # loop through dataframe to create all session column headers (names)
       j = 0
       for index, row in reportDf.iterrows():
          pdf.cell(int(reportDf.iloc[j]['width']),int(reportDf.iloc[j]['height']),str(reportDf.iloc[j]['name']),int(reportDf.iloc[j]['line']),int(reportDf.iloc[j]['return']),str(reportDf.iloc[j]['align']))
          j = j+1
       pdf.cell(reset)

       # Prepare to populate advisory list table
       # Autopopulated by lookup from dataframe, of column name and width, see header method
       # Mass report all data
       pdf.set_font('arial','', 6)

       h = 0
       for index, row in df.iterrows():
         j = 0
         for index, row in reportDf.iterrows():
            dataName = reportDf.iloc[j]['data']
            colName = reportDf.iloc[j]['name']
            data = str(df.iloc[h][dataName])
            #print(data+' '+str(data.isnumeric())+' '+str())
            insert = str(data)
            pdf.cell(int(reportDf.iloc[j]['width']),int(reportDf.iloc[j]['height']),insert,int(reportDf.iloc[j]['line']),int(reportDf.iloc[j]['return']),str(reportDf.iloc[j]['align']))
            j = j+1
         pdf.cell(reset)
         h = h+1

    # Write out report
    pdf.output(report, 'F')

    # Report
    print('Generated PDF report in '+report)
    print()

    #print(search_mics.randomMicPath)
    #print(search_mics.randomSquarePath)

def construct_process_report(report):
    print('Generating processed PDF report...')

    # You need to look into another PDF tool, this is getting silly
    # https://stackoverflow.com/questions/47263621/page-breaks-for-images-in-python-fpdf
    colwidth = 40
    col0width = 55
    col1width = 50
    col2width = 26 # 420 pixels
    col3width = 60
    colH = 4
    # https://towardsdatascience.com/how-to-create-pdf-reports-with-python-the-essential-guide-c08dd3ebf2ee
    # https://stackoverflow.com/questions/51864730/python-what-is-the-process-to-create-pdf-reports-with-charts-from-a-db

    pdf = FPDF()

    # New page for template
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(130, 18)
    pdf.cell(75, 18, "eBIC Session Report (EPU)", 0, 0, 'C')

    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(20,0)

    pdf.set_font('arial','B', 8)
    pdf.cell(0, 8, " ", 0, 2, 'L')
    pdf.cell(colwidth, colH, 'DLS Visit:', 0, 0, 'L')
    pdf.cell(0, colH, str(main.visitNameDLS), 0, 2, 'L')
    pdf.cell(-colwidth)
    pdf.cell(colwidth, colH, 'BAG (institute):', 0, 0, 'L')
    pdf.cell(0, colH, str(lookupBAG.name)+' ('+str(lookupBAG.institute)+')', 0, 2, 'L')
    pdf.cell(-colwidth, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-colwidth)
    pdf.cell(colwidth, colH, 'EPU session name:', 0, 0, 'L')
    pdf.cell(0, colH, str(xml_session.sessionName), 0, 2, 'L')
    pdf.cell(-colwidth)
    pdf.cell(colwidth, colH, 'Movie location:', 0, 0, 'L')
    pdf.cell(0, colH, str(findRaw.dir), 0, 2, 'L')
    pdf.cell(-colwidth) # Reset column position
    pdf.cell(colwidth, colH, 'Report creation time (user):', 0, 0, 'L')
    pdf.cell(0, colH, str(main.runDate)+' ('+str(main.user)+')', 0, 2, 'L')
    pdf.cell(-colwidth, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-colwidth) # Reset column position
    pdf.set_font('arial','', 8)

    # Pipeline
    pdf.cell(col0width, colH, 'Processed directory:', 0, 0, 'L')
    pdf.cell(0, colH, str(searchProcessed.path), 0, 2, 'L')
    pdf.cell(-col0width) # Reset column position
    pdf.cell(col0width, colH, 'Number of pipeline submissions:', 0, 0, 'L')
    pdf.cell(0, colH, str(searchProcessed.number), 0, 2, 'L')
    pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col0width) # Reset column position

    pdf.cell(col3width, colH, 'Number of motion corrected star lines:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.motionCorNo), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position

    pdf.cell(col3width, colH, 'Mean motion (A - early / late / total):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.earlyMotionMean)+' / '+str(analyseProcessed.lateMotionMean)+' / '+str(analyseProcessed.totMotionMean), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Min motion (A - early / late / total):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.earlyMotionMin)+' / '+str(analyseProcessed.lateMotionMin)+' / '+str(analyseProcessed.totMotionMin), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Max motion (A - early / late / total):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.earlyMotionMax)+' / '+str(analyseProcessed.lateMotionMax)+' / '+str(analyseProcessed.totMotionMax), 0, 2, 'L')
    pdf.cell(-col3width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col3width) # Reset column position

    pdf.cell(col3width, colH, 'Number of CTFFIND star lines:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.ctfFindNo), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Min CTF resolution fit (A):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.ctfMin), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Max CTF resolution fit (A):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.ctfMax), 0, 2, 'L')
    pdf.cell(-col3width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col3width) # Reset column position

    # Relion specific - particles
    pdf.cell(col3width, colH, 'Number of particle star lines:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.particleNo), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Minimum pick diameter (ang):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.pickMin), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Maximum pick diameter (ang):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.pickMax), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Extracted box (px):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.extractPx), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Extracted box (ang):', 0, 0, 'L')
    pdf.cell(0, colH, str(round(float(analyseProcessed.extractAng))), 0, 2, 'L')
    pdf.cell(-col3width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col3width) # Reset column position

    # Particle density and clustering analysis
    pdf.cell(col3width, colH, 'Average particles per micrograph:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.particlesPerMicAve), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'StDev particles per micrograph:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.particlesPerMicStd), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Average particles per micrograph (norm):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.particlesPerMicAveNorm), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'StDev particles per micrograph (norm):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.particlesPerMicStdNorm), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position

    pdf.cell(col3width, colH, 'Average pct of particles per mic in cluster:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.particlesPerMicClusterAve), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'StDev pct of particles per mic in cluster:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.particlesPerMicClusterStd), 0, 2, 'L')
    pdf.cell(-col3width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col3width) # Reset column position

    # Relion specific - 2D
    pdf.cell(col3width, colH, 'Mask (ang):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.ptclDAng), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Max probability distribution min:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.maxProbMin), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, 'Max probability distribution max:', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.maxProbMax), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position

    pdf.cell(col3width, colH, '2D class resolution (min):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.class2dresMin), 0, 2, 'L')
    pdf.cell(-col3width) # Reset column position
    pdf.cell(col3width, colH, '2D class resolution (max):', 0, 0, 'L')
    pdf.cell(0, colH, str(analyseProcessed.class2dresMax), 0, 2, 'L')
    pdf.cell(-col3width, colH, '', 0, 2, 'L') # Table space
    pdf.cell(-col3width) # Reset column position

    ############################
    ## Dataset statistics
    ############################

    # New page for template
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(130, 18)
    pdf.cell(75, 18, "eBIC Session Report (EPU)", 0, 0, 'C')

    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(20,0)

    # Insert motion and picking histograms
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)

    row = 80
    # row1 col0
    col = 5
    rowN = 0
    pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_motion_early.png', x = col, y = starty+(row*rowN), w = 100, h = 80, type = '', link = '')

    # row1 col1
    col = 105
    rowN = 0
    pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_motion_late.png', x = col, y = starty+(row*rowN), w = 100, h = 80, type = '', link = '')

    # row2 col0
    col = 5
    rowN = 1
    pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_motion_total.png', x = col, y = starty+(row*rowN), w = 100, h = 80, type = '', link = '')

    # row2 col1
    col = 105
    rowN = 1
    pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_ctf.png', x = col, y = starty+(row*rowN), w = 100, h = 80, type = '', link = '')

    # row3 col0
    col = 15
    rowN = 2
    pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate.png', x = col, y = starty+(row*rowN), w = 80, h = 80, type = '', link = '')

    # row3 col1
    col = 105
    rowN = 2
    #pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate.png', x = col, y = starty+(row*rowN), w = 80, h = 80, type = '', link = '')

    ############################
    ## Dataset statistics shown by shot
    ############################

    # New page for template
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(130, 18)
    pdf.cell(75, 18, "eBIC Session Report (EPU)", 0, 0, 'C')

    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(20,0)

    # Insert motion and picking histograms
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)

    shotTemplateW = 60
    shotTemplateH = 67
    shotBoxPlotW = 97
    shotBoxPlotH = 67

    # Row spacing
    row = 63

    # row1 col0
    col = 10
    rowN = 0
    #path = main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate.png'
    #if os.path.exists(path):
    #    pdf.image(path, x = col, y = starty+(row*rowN), w = 64, h = 64, type = '', link = '')
    path = main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate_'+'rlnCtfMaxResolution.png'
    if os.path.exists(path):
        pdf.image(path, x = col, y = starty+(row*rowN), w = shotTemplateW, h = shotTemplateH, type = '', link = '')

    # row1 col1
    col = 90
    rowN = 0
    path = main.plot_dir+'/shotAnalysis/'+xml_session.sessionName+'_shotAnalysis_'+'rlnCtfMaxResolution'+'_boxplot.png'
    if os.path.exists(path):
        pdf.image(path, x = col, y = starty+(row*rowN), w = shotBoxPlotW, h = shotBoxPlotH, type = '', link = '')

    # row2 col0
    #col = 10
    #rowN = 1
    #path = main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate_'+'rlnAccumMotionEarly.png'
    #if os.path.exists(path):
    #    pdf.image(path, x = col, y = starty+(row*rowN), w = shotTemplateW, h = shotTemplateH, type = '', link = '')

    # row2 col1
    #col = 90
    #rowN = 1
    #path = main.plot_dir+'/shotAnalysis/'+xml_session.sessionName+'_shotAnalysis_'+'rlnAccumMotionEarly'+'_boxplot.png'
    #if os.path.exists(path):
    #    pdf.image(path, x = col, y = starty+(row*rowN), w = shotBoxPlotW, h = shotBoxPlotH, type = '', link = '')

    # row3 col0
    #col = 10
    #rowN = 2
    #path = main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate_'+'rlnAccumMotionTotal.png'
    #if os.path.exists(path):
    #    pdf.image(path, x = col, y = starty+(row*rowN), w = shotTemplateW, h = shotTemplateH, type = '', link = '')

    # row3 col1
    #col = 90
    #rowN = 2
    #path = main.plot_dir+'/shotAnalysis/'+xml_session.sessionName+'_shotAnalysis_'+'rlnAccumMotionTotal'+'_boxplot.png'
    #if os.path.exists(path):
    #    pdf.image(path, x = col, y = starty+(row*rowN), w = shotBoxPlotW, h = shotBoxPlotH, type = '', link = '')

    # row4 col0
    col = 10
    rowN = 1
    path = main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate_'+'pct_clustered.png'
    if os.path.exists(path):
        pdf.image(path, x = col, y = starty+(row*rowN), w = shotTemplateW, h = shotTemplateH, type = '', link = '')

    # row4 col1
    col = 90
    rowN = 1
    path = main.plot_dir+'/shotAnalysis/'+xml_session.sessionName+'_shotAnalysis_'+'pct_clustered'+'_boxplot.png'
    if os.path.exists(path):
        pdf.image(path, x = col, y = starty+(row*rowN), w = shotBoxPlotW, h = shotBoxPlotH, type = '', link = '')

    # row5 col0
    col = 10
    rowN = 2
    path = main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate_'+'ideal_picking.png'
    if os.path.exists(path):
        pdf.image(path, x = col, y = starty+(row*rowN), w = shotTemplateW, h = shotTemplateH, type = '', link = '')

    # row5 col1
    col = 90
    rowN = 2
    path = main.plot_dir+'/shotAnalysis/'+xml_session.sessionName+'_shotAnalysis_'+'ideal_picking'+'_boxplot.png'
    if os.path.exists(path):
        pdf.image(path, x = col, y = starty+(row*rowN), w = shotBoxPlotW, h = shotBoxPlotH, type = '', link = '')

    # row5 col0
    col = 10
    rowN = 3
    path = main.plot_dir+'/'+xml_session.sessionName+'_shotTemplate_'+'Class2dRes_mean.png'
    if os.path.exists(path):
        pdf.image(path, x = col, y = starty+(row*rowN), w = shotTemplateW, h = shotTemplateH, type = '', link = '')

    # row5 col1
    col = 90
    rowN = 3
    path = main.plot_dir+'/shotAnalysis/'+xml_session.sessionName+'_shotAnalysis_'+'Class2dRes_mean'+'_boxplot.png'
    if os.path.exists(path):
        pdf.image(path, x = col, y = starty+(row*rowN), w = shotBoxPlotW, h = shotBoxPlotH, type = '', link = '')

    ############################
    ## Picking representations
    ############################

    # New page for template
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(130, 18)
    pdf.cell(75, 18, "eBIC Session Report (EPU)", 0, 0, 'C')

    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(20,0)

    # Insert some random micrographs
    micSizeWidth = 60
    micSizeHeight = micSizeWidth / float(getScopeNameMag.cameraRatio)

    # Calculate image insertion size for pick plots
    if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random1Name)+'_nnclustering_particles.png'):
        image = Image.open(main.plot_dir+'/particles/'+str(search_mics.random1Name)+'_nnclustering_particles.png')
        width, height = image.size
        ratio = width/height
    else:
        ratio = 1
    pickSizeWidth = 64
    pickSizeHeight = pickSizeWidth / ratio

    # Calculate image insertion size for Class2D graph picks
    if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random1Name)+'_2Dclass_res.png'):
        image = Image.open(main.plot_dir+'/particles/'+str(search_mics.random1Name)+'_2Dclass_res.png')
        width, height = image.size
        ratio = width/height
    else:
        ratio = 1
    class2dSizeWidth = 64
    class2dSizeHeight = 64 / ratio

    # Histogram graph size
    graphSizeHeight = 60
    graphSizeWidth = 75

    rowSpaceHeight = 62
    rowLineheight = 4

    # column x position
    col0 = 10
    col1 = 75
    col2 = 140

    pdf.cell(col0width, colH, 'Randomly selected acqusition images and particle picks:', 0, 2, 'L')

    startx = pdf.get_x()
    starty = pdf.get_y()+colH

    # Micrograph display
    # row1 col0
    if os.path.exists(main.plot_dir+'/'+xml_session.sessionName+'_ptcls_per_mic_normalised.png'):
        pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_ptcls_per_mic_normalised.png', x = 20, y = starty+(rowSpaceHeight*0), w = graphSizeWidth, h = graphSizeHeight, type = '', link = '')
    #if os.path.exists(main.plot_dir+'/particles/'+search_mics.random0JPG):
    #    pdf.image(main.plot_dir+'/particles/'+search_mics.random0JPG, x = col0, y = starty+(62*0), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
    # row1 col1
    if os.path.exists(main.plot_dir+'/'+xml_session.sessionName+'_ptcls_per_mic_clustered.png'):
        pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_ptcls_per_mic_clustered.png', x = 105, y = starty+(rowSpaceHeight*0), w = graphSizeWidth, h = graphSizeHeight, type = '', link = '')
    #if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random0Name)+'_nnclustering_particles.png'):
    #    pdf.image(main.plot_dir+'/particles/'+str(search_mics.random0Name)+'_nnclustering_particles.png', x = col1, y = starty+(rowSpaceHeight*0), w = pickSizeWidth, h = pickSizeHeight, type = '', link = '')

    # row2 col0
    if os.path.exists(main.plot_dir+'/particles/'+search_mics.random1JPG):
        pdf.image(main.plot_dir+'/particles/'+search_mics.random1JPG, x = col0, y = starty+(rowSpaceHeight*1), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
    # row2 col1
    if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random1Name)+'_nnclustering_particles.png'):
        pdf.image(main.plot_dir+'/particles/'+str(search_mics.random1Name)+'_nnclustering_particles.png', x = col1, y = starty+(62*1), w = pickSizeWidth, h = pickSizeHeight, type = '', link = '')
    # row2 col2
    if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random1Name)+'_2Dclass_res.png'):
        pdf.image(main.plot_dir+'/particles/'+str(search_mics.random1Name)+'_2Dclass_res.png', x = col2, y = starty+(62*1), w = class2dSizeWidth, h = class2dSizeHeight, type = '', link = '')

    # row3 col0
    if os.path.exists(main.plot_dir+'/particles/'+search_mics.random2JPG):
        pdf.image(main.plot_dir+'/particles/'+search_mics.random2JPG, x = col0, y = starty+(rowSpaceHeight*2), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
    # row3 col1
    if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random2Name)+'_nnclustering_particles.png'):
        pdf.image(main.plot_dir+'/particles/'+str(search_mics.random2Name)+'_nnclustering_particles.png', x = col1, y = starty+(62*2), w = pickSizeWidth, h = pickSizeHeight, type = '', link = '')
    # row3 col2
    if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random2Name)+'_2Dclass_res.png'):
        pdf.image(main.plot_dir+'/particles/'+str(search_mics.random2Name)+'_2Dclass_res.png', x = col2, y = starty+(62*2), w = class2dSizeWidth, h = class2dSizeHeight, type = '', link = '')

    # row4 col0
    if os.path.exists(main.plot_dir+'/particles/'+search_mics.random3JPG):
        pdf.image(main.plot_dir+'/particles/'+search_mics.random3JPG, x = col0, y = starty+(rowSpaceHeight*3), w = micSizeWidth, h = micSizeHeight, type = '', link = '')
    # row4 col1
    if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random3Name)+'_nnclustering_particles.png'):
        pdf.image(main.plot_dir+'/particles/'+str(search_mics.random3Name)+'_nnclustering_particles.png', x = col1, y = starty+(62*3), w = pickSizeWidth, h = pickSizeHeight, type = '', link = '')
    # row4 col2
    if os.path.exists(main.plot_dir+'/particles/'+str(search_mics.random3Name)+'_2Dclass_res.png'):
        pdf.image(main.plot_dir+'/particles/'+str(search_mics.random3Name)+'_2Dclass_res.png', x = col2, y = starty+(62*3), w = class2dSizeWidth, h = class2dSizeHeight, type = '', link = '')

    # Add labels over the top
    #pdf.set_font('arial','', 6)
    #pdf.cell(0, graphSizeHeight, '', 0, 2, 'L') # Table space
    #pdf.cell(micSizeWidth, rowLineheight, str('test'), 0, 2, 'L')

    #pdf.cell(0, micSizeHeight, '', 0, 2, 'L') # Table space
    #pdf.cell(micSizeWidth, rowLineheight, str(search_mics.random1Name), 0, 2, 'L')

    #pdf.cell(0, micSizeHeight, '', 0, 2, 'L') # Table space
    #pdf.cell(micSizeWidth, rowLineheight, str(search_mics.random2Name), 0, 2, 'L')

    #pdf.cell(0, micSizeHeight, '', 0, 2, 'L') # Table space
    #pdf.cell(micSizeWidth, rowLineheight, str(search_mics.random3Name), 0, 2, 'L')

    ############################
    ## Collection statistics
    ############################

    if os.path.exists(main.plot_dir+'/'+xml_session.sessionName+'_motioncorr_timings.png'):
       # New page for template
       pdf.add_page()
       pdf.set_xy(0, 0)
       pdf.set_font('arial','B', 18)
       pdf.cell(60)
       pdf.cell(75, 20, "eBIC Session Report (EPU)", 0, 2, 'C')
       pdf.cell(-35)
       pdf.set_font('arial','', 8)

       # Header
       pdf.set_font('arial','', 8)
       pdf.cell(col0width, colH, 'The following plot shows the motion correction over time:', 0, 0, 'L')
       pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
       pdf.cell(-col0width) # Reset column position

       # Motioncorrection statistics
       pdf.cell(col0width, colH, 'First collected (time):', 0, 0, 'L')
       pdf.cell(0, colH, str(countMotionCor.firstCollectedTime), 0, 2, 'L')
       pdf.cell(-col0width) # Reset column position

       pdf.cell(col0width, colH, 'First corrected (time):', 0, 0, 'L')
       pdf.cell(0, colH, str(countMotionCor.firstCorrectedTime), 0, 2, 'L')
       pdf.cell(-col0width) # Reset column position

       pdf.cell(col0width, colH, 'Final collected (time):', 0, 0, 'L')
       pdf.cell(0, colH, str(countMotionCor.lastCollectedTime), 0, 2, 'L')
       pdf.cell(-col0width) # Reset column position

       pdf.cell(col0width, colH, 'Final corrected (time):', 0, 0, 'L')
       pdf.cell(0, colH, str(countMotionCor.lastCorrectedTime), 0, 2, 'L')
       pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
       pdf.cell(-col0width) # Reset column position

       pdf.cell(col0width, colH, 'The following graph shows the corrected versus collected file timestamps', 0, 2, 'L')
       pdf.cell(col0width, colH, 'This is in alpha testing, take with a pinch of salt', 0, 2, 'L')
       pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
       pdf.cell(-col0width) # Reset column position

       # Insert exposure plot
       pdf.cell(col1width*-2, colH, '', 0, 2, 'L') # Table space
       pdf.image(main.plot_dir+'/'+xml_session.sessionName+'_motioncorr_timings.png', x = 20, y = None, w = 140, h = 80, type = '', link = '')

    # Write out report
    pdf.output(report, 'F')

    # Report
    print('Generated PDF report in '+report)
    print()

###
# Execute python functions here
###

main()
