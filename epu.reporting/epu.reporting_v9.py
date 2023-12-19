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

import os, time
import sys
from pathlib import Path
import argparse
parser = argparse.ArgumentParser()
import glob
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.transforms as mtrans
import seaborn as sns
import shutil
import math
from tqdm import tqdm

from functools import reduce

from fpdf import FPDF

# Suppress python warnings....
import warnings
warnings.filterwarnings("ignore")

# Argument parsing
# https://stackoverflow.com/questions/11604653/how-to-add-command-line-arguments-with-flags-in-python3
# DEV need to include in here option to search through all directoeis by time or one specific directory
parser.add_argument("-i", "--input", nargs='+', help="Required: Input directory", required=True) # Allow multiple inputs
parser.add_argument("-o", "--output", help="Required: Output directory for analysis", required=True)
parser.add_argument("--suffix", help="Optional: Suffix for PDF report")
parser.add_argument("--prefix", help="Optional: Prefix for PDF report")
parser.add_argument("-c", "--csv", help="Optional: Output csv name")
parser.add_argument("-cc", "--corcount", help="Optional: Y = Count motion correction mics and timestamps - slow (default: N)")
parser.add_argument("-s", "--session", help="Optional: Input search term for Visit")
parser.add_argument("--include", help="Optional: Exclude all but rows containing the include term, no spaces, single input i.e. KriosI")
parser.add_argument("--exclude", nargs='+', help="Exclude specific errors from report, multi input allowed (Epu_dm,Atl_dm, Tomo, Super, Error)") # Allow multiple inputs
parser.add_argument("-p", "--period", help="Optional: Directories created within n days will be analysed (where n = time in days)")
parser.add_argument("--test", help="Optional: Print out directories for analysis and exit")
parser.add_argument("-ra", "--reanalyse", help="Optional: Global analysis and reporting on existing data only")
parser.add_argument("-b", "--BAG", help="Optional: Show BAG PI, Institute or Anonymise (B/I/A)")
parser.add_argument("-t", "--timemethod", help="Optional: X = XML interal timestamp (slow/accurate), F = file timestamp (fast/inaccurate), M = MRC header timestamp (default = X)")
parser.add_argument("--screening", help="Optional: Perform per square screening analysis (default = N)")
parser.add_argument("--silent", help="Optional: Report less to the terminal")
parser.add_argument("--sizecalc", help="Optional: Y = Calculate sizes of Supervisor/raw/processed directories")
args = parser.parse_args()
print( "input: {}".format(args.input))
print( "output: {}".format(args.output))
print( "include: {}".format(args.include))
print( "exclude: {}".format(args.exclude))
        #args.session,
        #args.include,
        #args.period,
        #args.output,
        #args.reanalyse,
        #args.test,
        #args.timemethod
        #))

# Naughty global variables
#figHeight =
#figWidth =

def main():

    # Input directory to search for data directories
    main.dirInput = args.input
    # Time in days for analysis
    # Deal with no entry as meaning analyse all
    if not args.period:
        main.periodInput = 0
        main.periodReport = 'All'
    else:
        main.periodInput = args.period
        main.periodReport = args.period

    # BAG name, institute or anonymise
    if args.BAG == 'B':
        main.BAGname = 'BAG'
        main.xplot = 'BAG'
    elif args.BAG == 'I':
        main.BAGname = 'Institute'
        main.xplot = 'Institute'
    elif args.BAG == 'A':
        main.BAGname = 'Anonymous'
        main.xplot = 'BAG (anonymised)'
    else:
        main.BAGname = 'BAG'
        main.xplot = 'BAG'

    # Output report file names
    if not args.csv:
        main.outputName = 'epu_time_period_report'
    else:
        main.outputName = args.csv

    # Out suffix if supplied
    if not args.suffix:
        main.outputName = str(main.outputName)
    else:
        main.outputName = str(main.outputName)+'_'+str(args.suffix)

    # Out prefix if supplied
    if not args.prefix:
        main.outputName = str(main.outputName)
    else:
        main.outputName = str(args.prefix)+'_'+str(main.outputName)

    # include suffix handling
    if not args.include:
       args.include = ''

    # Get command used to execute this python script
    python = os.path.basename(sys.executable)
    # Join will concatenate list into single string
    script = ' '.join(sys.argv)
    main.command = str(python)+' '+str(script)

    # DEV need some kind of error checking based on no directories returned
    # This will find the visit directories in the input paths provided
    dirAnalyse = searchDirectories(main.dirInput, main.periodInput)
    main.visitNoAll = len(dirAnalyse)

    # If test dry run only print directories to analyse
    if args.test == 'Y':
        print(dirAnalyse)
    elif args.reanalyse == 'Y':
        # Gather data and analyse csv, write global reports
        globalAnalysis()  
    else:
        # Run analysis
        runAnalysis(dirAnalyse)
        # Gather data and analyse csv, write global reports
        globalAnalysis()

    print('Finished')

def globalAnalysis():
    csv_dir = Path(args.output+'/EMinsight/csv')
    global_report_dir = Path(args.output+'/global_report')
    global_plot_dir = Path(args.output+'/global_report/plots')
    global_csv_dir = Path(args.output+'/global_report/csv')

    # Only make report directory, if doesn't exist
    if global_plot_dir.exists() and global_plot_dir.is_dir():
        print('Output directory exists')
        print()
        #shutil.rmtree(global_report_dir)
        #os.mkdir(global_report_dir)
    else:
        # https://stackoverflow.com/questions/5210778/elegant-way-to-make-all-dirs-in-a-path
        os.makedirs(global_plot_dir)
        os.makedirs(global_csv_dir)

    # Combine data in csv files into a single table for comparison and analysis
    # Sanity check whether there are csv and then combine them
    if glob.glob(str(csv_dir)+'/*_optics.csv'):
        print('Found individual session *_optics.csv, combining for global reports.')
        combineCsv(str(global_csv_dir), str(csv_dir), str('optics'))
    else:
        print('No *_optics.csv found!')

    if glob.glob(str(csv_dir)+'/*_session.csv'):
        print('Found individual session *_session.csv, combining for global reports.')
        combineCsv(str(global_csv_dir), str(csv_dir), str('session'))
        # Count rows
        main.sessionNo = countCsvRow(str(global_csv_dir), str(csv_dir), str('session'))
        # Count errors
        errors = countCsvError(str(global_csv_dir), str(csv_dir), str('session'))
        main.sessionErrorNone = float(errors[0])
        main.sessionErrorSuper = float(errors[1])
        main.sessionErrorAtlas = float(errors[2])
        main.sessionErrorEpu = float(errors[3])
        main.sessionErrorTomo = float(errors[4])
        main.sessionErrorEmpty = float(errors[5])
        # Calculate session error totals for reporting
        main.sessionErrorAll = main.sessionErrorSuper + main.sessionErrorAtlas + main.sessionErrorEpu
        # Count processed
        main.processedNo = countCsvRowTerm(str(global_csv_dir), str(csv_dir), str('session'), str('Processed'), str('Y'))
    else:
        print('No *_session.csv found!')

    if glob.glob(str(csv_dir)+'/*_processed.csv'):
        print('Found individual session *_processed.csv, combining for global reports.')
        combineCsv(str(global_csv_dir), str(csv_dir), str('processed'))
        #main.processedNo = countCsvRow(str(global_csv_dir), str(csv_dir), str('processed'))
    else:
        print('No *_processed.csv found!')

    if glob.glob(str(csv_dir)+'/*_shots.csv'):
        print('Found individual session *_shots.csv, combining for global reports.')
        combineCsv(str(global_csv_dir), str(csv_dir), str('shots'))
        #main.processedNo = countCsvRow(str(global_csv_dir), str(csv_dir), str('processed'))
    else:
        print('No *_shots.csv found!')

    if glob.glob(str(csv_dir)+'/*_advisory.csv'):
        print('Found individual session *_advisory.csv, combining for global reports.')
        combineCsv(str(global_csv_dir), str(csv_dir), str('advisory'))
        #main.processedNo = countCsvRow(str(global_csv_dir), str(csv_dir), str('processed'))
    else:
        print('No *_advisory.csv found!')

    # DEV DEV DEV
    # Need to combine shots.csv so we can find sessions with interesting shot analyses
    #if glob.glob(str(csv_dir)+'/*_shots.csv'):
    #    print('Found individual session *_shots.csv, combining for global reports.')
    #    combineCsv(str(global_csv_dir), str(csv_dir), str('shots'))
        #main.processedNo = countCsvRow(str(global_csv_dir), str(csv_dir), str('processed'))
    #else:
    #    print('No *_shots.csv found!')

    print('Creating global session statistic plots for reports:')

    # Individual session bar charts to show session timings - superceeded by histogram plotting
    plotTimings(str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv', global_plot_dir)
    # Make a bar chart showing mics per session/day
    plotMicCount(str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv', global_plot_dir)
    plotMicArea(str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv', global_plot_dir)
    plotAreaRate(str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv', global_plot_dir)
    # Make a scatter plot showing mag - data rate correlation
    plotAreaRateAnalysis(str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv', global_plot_dir)
    # Plot optics set up
    plotOptics(str(global_csv_dir)+'/'+str(main.outputName)+'_optics.csv', global_plot_dir)
    # Plot sessions set up
    plotSessions(str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv', global_plot_dir)
    # Plot processed set up
    plotProcessed(str(global_csv_dir)+'/'+str(main.outputName)+'_processed.csv', global_plot_dir)

    # Make histograms plots for reports, these functions are seaborn and interfere with bar charts in prior section so place after....
    sessionCsv = str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv'
    # Global statistic histogram - session length - metadata indepdendent
    #plotTimingsHistogram(sessionCsv, global_plot_dir)
    plotSessionHistograms(sessionCsv, "Collection time (hrs)", 40, 'gold', 'sienna', 'Collection time (hrs)', 'EPU Session Count', 'Collection time (hrs)', global_plot_dir, 'global_session_time_collection_histogram')
    plotSessionHistograms(sessionCsv, "Setup time (hrs)", 40, 'c', 'darkslategray', 'Setup time (hrs)', 'EPU Session Count', 'Setup time (hrs)', global_plot_dir, 'global_session_time_setup_histogram')
    # Global statistics on sessions in time period
    plotSessionDates(sessionCsv, "EPU session date", "Session count per time period", 'EPU Session Count', 'Year week number', global_plot_dir, 'global_session_counts_per_year_week')
    # Global statistic histogram - micrograph count
    plotSessionHistograms(sessionCsv, "Total_EPU_mics", 40, 'coral', 'darkred', 'EPU session: Total micrographs', 'EPU Session Count', 'Total micrographs', global_plot_dir, 'global_session_mic_count_histogram')
    # Global statistic histogram - micrograph area
    plotSessionHistograms(sessionCsv, "Total area (um^2)", 40, 'mediumslateblue', 'midnightblue', 'EPU session: Total area (um^2)', 'EPU Session Count', 'Total area (um^2)', global_plot_dir, 'global_session_area_total_histogram')
    # Global statistic histogram - micrograph area
    plotSessionHistograms(sessionCsv, "Rate (um^2/hr)", 40, 'steelblue', 'midnightblue', 'EPU session: Rate (um^2/hr)', 'EPU Session Count', 'Rate (um^2/hr)', global_plot_dir, 'global_session_area_rate_histogram')

    # Make scatter plots based on processed analysis
    input = str(global_csv_dir)+'/'+str(main.outputName)+'_processed.csv'
    if os.path.exists(input):
        plotProcessedAnalysis(input, global_plot_dir)
    # BAG performance
    BAGperformance(str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv', global_plot_dir)
    plotBAGProc(str(global_csv_dir)+'/'+str(main.outputName)+'_session_BAG_report.csv', global_plot_dir)
    # Visit types
    plotVisitType(str(global_csv_dir)+'/'+str(main.outputName)+'_session_visit_type.csv', global_plot_dir)
    print()

    # Print useful stuff on csv files
    print('csv file contents')
    csv = str(global_csv_dir)+'/'+str(main.outputName)+'_session.csv'
    print(csv)
    df = pd.read_csv(csv)
    num_lines = len(df)
    print("Number of lines in the CSV file: ", num_lines)
    print()
    csv = str(global_csv_dir)+'/'+str(main.outputName)+'_processed.csv'
    print(csv)
    df = pd.read_csv(csv)
    num_lines = len(df)
    print("Number of lines in the CSV file: ", num_lines)
    print()

    print('Creating reports:')
    # Create PDF report
    csv = main.outputName+'_processed.csv'
    if os.path.exists(str(global_csv_dir)+'/'+csv):
        construct_global_process_report(str(global_plot_dir), str(global_csv_dir), str(global_report_dir))
    csv = main.outputName+'_session.csv'
    if os.path.exists(str(global_csv_dir)+'/'+csv):  
        construct_global_session_report(str(global_plot_dir), str(global_csv_dir), str(global_report_dir)) 

def combineCsv(path, csv_in, suffix):
    # Delete any previous combined CSV
    combined_csv_path = os.path.join(path, f'{main.outputName}_{suffix}.csv')
    if os.path.exists(combined_csv_path):
        os.remove(combined_csv_path)
        print(f'Deleting old global CSV report: {main.outputName}_{suffix}.csv')
    else:
        print(f'No prior global CSV report found, creating new: {main.outputName}_{suffix}.csv')

    # Initialize an empty DataFrame to store the combined data
    combined_df = None
    excluded_count = 0  # Initialize a count for excluded rows

    # Iterate over all CSV files in the input directory
    csv_files = glob.glob(os.path.join(csv_in, f'*{suffix}.csv'))
    for csv_file in csv_files:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(csv_file)

        # Apply filtering based on args.exclude
        if args.exclude:
            for exclude in args.exclude:
                excluded_rows = df[df.apply(lambda row: row.astype(str).str.contains(exclude, case=True).any(), axis=1)]
                df = df[~df.apply(lambda row: row.astype(str).str.contains(exclude, case=True).any(), axis=1)]
                excluded_count += len(excluded_rows)
        else:
            excluded_count = np.nan

        # Apply filtering based on args.include
        if args.include:
            df = df[df.apply(lambda row: row.astype(str).str.contains(args.include, case=True).any(), axis=1)]
            included_count = len(df)
        else:
            included_count = np.nan

        # If Visit search term is provided, filter the DataFrame
        if args.session:
            df = df[df['Visit'].str.contains(args.session)]

        # Append the DataFrame to the combined DataFrame
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    
    # Calculate the excludeFilterNo variable
    combineCsv.excludeFilterNo = excluded_count
    combineCsv.includeFilterNo = included_count

    # Sort the combined DataFrame by date if "EPU session date" column exists
    if 'EPU session date' in combined_df.columns:
        combined_df.sort_values(by=["EPU session date"], inplace=True)

    # Save the combined DataFrame to a new CSV file
    combined_df.to_csv(combined_csv_path, mode='w', header=True, index=False)

    # Report to the terminal
    print(f'Created new global CSV report: {main.outputName}_{suffix}.csv')
    print()

    ## Totals and runtimes are subject to above filtering
    if suffix == 'session':
        df = pd.read_csv(path+'/'+str(main.outputName)+'_'+suffix+'.csv', parse_dates=["EPU session date"])

        # Get total number of EPU sessions pre-filtering
        row = df.shape[0]
        combineCsv.totalEPU = row

        # Remove duplicates of visits arising due to multiple EPU sessions, to only see visit totals
        df.drop_duplicates('Visit', inplace = True)
        # Group by Visit type bi, nt, nr, cm, etc to know what was scheduled
        dfVisitType = df.groupby(["Visit type"]).size().to_frame('Visits')
        dfVisitType = dfVisitType.reset_index()

        # Get total session lengths per visit type
        dfVisitTime = df.groupby(["Visit type"]).agg(setup = ('Setup time (hrs)','sum'), collect = ('Collection time (hrs)','sum')).reset_index()
        dfVisitTime['runtime'] = dfVisitTime.sum(axis=1)
        dfVisitTime = dfVisitTime.drop(columns=['setup','collect'])
        #print(dfTime)

        # Combine dataframes
        data_frames = [dfVisitType, dfVisitTime]
        dfVisit = reduce(lambda  left,right: pd.merge(left,right,on=['Visit type'],how='outer'), data_frames).fillna('0')

        #print(dfVisitType)
        dfVisit.to_csv(path+'/'+str(main.outputName)+'_'+suffix+'_visit_type.csv', mode='w', header=True, index=False)

    ## Totals and runtimes are subject to above filtering
    if suffix == 'session':
        df = pd.read_csv(path+'/'+str(main.outputName)+'_'+suffix+'.csv', parse_dates=["EPU session date"])

        # Get total number of EPU sessions pre-filtering
        row = df.shape[0]
        combineCsv.totalEPU = row

        # Group by Visit type bi, nt, nr, cm, etc to know what was scheduled
        dfVisitType = df.groupby(["Visit type"]).size().to_frame('Visits')
        dfVisitType = dfVisitType.reset_index()
        print(dfVisitType)

        # Get total session lengths per visit type
        dfVisitTime = df.groupby(["Visit type"]).agg(setup = ('Setup time (hrs)','sum'), collect = ('Collection time (hrs)','sum')).reset_index()
        dfVisitTime['runtime'] = dfVisitTime.sum(axis=1)
        dfVisitTime = dfVisitTime.drop(columns=['setup','collect'])
        print(dfVisitTime)
        
        # Get total session types SPA/Tomo
        #dfSessionType = df.groupby(["Visit type"])['Type'].transform(lambda x: x[x.str.contains('SPA')].count())
        #print(dfSessionType)
        #dfSessionType['SPA'] = dfSessionType.sum(axis=1)

        # Combine dataframes
        data_frames = [dfVisitType, dfVisitTime]
        dfVisit = reduce(lambda  left,right: pd.merge(left,right,on=['Visit type'],how='outer'), data_frames).fillna('0')

        #print(dfVisitType)
        dfVisit.to_csv(path+'/'+str(main.outputName)+'_'+suffix+'_session_type.csv', mode='w', header=True, index=False)

    ## Perform additional separate operations on optics.csv to get optics settings per grid type
    if suffix == 'optics':
        # Read in optics data
        df = pd.read_csv(path+'/'+str(main.outputName)+'_'+suffix+'.csv', parse_dates=["EPU session date"])

        ## Pull out the most common configuration for the presets
        ## MOST COMMON PRESET DATA RANGES per grid type
        # Master dataframe to put all presets into
        dfOpticsAll = pd.DataFrame(columns=['Hole (um)', 'Space (um)'])

        ## Loop to go through all presets
        presetList = ['Atlas_mag', 'Atlas_spot', 'Atlas_C2', 'Atlas_beam', 'Sq_mag', 'Sq_spot', 'Sq_C2', 'Sq_beam', 'Hole_mag', 'Hole_spot', 'Hole_C2', 'Hole_beam', 'Data_mag', 'Data_spot', 'Data_C2', 'Data_beam']
        for preset in presetList:
            # Get min max values
            # Groups by Hole/Space size, then determines the min and max values on the preset currently being queried by loop
            dfMinMax = df.groupby(['Hole (um)','Space (um)'])[preset].aggregate(['min','max'])
            dfMinMax.rename(columns = {'max':preset+'_max'}, inplace = True)
            dfMinMax.rename(columns = {'min':preset+'_min'}, inplace = True)
            
            # Combine dataframes
            #data_frames = [dfMinMax]
            #dfOpticsMag = reduce(lambda  left,right: pd.merge(left,right,on=["Hole (um)","Space (um)"],how='outer'), data_frames).fillna('0')
            dfOpticsMinMax = dfMinMax

            # Add dataframes to master dataframe
            data_frames = [dfOpticsAll, dfOpticsMinMax]
            dfOpticsAll = reduce(lambda  left,right: pd.merge(left,right,on=["Hole (um)","Space (um)"],how='outer'), data_frames).fillna('0') 

        # Save preset usage to csv file
        #print(dfOpticsAll)
        dfOpticsAll.to_csv(path+'/'+str(main.outputName)+'_'+suffix+'_preset_ranges.csv', mode='w', header=True, index=False)

        ## ALTERNATIVE METHOD
        ## ONLY PULL OUT MOST COMMON MODAL VALUES FOR PRESETS per grid type
        ## Pull out the most common configuration for the presets
        # Master dataframe to put all presets into
        dfOpticsAll = pd.DataFrame(columns=['Hole (um)', 'Space (um)'])

        presetList = ['Atlas', 'Sq', 'Hole', 'Data']
        for preset in presetList:
            # Preset magnification - most common
            query = preset+'_mag'
            dfMagMode = df.groupby(['Hole (um)','Space (um)'])[query].agg(lambda x: pd.Series.mode(x)[0])
            # Preset spot size - most common
            query = preset+'_spot'
            dfSpotMode = df.groupby(['Hole (um)','Space (um)'])[query].agg(lambda x: pd.Series.mode(x)[0])
            # Preset C2 size - most common
            query = preset+'_C2'
            dfC2Mode = df.groupby(['Hole (um)','Space (um)'])[query].agg(lambda x: pd.Series.mode(x)[0])
            # Preset probe mode - most common
            query = preset+'_probe'
            dfProbeMode = df.groupby(['Hole (um)','Space (um)'])[query].agg(lambda x: pd.Series.mode(x)[0])
            # Preset defocus - most common
            query = preset+'_DF'
            dfDefocusMode = df.groupby(['Hole (um)','Space (um)'])[query].agg(lambda x: pd.Series.mode(x)[0])
            # Preset beam size - range (include column renaming to make column name unique to preset)
            query = preset+'_beam'
            dfFilt = df.loc[df[query] > 0]
            dfBeamMinMax = dfFilt.groupby(['Hole (um)','Space (um)'])[query].aggregate(['min','max'])
            dfBeamMinMax.rename(columns = {'max':query+'_max'}, inplace = True)
            dfBeamMinMax.rename(columns = {'min':query+'_min'}, inplace = True)

            # Add dataframes to master dataframe
            data_frames = [dfOpticsAll, dfMagMode, dfSpotMode, dfC2Mode, dfProbeMode, dfDefocusMode, dfBeamMinMax]
            dfOpticsAll = reduce(lambda  left,right: pd.merge(left,right,on=["Hole (um)","Space (um)"],how='outer'), data_frames).fillna('0') 

        # Save preset usage to csv file
        #print(dfOpticsAll)
        dfOpticsAll.to_csv(path+'/'+str(main.outputName)+'_'+suffix+'_preset_modes.csv', mode='w', header=True, index=False)

    ## Perform additional separate operations on processed.csv to characterise grid specific outcomes
    if suffix == 'processed':
        df = pd.read_csv(path+'/'+str(main.outputName)+'_'+suffix+'.csv', parse_dates=["EPU session date"])

        ## Pull out the most common configuration for the presets
        ## MOST COMMON PRESET DATA RANGES per grid type
        # Master dataframe to put all presets into
        dfSessionAll = pd.DataFrame(columns=['Hole (um)', 'Space (um)','Shots per hole'])

        paramList = ['Beam (um)','apix','Early_motion_mean', 'Late_motion_mean', 'Total_motion_mean']
        last = paramList[-1]
        for param in paramList:
            # Get min max values
            # Groups by Hole/Space size, then determines the min and max values on the preset currently being queried by loop
            if param == last:
                dfMean = df.groupby(['Hole (um)','Space (um)','Shots per hole'])[param].aggregate(['mean','count']).round(2)
                dfMean.rename(columns = {'mean':param+'_mean'}, inplace = True)
            else:
                dfMean = df.groupby(['Hole (um)','Space (um)','Shots per hole'])[param].aggregate(['mean']).round(2)
                dfMean.rename(columns = {'mean':param+'_mean'}, inplace = True)
            
            # Combine dataframes
            dfMotionMean = dfMean

            # Add dataframes to master dataframe
            data_frames = [dfSessionAll, dfMotionMean]
            dfSessionAll = reduce(lambda  left,right: pd.merge(left,right,on=["Hole (um)","Space (um)","Shots per hole"],how='outer'), data_frames).fillna('0') 

        # Save preset usage to csv file
        #print(dfOpticsAll)
        dfSessionAll.to_csv(path+'/'+str(main.outputName)+'_'+suffix+'_motion.csv', mode='w', header=True, index=False)

    ## Perform additional separate operations on session.csv to characterise and report on BAGs
    if suffix == 'session':
       df = pd.read_csv(path+'/'+str(main.outputName)+'_'+suffix+'.csv', parse_dates=["EPU session date"])
       ## Calculate a normalisation of the area rate (performance per BAG) and put into main session csv
       # As range of values min to max
       normalized_df_umrate=(df['Rate (um^2/hr)']-df['Rate (um^2/hr)'].min())/(df['Rate (um^2/hr)'].max()-df['Rate (um^2/hr)'].min())
       # As range from 0 to max
       df['~rate'] = df['Rate (um^2/hr)'] 
       df['~rate']=(df['Rate (um^2/hr)']-0)/(df['Rate (um^2/hr)'].max()-0)
       df['~rate'] = df['~rate'].astype(float).round(2)
       
       ### The following will set up a new csv as a BAG report
       ## Pull out nromalised data rates for each BAG
       #dfRate = df[['BAG', 'Institute', '~rate']]
       #dfRate = df.groupby(["BAG", "Institute"])['~rate'].mean()s

       ## Get Microscope totals, for visits (note drop of duplicates)
       # Count number of sessions with Krios I
       dfKrios1Visit = df.query("Instrument == 'Krios1'")
       #dfKrios1.drop_duplicates('Visit', inplace = True)
       dfKrios1 = dfKrios1Visit.drop_duplicates(["Visit"]).copy()
       dfKrios1 = dfKrios1[['BAG', 'Institute', 'Anonymous', 'Instrument']].value_counts().to_frame('Krios1').reset_index().drop(columns=['Instrument'])

       dfKrios2Visit = df.query("Instrument == 'Krios2'")
       #dfKrios2.drop_duplicates('Visit', inplace = True)
       dfKrios2 = dfKrios2Visit.drop_duplicates(["Visit"]).copy()
       dfKrios2 = dfKrios2[['BAG', 'Institute', 'Anonymous', 'Instrument']].value_counts().to_frame('Krios2').reset_index().drop(columns=['Instrument'])

       dfKrios3Visit = df.query("Instrument == 'Krios3'")
       #dfKrios3.drop_duplicates('Visit', inplace = True)
       dfKrios3 = dfKrios3Visit.drop_duplicates(["Visit"]).copy()
       dfKrios3 = dfKrios3[['BAG', 'Institute', 'Anonymous', 'Instrument']].value_counts().to_frame('Krios3').reset_index().drop(columns=['Instrument'])

       dfKrios4Visit = df.query("Instrument == 'Krios4'")
       #dfKrios4.drop_duplicates('Visit', inplace = True)
       dfKrios4 = dfKrios4Visit.drop_duplicates(["Visit"]).copy()
       dfKrios4 = dfKrios4[['BAG', 'Institute', 'Anonymous', 'Instrument']].value_counts().to_frame('Krios4').reset_index().drop(columns=['Instrument'])

       dfTalosVisit = df.query("Instrument == 'Talos'")
       #dfTalos.drop_duplicates('Visit', inplace = True)
       dfTalos = dfTalosVisit.drop_duplicates(["Visit"]).copy()
       dfTalos = dfTalos[['BAG', 'Institute', 'Anonymous', 'Instrument']].value_counts().to_frame('Talos').reset_index().drop(columns=['Instrument'])

       dfGlaciosVisit = df.query("Instrument == 'GlcsII'")
       #dfGlacios.drop_duplicates('Visit', inplace = True)
       dfGlacios = dfGlaciosVisit.drop_duplicates(["Visit"]).copy()
       dfGlacios = dfGlacios[['BAG', 'Institute', 'Anonymous', 'Instrument']].value_counts().to_frame('GlaciosII').reset_index().drop(columns=['Instrument'])

       # Combine dataframes
       data_frames = [dfKrios1, dfKrios2, dfKrios3, dfKrios4, dfTalos, dfGlacios]
       dfScope = reduce(lambda  left,right: pd.merge(left,right,on=['BAG', 'Institute', 'Anonymous'],how='outer'), data_frames).fillna('0')

       # merge removes float
       dfScope['Krios1'] = pd.to_numeric(dfScope['Krios1'], errors='coerce')
       dfScope['Krios2'] = pd.to_numeric(dfScope['Krios2'], errors='coerce')
       dfScope['Krios3'] = pd.to_numeric(dfScope['Krios3'], errors='coerce')
       dfScope['Krios4'] = pd.to_numeric(dfScope['Krios4'], errors='coerce')
       dfScope['Talos'] = pd.to_numeric(dfScope['Talos'], errors='coerce')
       dfScope['GlaciosII'] = pd.to_numeric(dfScope['GlaciosII'], errors='coerce')

       # Apply normalisations 
       column_list = ['Krios1', 'Krios2', 'Krios3', 'Krios4', 'Talos', 'GlaciosII']
       #dfScope['Total'] = dfScope.sum(axis=1) 
       dfScope['Total'] = dfScope.loc[:,column_list].sum(axis=1)
       dfScope['Krios1'] = dfScope.apply(lambda row : roundup(row['Krios1']/row['Total'],5), axis = 1)
       dfScope['Krios2'] = dfScope.apply(lambda row : roundup(row['Krios2']/row['Total'],5), axis = 1)
       dfScope['Krios3'] = dfScope.apply(lambda row : roundup(row['Krios3']/row['Total'],5), axis = 1)
       dfScope['Krios4'] = dfScope.apply(lambda row : roundup(row['Krios4']/row['Total'],5), axis = 1)
       dfScope['Talos'] = dfScope.apply(lambda row : roundup(row['Talos']/row['Total'],5), axis = 1)
       dfScope['GlaciosII'] = dfScope.apply(lambda row : roundup(row['GlaciosII']/row['Total'],5), axis = 1)

       # Drop Total now the normalisation is done
       dfScope.drop(columns=['Total'])
       #print(dfScope)
 
       ## Get accumulated totals for processed Y/N per BAG
       # Count number of sessions with Processed
       dfProcY = df[['BAG', 'Institute', 'Anonymous', 'Processed']].replace('N', np.nan).value_counts().to_frame('ProcY').reset_index()
       # Count number of sessions without Processed
       dfProcN = df[['BAG', 'Institute', 'Anonymous', 'Processed']].replace('Y', np.nan).value_counts().to_frame('ProcN').reset_index()
       # Merge into single data frame, drop original processed column
       dfProc = pd.merge(left=dfProcY.drop('Processed', axis=1), right=dfProcN.drop('Processed', axis=1), on=['BAG', 'Institute', 'Anonymous'], how='outer')
       #print(dfProc)

       ## Get total session lengths
       dfTime = df.groupby(["BAG", "Institute", "Anonymous"]).agg(setup = ('Setup time (hrs)','sum'), collect = ('Collection time (hrs)','sum')).reset_index()
       dfTime['runtime'] = dfTime.sum(axis=1)
       dfTime = dfTime.drop(columns=['setup','collect'])
       #print(dfTime)

       ## Get number of unique visit, sessions and session type per BAG
       # Number of sessions
       dfBAGsession = df.groupby(["BAG", "Institute", "Anonymous"]).size().to_frame('Sessions')
       #print(dfBAGsession)

       # Types of sessions
       dfSPA = df.query("Type == 'SPA'")
       #dfSPA.drop_duplicates('Visit', inplace = True)
       dfSPA = dfSPA[['BAG', 'Institute', 'Anonymous', 'Type']].value_counts().to_frame('SPA').reset_index().drop(columns=['Type'])
       dfTomo = df.query("Type == 'Tomo'")
       #dfTomo.drop_duplicates('Visit', inplace = True)
       dfTomo = dfTomo[['BAG', 'Institute', 'Anonymous', 'Type']].value_counts().to_frame('Tomo').reset_index().drop(columns=['Type'])
       dfUnknown = df.query("Type == 'Unknown'")
       #dfUnknown.drop_duplicates('Visit', inplace = True)
       dfUnknown = dfUnknown[['BAG', 'Institute', 'Anonymous', 'Type']].value_counts().to_frame('Unknown').reset_index().drop(columns=['Type'])
       # Combine dataframes
       data_frames = [dfSPA, dfTomo, dfUnknown]
       dfSessionType = reduce(lambda  left,right: pd.merge(left,right,on=['BAG', 'Institute', 'Anonymous'],how='outer'), data_frames).fillna('0')
       #print(dfSessionType)

       ## Get the number of visits per BAG
       # First group to remove multiple EPU session per visit
       dfBAGvisit = df.groupby(["Visit", "BAG", "Institute", "Anonymous"]).size().to_frame('Visits')
       # Then group to get number of visits
       dfBAGvisits = dfBAGvisit.groupby(["BAG", "Institute", "Anonymous"]).size().to_frame('Visits')
       combineCsv.BAGvisits = dfBAGvisits['Visits'].sum()
       #print(dfBAGvisits)

       ## Count errors in sessions
       dfError = df.value_counts(subset=['BAG', 'Institute', 'Anonymous', 'Error']).to_frame('Errors').reset_index()
       dfError = dfError[dfError["Error"].str.contains("None|EPU_dm|Atl_dm|Tomo") == False]
       dfErrors = dfError.groupby(["BAG", "Institute", "Anonymous"])['Errors'].sum().to_frame('Errors')

       ## Merge dataframes for BAG report csv
       #data_frames = [dfBAG, dfRate, dfErrors, dfProc]
       data_frames = [dfBAGvisits, dfBAGsession, dfErrors, dfTime, dfSessionType, dfScope, dfProc]
       dfFinal = reduce(lambda  left,right: pd.merge(left,right,on=['BAG', 'Institute', 'Anonymous'],how='outer'), data_frames).fillna('0')

       # Infer sessions that you don't known processing state
       dfFinal['ProcUnknown'] = dfFinal.apply(lambda row : (float(row['Sessions'])-(float(row['ProcN'])+float(row['ProcY']))), axis = 1)

       # Apply some normalisations
       dfFinal['ProcY'] = dfFinal.apply(lambda row : float(row['ProcY'])/(float(row['Sessions'])), axis = 1)
       dfFinal['ProcN'] = dfFinal.apply(lambda row : float(row['ProcN'])/(float(row['Sessions'])), axis = 1)
       dfFinal['ProcUnknown'] = dfFinal.apply(lambda row : float(row['ProcUnknown'])/(float(row['Sessions'])), axis = 1)

       # Apply normalisations
       dfFinal['SPA'] = dfFinal.apply(lambda row : float(row['SPA'])/(float(row['Sessions'])), axis = 1)
       dfFinal['Tomo'] = dfFinal.apply(lambda row : float(row['Tomo'])/(float(row['Sessions'])), axis = 1)
       dfFinal['Unknown'] = dfFinal.apply(lambda row : float(row['Unknown'])/(float(row['Sessions'])), axis = 1)

       # Apple sorting based on what you are plotting
       column = main.BAGname
       dfFinal.sort_values(by=[column], inplace=True)

       # Save merged dataframe to csv
       #print(dfFinal)
       dfFinal.to_csv(path+'/'+str(main.outputName)+'_'+suffix+'_BAG_report.csv', mode='w', header=True, index=False)

def countCsvRow(path, csv_in, suffix):
    # Get number of rows in final global csv
    df = pd.read_csv(path+'/'+str(main.outputName)+'_'+suffix+'.csv')
    return len(df.index)

def countCsvRowTerm(path, csv_in, suffix, column, term):
    # Get number of rows in final global csv
    df = pd.read_csv(path+'/'+str(main.outputName)+'_'+suffix+'.csv')
    try:
       count = df[column].value_counts()[term]
    except:
       count = '0'
    return count

def countCsvError(path, csv_in, suffix):
    # Get number of rows in final global csv
    df = pd.read_csv(path+'/'+str(main.outputName)+'_'+suffix+'.csv')
    try:
       none = df['Error'].value_counts()['None']
    except:
       none = '0'
    try:
       supervisor = df['Error'].value_counts()['Super']
    except:
       supervisor = '0'
    try:
       atlas = df['Error'].value_counts()['Atl_dm']
    except:
       atlas = '0'
    try:
       epu = df['Error'].value_counts()['EPU_dm']
    except:
       epu = '0'
    try:
       tomo = df['Error'].value_counts()['Tomo']
    except:
       tomo = '0'
    try:
       empty = df['Error'].value_counts()['Empty']
    except:
       empty = '0'

    return none, supervisor, atlas, epu, tomo, empty

def runAnalysis(dirAnalyse):
    j = 1
    k = 0

    length = main.visitNoAll

    # Loop through each session folder and then look inside each Supervisor directory, if an EPU directory then process it
    # Each session or visit is in dirAnalyse as a list
    # The epu.xml_summary.py will then loop through supervisor directories in the current session or visit being analysed
    for i in dirAnalyse:
        print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        print('Directory '+i)
        print('Working on '+str(j)+'/'+str(length)+': '+i)
        print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        j = j+1
        # Search supervisor returns a dataframe with column names [Directory, Epu, Atlas] and thus what their type is
        #print(searchSupervisor(i))
        #dirSupervisor = searchSupervisor(i)

        ## If epu or atlas files not found, something wrong or it's tomo, iterate to next directory

        ## NOW you can have a way to look through the dataframe, if epu session then use it, if not then skip
        # Remember you'll need code to skip if no atlas found
        # You can then call an error function to write out a csv describing if a directory is missing metadata
        #print(dirSupervisor.at[k, 'Epu'])

        # No need then to have skipping at the visit directory level, want skipping at Supervisor directory level in the epu.xml_summary.py script
        #if not searchSupervisorEPU(i):
            #print('No epu session xml file found, skipping')
            #continue
        #if not searchSupervisorAtlas(i):
            #print('No atlas session xml file found, skipping')
            #continue

        # Characterise motion correction timings, it is slow, default N = no
        if args.corcount == 'Y':
           cc = 'Y'
        elif args.corcount == 'N':
           cc = 'N'
        else:
           cc = 'N'

        # Characterise directory sizes
        if args.sizecalc == 'Y':
           sc = 'Y'
        elif args.sizecalc == 'N':
           sc = 'N'
        else:
           sc = 'N'

        # Which time method to use, xml or file timestamps, default X = xml
        if args.timemethod == 'X':
           time = 'X'
        elif args.timemethod == 'F':
           time = 'F'
        else:
           time = 'X'

        # Silent output, default Y
        if args.silent == 'Y':
           silent = 'Y'
        elif not args.silent:
           silent = 'Y'
        elif args.silent == 'N':
           silent = 'N'
        else:
           silent = 'Y'

        # Square screening report
        if not args.screening:
            args.screening = 'N'
        elif args.screening == 'N':
            args.screening = 'N'
        else:
            args.screening = 'Y'

        # Run epu.xml_summary.py to parse out session statistics
        subprocess.call("epu.xml_summary.py -i "+i+" -r Y -t "+time+" -r N -cc "+cc+" --sizecalc "+sc+" -o "+args.output+" --silent "+silent+" --screening "+args.screening, shell=True)

def searchDirectories(dirInput, timeInput):
    # Get date sorted list of directories
    allDir = []
    # Loop through the inputs directories provided
    for i in dirInput:
        # Get the sub directories of each input directory
        # Note disks may be offline or bad path specified, so skip if not found
        try:
            subDir = (next(os.walk(i))[1])
        except StopIteration:
            pass
        # Loop through sub directories and add to the allDir list
        # Same error checking as above is paths not found
        try:
            subDir
        except NameError:
            pass
        else:
            for j in subDir:
                # Make sure the parent path is prefixed
                allDir.append(i+'/'+j)

    # Sort directories by modification time in ascending order
    dirSearch = allDir
    dirSearch.sort(key=lambda x: os.path.getmtime(x)) # getmtime for modified date, getctime for 'created' date

    #print('Directories to analyse: ')
    #print(dirSearch)

    # Loop through directories and check how old
    # https://stackoverflow.com/questions/12485666/python-deleting-all-files-in-a-folder-older-than-x-days
    now = time.time()

    dirAnalyse = []
    # If no period entered then work on all directories in input directory
    if timeInput == 0:
        for f in dirSearch:
            dirAnalyse.append(f)
    else:
        for f in dirSearch:
            filestamp = os.path.getmtime(f)
            filecompare = now - float(timeInput) * 86400
            if  filestamp > filecompare:
                dirAnalyse.append(f)
                #print(f)

    # Sort directories by modification time in descending order
    dirAnalyse = sorted(dirAnalyse, key=lambda t: -os.stat(t).st_mtime)
    print(dirAnalyse[0])

    # If search term provided then filter the directory list at the point of running analysis
    if not args.session:
        directories = dirAnalyse
    else:
        directories = [k for k in dirAnalyse if str(args.session) in k]

    searchDirectories.filteredNo = len(directories)

    return directories

def searchSupervisor(path):
    # Get list of Supervisor directories
    supervisorList = glob.glob(path+'/Supervisor*')

    DirectoryList = pd.DataFrame({"Directory":[],"Epu":[],"Atlas":[]})
    #epuDirectoryList.append(path) # This could be a way to get the Visit
    # Add dir to directorylist if it contains EpuSession.dm files
    for f in supervisorList:

        # Is EPU directory
        epuSession = glob.glob(f+'/EpuSession.dm')
        if epuSession:
            #print(epuSession)
            Epu = 'EpuSession.dm'
        else:
            Epu = 'None'

        # Is Atlas directory
        epuSession = glob.glob(f+'/ScreeningSession.dm')
        if epuSession:
            #print(epuSession)
            Atlas = 'ScreeningSession.dm'
        else:
            Atlas = 'None'

        # Setup table describing Supervisor directories
        #DirectoryList = DirectoryList.append({'Directory':f,'Epu':Epu,'Atlas':Atlas}, ignore_index=True)
        new_row = pd.DataFrame({'Directory':f,'Epu':Epu,'Atlas':Atlas}, index=[0])
        DirectoryList = pd.concat([new_row,DirectoryList.loc[:]]).reset_index(drop=True)

    # Return list of Supervisor directories that contain EpuSession.dm xmls
    return DirectoryList

def searchSupervisorEPU(path):
    # Get list of Supervisor directories
    supervisorList = glob.glob(path+'/Supervisor*')

    epuDirectoryList = []
    #epuDirectoryList.append(path) # This could be a way to get the Visit
    # Add dir to directorylist if it contains EpuSession.dm files
    for f in supervisorList:
        #print(f)
        epuSession = glob.glob(f+'/EpuSession.dm')
        if epuSession:
            #print(epuSession)
            epuDirectoryList.append(epuSession[0])

    # Return list of Supervisor directories that contain EpuSession.dm xmls
    return epuDirectoryList

def searchSupervisorAtlas(path):
    # Get list of Supervisor directories
    supervisorList = glob.glob(path+'/Supervisor*')

    atlasDirectoryList = []
    #atlasDirectoryList.append(path) # This could be a way to get the Visit
    # Add dir to directorylist if it contains ScreeningSession.dm files
    for f in supervisorList:
        atlasSession = glob.glob(f+'/ScreeningSession.dm')
        #print(atlasSession)
        if atlasSession:
            #print(atlasSession)
            atlasDirectoryList.append(atlasSession[0])

    # Return list of Supervisor directories that contain .dm xmls
    return atlasDirectoryList

def roundup(n, decimals=0):
    # https://realpython.com/python-rounding/
    multiplier = 10 ** decimals
    return math.ceil(n * multiplier) / multiplier

def plotTimings(f, output):
    print('plot timings')

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    df = pd.read_csv(f, parse_dates=['EPU session date','Session date'])
    fig0 = plt.figure(0)
    # https://medium.com/@jb.ranchana/easy-way-to-create-stacked-bar-graphs-from-dataframe-19cc97c86fe3

    figHeight = 3
    figWidth = 12
    plotTimings.figRatio = figWidth / figHeight

    ax = df.plot.bar(x='EPU session date', y=["Setup time (hrs)","Collection time (hrs)"], stacked=True, color=['lightblue','orange'], title='Session Timings', figsize=(figWidth,figHeight))
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d")

    plt.axhline(y=24,linewidth=1,color='grey',linestyle='dashed')
    plt.axhline(y=48,linewidth=1,color='grey',linestyle='dashed')
    plt.axhline(y=72,linewidth=1,color='grey',linestyle='dashed')
    plt.axhline(y=96,linewidth=1,color='grey',linestyle='dashed')
    plt.title("Session timings")
    plt.xlabel("Date")
    plt.ylabel("Hours")
    ax.set_xticklabels(labels=df['EPU session date'], rotation=70, rotation_mode="anchor", ha="right", fontsize=6)

    plt.tight_layout()
    fig0 = ax.get_figure()
    fig0.savefig(str(output)+'/session_timings.png', dpi=300)
    plt.figure(0).clear()
    #plt.show()

def plotTimingsHistogram(f, output):
    df = pd.read_csv(f, parse_dates=['EPU session date','Session date'])
    
    # Plot as histogram
    col = "Collection time (hrs)"
    filename = 'session_timings_collection_histogram'
    yaxis = 'Session count'
    xaxis = 'Collection time (hrs)'
    title = 'Collection time by EPU session'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='orange',
            edgecolor="brown", linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(output)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)
    plt.close(11)

    # Plot as histogram
    col = "Setup time (hrs)"
    filename = 'session_timings_setup_histogram'
    yaxis = 'Session count'
    xaxis = 'Setup time (hrs)'
    title = 'Setup time by EPU session'
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='lightblue',
            edgecolor="navy", linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(output)+'/'+str(filename)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)
    plt.close(11)

def plotVisitType(f, output):
    print('Visit type')
    df = pd.read_csv(f)

    # Make some changes to dataframe so x axis catagories make sense to reader
    df['Visit type'] = df['Visit type'].str.replace('bi','User (bi)')
    df['Visit type'] = df['Visit type'].str.replace('nt','in-house (nt)')
    df['Visit type'] = df['Visit type'].str.replace('nr','in-house (nr)')
    df['Visit type'] = df['Visit type'].str.replace('cm','Comission (cm)')
    # The following are oddities, so rename Other
    df['Visit type'] = df['Visit type'].str.replace('Su','Other')
    df['Visit type'] = df['Visit type'].str.replace('ra','Other')
    df['Visit type'] = df['Visit type'].str.replace('sm','Other')
    df['Visit type'] = df['Visit type'].str.replace('sw','Other')

    # Do some grouping to sum Other's together, re-index to deal with header issue
    #df2 = df.groupby(['Visits'], as_index=False).sum().to_frame('Visits')
    df = df.groupby(["Visit type"]).agg(Visits = ('Visits','sum'), runtime = ('runtime','sum')).reset_index()
    df = df.reset_index().sort_values(['Visits'],ascending=False)

    fig0 = plt.figure(0)
    # https://medium.com/@jb.ranchana/easy-way-to-create-stacked-bar-graphs-from-dataframe-19cc97c86fe3

    figHeight = 4
    figWidth = 5

    ## Number of visits by type
    ax = df.plot.bar(x='Visit type', y=["Visits"], color=['pink'], label=['Visits'], figsize=(figWidth,figHeight))
    # line
    plt.grid(visible='True', which='both', axis='y', color='lightgrey', linestyle='--', linewidth=0.5)
    #plt.minorticks_on()
    # x-axis
    ax.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax.get_xticklabels():
       t.set_transform(t.get_transform()+trans)
    # y-axis
    ax.yaxis.set_label_position("left")
    ax.yaxis.tick_left()
    # Legend
    # Legend
    plt.legend(bbox_to_anchor=(0, 1.15, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, frameon=False, framealpha=1.0)
    # Title
    plt.title("Total visit types")
    plt.xlabel("Visit type")
    plt.ylabel("Count")
    # Save fig
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig0 = ax.get_figure()
    fig0.savefig(str(output)+'/visit_type_counts.png', dpi=300)
    plt.figure(0).clear()
    #plt.show()

    ## Runtime of visits by type
    ax = df.plot.bar(x='Visit type', y=["runtime"], color=['purple'], label=['Runtime'], figsize=(figWidth,figHeight))
    # line
    plt.grid(visible='True', which='both', axis='y', color='lightgrey', linestyle='--', linewidth=0.5)
    #plt.minorticks_on()
    # x-axis
    ax.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax.get_xticklabels():
       t.set_transform(t.get_transform()+trans)
    # y-axis
    ax.yaxis.set_label_position("left")
    ax.yaxis.tick_left()
    # Legend
    # Legend
    plt.legend(bbox_to_anchor=(0, 1.15, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, frameon=False, framealpha=1.0)
    # Title
    plt.title("Total visit times")
    plt.xlabel("Visit type")
    plt.ylabel("Time (hrs)")
    # Save fig
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig0 = ax.get_figure()
    fig0.savefig(str(output)+'/visit_type_runtime.png', dpi=300)
    plt.figure(0).clear()
    #plt.show()

def plotBAGProc(f, output):
    print('Session performance plots: processed statistics')
    df = pd.read_csv(f)
    fig0 = plt.figure(0)
    # https://medium.com/@jb.ranchana/easy-way-to-create-stacked-bar-graphs-from-dataframe-19cc97c86fe3

    figHeight = 4
    figWidth = 5

    #xplot = 'BAG'
    #xplot = 'Institute'
    #xplot = 'Anonymous'
    xplot = main.BAGname

    ## Pipeline initiation
    ax = df.plot.bar(x=xplot, y=["ProcY", "ProcN", "ProcUnknown"], stacked=True, color=['lightblue','orange','lightgrey'], edgecolor = "silver", label=['Yes', 'No', 'Unknown'], figsize=(figWidth,figHeight))
    # x-axis
    ax.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax.get_xticklabels():
       t.set_transform(t.get_transform()+trans)
    # y-axis
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    # Legend
    plt.legend(bbox_to_anchor=(0, 1.15, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, frameon=False, framealpha=1.0)
    # Title
    plt.title("Pipeline started for EPU sessions")
    plt.xlabel(main.xplot)
    plt.ylabel("Percentage")
    # Save fig
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig0 = ax.get_figure()
    fig0.savefig(str(output)+'/BAG_processed_counts.png', dpi=300)
    plt.figure(0).clear()

    ## Number of visits/sessions
    ax = df.plot.bar(x=xplot, y=["Visits", "Sessions"], color=['pink', 'lightblue'], label=['Visits', 'EPU Sessions'], figsize=(figWidth,figHeight))
    # line
    plt.grid(visible='True', which='both', axis='y', color='lightgrey', linestyle='--', linewidth=0.5)
    #plt.minorticks_on()
    # x-axis
    ax.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax.get_xticklabels():
       t.set_transform(t.get_transform()+trans)
    # y-axis
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    # Legend
    # Legend
    plt.legend(bbox_to_anchor=(0, 1.15, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, frameon=False, framealpha=1.0)
    # Title
    plt.title("Total EPU Sessions initiated during Visits")
    plt.xlabel("BAG")
    plt.ylabel("Count")
    # Save fig
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig0 = ax.get_figure()
    fig0.savefig(str(output)+'/BAG_session_counts.png', dpi=300)
    plt.figure(0).clear()
    #plt.show()

    ## Session types
    ax = df.plot.bar(x=xplot, y=["SPA", "Tomo", "Unknown"], stacked=True, color=['lightblue','lightgreen','lightgrey'], edgecolor = "silver", label=['SPA','Tomo','Unknown'], figsize=(figWidth,figHeight))
    # x-axis
    ax.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax.get_xticklabels():
       t.set_transform(t.get_transform()+trans)
    # y-axis
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    # Legend
    plt.legend(bbox_to_anchor=(0, 1.15, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, frameon=False, framealpha=1.0)
    # Title
    plt.title("Session Type")
    plt.xlabel("BAG")
    plt.ylabel("Percentage")
    # Save fig
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig0 = ax.get_figure()
    fig0.savefig(str(output)+'/BAG_session_type.png', dpi=300)
    plt.figure(0).clear()

    ## Runtime of visits
    ax = df.plot.bar(x=xplot, y=["runtime"], color=['purple'], label=['Runtime'], figsize=(figWidth,figHeight))
    # line
    plt.grid(visible='True', which='both', axis='y', color='lightgrey', linestyle='--', linewidth=0.5)
    #plt.minorticks_on()
    # x-axis
    ax.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax.get_xticklabels():
       t.set_transform(t.get_transform()+trans)
    # y-axis
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    # Legend
    # Legend
    plt.legend(bbox_to_anchor=(0, 1.15, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, frameon=False, framealpha=1.0)
    # Title
    plt.title("Total runtime of EPU Sessions during Visits")
    plt.xlabel("BAG")
    plt.ylabel("Hours (hrs)")
    # Save fig
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig0 = ax.get_figure()
    fig0.savefig(str(output)+'/BAG_session_runtime.png', dpi=300)
    plt.figure(0).clear()
    #plt.show()

    ## Krios assignment within BAG
    ax = df.plot.bar(x=xplot, y=["Krios1", "Krios2", "Krios3", "Krios4", "Talos", "GlaciosII"], stacked=True, color=['black','dimgrey','lightgrey','whitesmoke','royalblue','teal'], edgecolor = "silver", label=['Krios1', 'Krios2', 'Krios3', 'Krios4', 'Talos', 'Glacios2'], figsize=(figWidth,figHeight))
    # x-axis
    ax.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax.get_xticklabels():
       t.set_transform(t.get_transform()+trans)
    # y-axis
    ax.yaxis.set_label_position("left")
    ax.yaxis.tick_left()
    # Legend
    plt.legend(bbox_to_anchor=(0, 1.15, 1, 0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, frameon=False, framealpha=1.0)
    # Title
    plt.title("Instrument allocation distribution")
    plt.xlabel("BAG")
    plt.ylabel("Percentage")
    # Save fig
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig0 = ax.get_figure()
    fig0.savefig(str(output)+'/BAG_instrument_allocation.png', dpi=300)
    plt.figure(0).clear()

def plotSessionHistograms(csv, column, bins, fill, edge, title, yaxis, xaxis, outdir, filename):
    df = pd.read_csv(csv, parse_dates=['EPU session date','Session date'])

    print('Plotting histogram analysis of column: '+str(column))

    col = column
    suffix = str(filename)
    yaxis = str(yaxis)
    xaxis = str(xaxis)
    title = str(title)
    # Numeric
    df_plot = df
    df_plot[col] = pd.to_numeric(df_plot[col])
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=bins, color=fill,
            edgecolor=edge, linewidth=0.5).set(title=title)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.set_xticklabels(rotation=55)
    plt.tight_layout()
    ax10.figure.savefig(str(outdir)+'/'+str(suffix)+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)
    plt.close(10)

def plotSessionDates(csv, column, title, yaxis, xaxis, outdir, filename):
    df = pd.read_csv(csv, parse_dates=['EPU session date','Session date'])
    df["EPU session date"] = df["EPU session date"].astype("datetime64")
    df["Session date"] = df["Session date"].astype("datetime64")

    print('Plotting bar chart analysis of column: '+str(column))

    col = column
    suffix = str(filename)
    yaxis = str(yaxis)
    xaxis = str(xaxis)
    title = str(title)

    # Plot
    fig10 = plt.figure(10)
    # Convert session date to week of the year DEV DEV DEV
    df["EPU session week"] = df["EPU session date"].dt.isocalendar().week
    # Now group by new year week number
    df2 = df.groupby(["EPU session week"]).size().to_frame('Sessions per week').reset_index()

    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    sns.despine()
    # For bar charts go through this: 
    # https://learningactors.com/how-to-make-better-looking-charts-in-python/
    
    ax10 = sns.barplot(x = 'EPU session week', y = 'Sessions per week', data = df2)

    ax10.set(xlabel=xaxis, ylabel=yaxis)
    plt.xticks(rotation=70)
    plt.tight_layout()
    ax10.figure.savefig(str(outdir)+'/'+str(suffix)+'.png', dpi=300)

    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)
    plt.close(11)

def plotMicCount(f, output):
    print('plot mic count')

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    #df = pd.read_csv(f, parse_dates=['EPU session date'])
    df = pd.read_csv(f, parse_dates=['EPU session date'])
    # Dates are indeed annoying, convert to datetime and then to formatted string
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d")

    figHeight = 3
    figWidth = 12
    plotMicCount.figRatio = figWidth / figHeight

    fig1 = plt.figure(1)
    #ax0 = df.plot.bar(x="EPU session date", y="Total_EPU_mics", color='coral', figsize=(figWidth,figHeight))
    ax0 = df.plot.bar(x="EPU session date", y="Total_EPU_mics", color='coral', figsize=(figWidth,figHeight))

    ax0.set_xticklabels(labels=df['EPU session date'], rotation=70, rotation_mode="anchor", ha="right", fontsize=6)
    #ax0.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    #ax0.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    #ax0.set_ylabel("Total micrographs (x 1e3)")
    ax0.set_ylabel("Total micrographs")
    #fig1.add_axes([0.05, 0.1, 1, 0.8])

    plt.tight_layout()
    fig1 = ax0.get_figure()
    fig1.savefig(str(output)+'/micrograph_session_count.png', dpi=300)
    plt.figure(1).clear()
    plt.close(1)
    #plt.show()

def plotMicArea(f, output):
    print('plot mic area')
    
    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    #df = pd.read_csv(f, parse_dates=['EPU session date'])
    df = pd.read_csv(f, parse_dates=['EPU session date'])
    # Dates are indeed annoying, convert to datetime and then to formatted string
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d")

    figHeight = 3
    figWidth = 12
    plotMicCount.figRatio = figWidth / figHeight

    fig4 = plt.figure(4)
    ax3 = df.plot.bar(x="EPU session date", y="Total area (um^2)", color='coral', figsize=(figWidth,figHeight))
    ax3.set_xticklabels(labels=df['EPU session date'], rotation=70, rotation_mode="anchor", ha="right", fontsize=6)
    #ax=ax3.twinx()
    #ax3 = df.plot.scatter(x="EPU session date", y="Rate (um^2/hr)", color='purple', figsize=(figWidth,figHeight), ax=ax)
    ax3.set_ylabel("Total imaged area in microns")

    plt.tight_layout()
    fig4 = ax3.get_figure()
    fig4.savefig(str(output)+'/micrograph_session_area.png', dpi=300)
    plt.figure(4).clear()
    plt.close(4)
    #plt.show()

def BAGperformance(f, output):
    print('Session performance plots: session statistics')
    df = pd.read_csv(f, parse_dates=['EPU session date'])
    # Use datetime YY-MM-DD HH:MM to seperate sessions but plot with session name
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d %H:%M")

    figHeight = 4
    figWidth = 5
    plotAreaRateAnalysis.figRatio = figWidth / figHeight

    # Box plots
    fig10 = plt.figure(10)
    ax10 = df.boxplot(by = main.BAGname, column =['Grids'], grid = False, figsize=(5,4))
    ax10.set_ylabel("Number of Atlased grids")
    ax10.set_xlabel(main.xplot)
    ax10.yaxis.set_label_position("left")
    ax10.yaxis.tick_left()
    # x-axis
    ax10.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax10.get_xticklabels():
       t.set_transform(t.get_transform()+trans)

    title_boxplot = 'BAG report: Atlased grids'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    plt.xticks(ha='right')
    fig10 = ax10.get_figure()
    fig10.savefig(str(output)+'/micrograph_session_BAG_grids.png', dpi=300)
    plt.figure(10).clear()
    plt.close(10)
    #plt.show()

    # Box plots
    fig10 = plt.figure(10)
    ax10 = df.boxplot(by = main.BAGname, column =['Rate (um^2/hr)'], grid = False, figsize=(5,4))
    ax10.set_ylabel("Collection rate um^2/hr")
    ax10.set_xlabel(main.xplot)
    ax10.yaxis.set_label_position("right")
    ax10.yaxis.tick_right()
    # x-axis
    ax10.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax10.get_xticklabels():
       t.set_transform(t.get_transform()+trans)

    title_boxplot = 'BAG report: collection rate'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    plt.xticks(ha='right')
    fig10 = ax10.get_figure()
    fig10.savefig(str(output)+'/micrograph_session_BAG_rate.png', dpi=300)
    plt.figure(10).clear()
    plt.close(10)
    #plt.show()

    # Box plots
    fig10 = plt.figure(10)
    ax10 = df.boxplot(by = main.BAGname, column =['Collection time (hrs)'], grid = False, figsize=(5,4))
    ax10.set_ylabel("Collection run time")
    ax10.set_xlabel(main.xplot)
    ax10.yaxis.set_label_position("left")
    ax10.yaxis.tick_left()
    # x-axis
    ax10.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax10.get_xticklabels():
       t.set_transform(t.get_transform()+trans)

    title_boxplot = 'BAG report: collection run time'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    plt.xticks(ha='right')
    fig10 = ax10.get_figure()
    fig10.savefig(str(output)+'/micrograph_session_BAG_runtime.png', dpi=300)
    plt.figure(10).clear()
    plt.close(10)
    #plt.show()

    # Box plots
    fig10 = plt.figure(10)
    ax10 = df.boxplot(by = main.BAGname, column =['Total_EPU_mics'], grid = False, figsize=(5,4))
    ax10.set_ylabel("Total_EPU_mics")
    ax10.set_xlabel(main.xplot)
    ax10.yaxis.set_label_position("right")
    ax10.yaxis.tick_right()
    # x-axis
    ax10.tick_params(axis='x', labelrotation = 60)
    plt.xticks(ha='right')
    trans = mtrans.Affine2D().translate(20, 0)
    for t in ax10.get_xticklabels():
       t.set_transform(t.get_transform()+trans)

    title_boxplot = 'BAG report: total mics'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    plt.xticks(ha='right')
    fig10 = ax10.get_figure()
    fig10.savefig(str(output)+'/micrograph_session_BAG_mics.png', dpi=300)
    plt.figure(10).clear()
    plt.close(10)
    #plt.show()

def plotAreaRate(f, output):
    print('plot rate (um^2/hr)')

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    df = pd.read_csv(f, parse_dates=['EPU session date'])
    # Use datetime YY-MM-DD HH:MM to seperate sessions but plot with session name
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d")

    figHeight = 3
    figWidth = 12
    plotAreaRate.figRatio = figWidth / figHeight

    fig5 = plt.figure(5)
    ax4 = df.plot.bar(x="EPU session date", y="Rate (um^2/hr)", color='coral', figsize=(figWidth,figHeight))
    #ax4 = df.plot.scatter(x="EPU session date" , y="Rate (um^2/hr)", s=50, figsize=(figWidth,figHeight))
    #ax4 = df.plot(x="EPU session date", y="Rate (um^2/hr)", color='coral', figsize=(figWidth,figHeight), marker='o', linestyle='none')
    #ax2 = df.plot.bar(x="EPU session date" , y="Rate (p/hr)" , c="Mag (X)" , figsize=(8,4))

    ax4.set_xticklabels(labels=df['EPU session date'], rotation=70, rotation_mode="anchor", ha="right", fontsize=6)
    #ax4.set_xticklabels(labels=df['Visit'], rotation=70, rotation_mode="anchor", ha="right")
    ax4.set_ylabel("Collection rate (um^2/hr)")
    #fig5.add_axes([0.05, 0.1, 0.9, 0.8])

    plt.tight_layout()
    fig5 = ax4.get_figure()
    fig5.savefig(str(output)+'/micrograph_session_rate_area.png', dpi=300)
    plt.figure(5).clear()
    plt.close(5)
    #plt.show()

def plotProcessed(f, output):
    print('plot processed analyses')

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    df = pd.read_csv(f, parse_dates=['EPU session date'])
    # Use datetime YY-MM-DD HH:MM to seperate sessions but plot with session name
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d %H:%M")

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Ctf_min_A'
    # Filter data if required to remove outliers, quantile based filtering
    q_low = df[col].quantile(0.01)
    q_hi  = df[col].quantile(0.99)
    df_plot = df[(df[col] < q_hi) & (df[col] > q_low)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='plum',
            #binrange=(0.25,1.5),
            edgecolor="m", linewidth=2)
    ax10.set(xlabel='CTF best resolutions (Angstroms)', ylabel='Session count')
    ax10.figure.savefig(str(output)+'/processed_ctf_best_histogram.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Particles'
    filename = 'processed_particle_total_histogram'
    xaxis = 'Total number of particle images (1e6)'
    yaxis = 'Session count'
    # Filter data if required to remove outliers, quantile based filtering
    q_low = df[col].quantile(0.01)
    q_hi  = df[col].quantile(0.99)
    df_plot = df[(df[col] < q_hi) & (df[col] > q_low)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='paleturquoise',
            #binrange=(0.25,1.5),
            edgecolor="darkcyan", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    #col = 'Mask_ang'
    col = 'Ptcl_d_A'
    filename = 'processed_particle_size_histogram'
    xaxis = 'Particle size (Angstroms)'
    yaxis = 'Session count'
    # Filter data if required to remove outliers, quantile based filtering
    q_low = df[col].quantile(0.01)
    q_hi  = df[col].quantile(0.99)
    df_plot = df[(df[col] < q_hi) & (df[col] > q_low)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='paleturquoise',
            #binrange=(0.25,1.5),
            edgecolor="darkcyan", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Norm_av_ptcls_per_mic'
    filename = 'processed_particle_dist_histogram'
    xaxis = 'Average particle density (Normalised)'
    yaxis = 'Session count'
    # Filter data if required to remove outliers, quantile based filtering
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='paleturquoise',
            #binrange=(0.25,1.5),
            edgecolor="darkcyan", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Mean_ptcls_pct_clustered'
    filename = 'processed_particle_clustered_histogram'
    xaxis = 'Average particle clustering (percentage)'
    yaxis = 'Session count'
    # Filter data if required to remove outliers, quantile based filtering
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='paleturquoise',
            #binrange=(0.25,1.5),
            edgecolor="darkcyan", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Total_motion_max'
    filename = 'processed_total_motion_max_histogram'
    xaxis = 'Maximum total motion during exposure (Angstroms)'
    yaxis = 'Session count'
    # Filter data if required to remove outliers, quantile based filtering
    q_low = df[col].quantile(0.01)
    q_hi  = df[col].quantile(0.95)
    df_plot = df[(df[col] < q_hi) & (df[col] > q_low)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='tomato',
            #binrange=(0.25,1.5),
            edgecolor="darkred", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel=yaxis)
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

def plotSessions(f, output):
    print('plot session analyses')

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    df = pd.read_csv(f, parse_dates=['EPU session date'])
    # Use datetime YY-MM-DD HH:MM to seperate sessions but plot with session name
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d %H:%M")

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'apix (A/px)'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlabel='Detector pixel size (Angstroms)', ylabel='Count')
    ax10.figure.savefig(str(output)+'/session_pixel_size_histogram.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Total_EPU_mics'
    xaxis = 'Total number of micrographs collected'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/session_mics_histogram.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Supervisor size (GB)'
    filename = 'session_size_supervisor_histogram'
    xaxis = 'Supervisor data size (GB)'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='gold',
            #binrange=(0.25,1.5),
            edgecolor="orange", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Size (TB)'
    filename = 'session_size_mics_histogram'
    xaxis = 'Total dataset size (TB)'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Processed size (TB)'
    filename = 'session_size_processed_histogram'
    xaxis = 'Processed data size (TB)'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=40, color='green',
            #binrange=(0.25,1.5),
            edgecolor="darkgreen", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Rate (mic/hr)'
    filename = 'session_speed_histogram'
    xaxis = 'Collection rate (mic/hr)'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Average Foils per Square'
    filename = 'session_average_foils_per_square_histogram'
    xaxis = 'Average foil holes per square'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Shots per hole'
    filename = 'session_shots_per_hole_histogram'
    xaxis = 'Shots per foil hole'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='green',
            #binrange=(0.25,1.5),
            edgecolor="darkgreen", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)


    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Hole (um)'
    filename = 'session_foil_hole_size_histogram'
    xaxis = 'Specimen foil hole diameter (micron)'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='grey',
            #binrange=(0.25,1.5),
            edgecolor="dimgrey", linewidth=2)
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Total squares'
    filename = 'session_squares_total_histogram'
    xaxis = 'Number of squares available'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlim=(0, 600))
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Total squares'
    filename = 'session_squares_total_histogram'
    xaxis = 'Number of squares available'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlim=(0, 600))
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Targeted squares'
    filename = 'session_squares_targeted_histogram'
    xaxis = 'Number of targeted squares'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlim=(0, 300))
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Collected squares'
    filename = 'session_squares_collected_histogram'
    xaxis = 'Number of collected squares'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlim=(0, 300))
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Total dose (e/A2)'
    filename = 'session_total_dose_histogram'
    xaxis = 'Total dose on specimen (e-/Angstrom^2)'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='False',
            height=6, aspect=1.4, bins=18, color='green',
            #binrange=(0.25,1.5),
            edgecolor="darkgreen", linewidth=2)
    #ax10.set(xlim=(0, 300))
    ax10.set(xlabel=xaxis, ylabel='Count')
    ax10.figure.savefig(str(output)+'/'+filename+'.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

def plotOptics(f, output):
    print('plot optics analyses')

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    df = pd.read_csv(f, parse_dates=['EPU session date'])
    # Use datetime YY-MM-DD HH:MM to seperate sessions but plot with session name
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d %H:%M")

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    # Box plots
    fig10 = plt.figure(10)
    ax10 = df.boxplot(column =['Beam_um'], grid = False, figsize=(5,1), vert=False)
    ax10.set_ylabel("Microns")
    #ax10.yaxis.set_label_position("right")
    ax10.yaxis.tick_right()
    ax10.tick_params(axis='x', labelrotation = 45)
    title_boxplot = 'Beam diameter (um)'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    fig10 = ax10.get_figure()
    fig10.savefig(str(output)+'/optics_beam_diameter_box.png', dpi=300)
    plt.figure(10).clear()
    plt.close(10)
    #plt.show()

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Beam_um'
    # Filter data if required to remove outliers, quantile based filtering
    q_low = df[col].quantile(0.01)
    q_hi  = df[col].quantile(0.99)
    df_plot = df[(df[col] < q_hi) & (df[col] > q_low)]
    # Filter by absolute range
    #df_plot = df[df[col].between(0.1, 10)]
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='True',
            height=6, aspect=1.4, bins=18, color='slateblue',
            #binrange=(0.25,1.5),
            edgecolor="navy", linewidth=2)
    ax10.set(xlabel='Beam diameter (microns)', ylabel='Count')
    ax10.figure.savefig(str(output)+'/optics_beam_diameter_histogram.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Defocus max (um)'
    # Filter out strings
    #df_plot = df[df[col].str.contains('xml')==False]
    #df_plot[col] = pd.to_numeric(df_plot[col])
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='True',
            height=6, aspect=1.4, bins=20, color='orange',
            binrange=(-4,1),
            edgecolor="chocolate", linewidth=2)
    ax10.set(xlim=(-4, 1))
    ax10.set(xlabel='Defocus max (microns)', ylabel='Count')
    ax10.set_xticklabels(rotation=55)
    ax10.figure.savefig(str(output)+'/optics_defocus_max_histogram.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'Defocus min (um)'
    # Filter out strings 
    # DEV DEV DEV, if no strings are contained in the dataframe then this filter fails, need generic 'remove any string'
    #df_plot = df[df[col].str.contains('xml')==False]
    #df_plot[col] = pd.to_numeric(df_plot[col])
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='True',
            height=6, aspect=1.4, bins=20, color='orange',
            binrange=(-4,0),
            edgecolor="chocolate", linewidth=2)
    ax10.set(xlim=(-4, 1))
    ax10.set(xlabel='Defocus min (microns)', ylabel='Count')
    ax10.set_xticklabels(rotation=55)
    ax10.figure.savefig(str(output)+'/optics_defocus_min_histogram.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

    ## Histogram seaborn displot, be careful seaborn leaves parameters that affect all matplotlib 
    # Data
    col = 'apix (A/px)'
    df_plot = df
    # Plot
    fig10 = plt.figure(10)
    sns.set_style('white')
    sns.set_style("ticks")
    sns.set_context("talk") # This parameter appears to affect all matplotlib plots
    ax10 = sns.displot(data=df_plot, x=col, kind='hist', kde='True',
            height=6, aspect=1.4, bins=20, color='slateblue',
            #binrange=(-4,0),
            edgecolor="navy", linewidth=2)
    #ax10.set(xlim=(-4, 1))
    ax10.set(xlabel='Pixel size (Angstrom/px)', ylabel='Count')
    ax10.set_xticklabels(rotation=55)
    ax10.figure.savefig(str(output)+'/optics_apix_histogram.png', dpi=300)
    sns.set_context("paper") # This parameter appears to affect all matplotlib plots
    plt.figure(10).clear()
    plt.close(10)

def plotAreaRateAnalysis(f, output):
    print('plot area rate analyses')

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    df = pd.read_csv(f, parse_dates=['EPU session date'])
    # Use datetime YY-MM-DD HH:MM to seperate sessions but plot with session name
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d %H:%M")

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    figHeight = 2
    figWidth = 8.5
    plotAreaRateAnalysis.figRatio = figWidth / figHeight

    # Box plots

    df['apix (A/px)'] = df['apix (A/px)'].round(1)

    fig10 = plt.figure(10)
    ax10 = df.boxplot(by ='apix (A/px)', column =['Rate (um^2/hr)'], grid = False, figsize=(5,4))
    ax10.set_ylabel("Collection rate (um^2/hr)")
    ax10.yaxis.set_label_position("right")
    ax10.yaxis.tick_right()
    ax10.tick_params(axis='x', labelrotation = 45)
    title_boxplot = 'Rate (um^2/hr)'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    fig10 = ax10.get_figure()
    fig10.savefig(str(output)+'/micrograph_session_area_rate_box_plot.png', dpi=300)
    plt.figure(10).clear()
    plt.close(10)
    #plt.show()

    fig11 = plt.figure(11)
    ax11 = df.boxplot(by ='apix (A/px)', column =['Rate (mic/hr)'], grid = False, figsize=(5,4))
    ax11.set_ylabel("Collection rate (mic/hr)")
    ax11.yaxis.set_label_position("left")
    ax11.yaxis.tick_left()
    ax11.tick_params(axis='x', labelrotation = 45)
    title_boxplot = 'Rate (mic/hr)'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    fig11 = ax11.get_figure()
    fig11.savefig(str(output)+'/micrograph_session_mic_rate_box_plot.png', dpi=300)
    plt.figure(11).clear()
    plt.close(11)
    #plt.show()

    fig11 = plt.figure(11)
    ax11 = df.boxplot(by ='Shots per hole', column =['Rate (mic/hr)'], grid = False, figsize=(5,4))
    ax11.set_ylabel("Collection rate (mic/hr)")
    ax11.yaxis.set_label_position("left")
    ax11.yaxis.tick_left()
    ax11.tick_params(axis='x', labelrotation = 45)
    title_boxplot = 'Rate (mic/hr)'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    fig11 = ax11.get_figure()
    fig11.savefig(str(output)+'/micrograph_session_mic_rate_shots_box_plot.png', dpi=300)
    plt.figure(11).clear()
    plt.close(11)

    fig11 = plt.figure(11)
    ax11 = df.boxplot(by ='Shots per hole', column =['Rate (um^2/hr)'], grid = False, figsize=(5,4))
    ax11.set_ylabel("Collection rate (um^2/hr)")
    ax11.yaxis.set_label_position("right")
    ax11.yaxis.tick_right()
    ax11.tick_params(axis='x', labelrotation = 45)
    title_boxplot = 'Rate (um^2/hr)'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    fig11 = ax11.get_figure()
    fig11.savefig(str(output)+'/micrograph_session_area_rate_shots_box_plot.png', dpi=300)
    plt.figure(11).clear()
    plt.close(11)

    fig12 = plt.figure(12)
    ax12 = df.boxplot(by ='apix (A/px)', column =['Shots per hole'], grid = False, figsize=(5,4))
    ax12.set_ylabel("Shots per hole")
    ax12.yaxis.set_label_position("left")
    ax12.yaxis.tick_left()
    ax12.tick_params(axis='x', labelrotation = 45)
    title_boxplot = 'Shots per hole'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    fig12 = ax12.get_figure()
    fig12.savefig(str(output)+'/micrograph_session_mic_shotno_box_plot.png', dpi=300)
    plt.figure(12).clear()
    plt.close(12)
    #plt.show()

    fig13 = plt.figure(13)
    ax13 = df.boxplot(by ='apix (A/px)', column =['Hole (um)'], grid = False, figsize=(5,4))
    ax13.set_ylabel("Hole size (um)")
    ax13.yaxis.set_label_position("right")
    ax13.yaxis.tick_right()
    ax13.tick_params(axis='x', labelrotation = 45)
    title_boxplot = 'Hole (um)'
    plt.title( title_boxplot )
    plt.suptitle('')
    plt.tight_layout()
    fig13 = ax13.get_figure()
    fig13.savefig(str(output)+'/micrograph_session_mic_holesize_box_plot.png', dpi=300)
    plt.figure(13).clear()
    plt.close(13)
    #plt.show()

def construct_global_session_report(plots, csv, report):
    print('Generating global session PDF report...')

    # Column widths and heights
    colwidth = 50
    col0width = 65
    col1width = 50
    col2width = 44
    colH = 4
    col2H = 4.5
    col0H = 2 # This is a table spacer cell

    # https://towardsdatascience.com/how-to-create-pdf-reports-with-python-the-essential-guide-c08dd3ebf2ee
    # https://stackoverflow.com/questions/51864730/python-what-is-the-process-to-create-pdf-reports-with-charts-from-a-db
    pdf = FPDF()

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Global Session Report", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "All sessions", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    # Figure size params
    figHeight = 45
    figWidth = 180
    #figWidth = figHeight * plotTimings.figRatio

    # Basic reporting information
    pdf.cell(colwidth, colH, 'Input path:', 0, 0, 'L')
    # Will convert list into string with comma seperator
    path = ', '.join(args.input)
    pdf.multi_cell(0, colH, str(path), 0) # Use multicell so test will be wrapped

    pdf.cell(colwidth, colH, 'Output path:', 0, 0, 'L')
    path = os.path.realpath(args.output)
    pdf.multi_cell(0, colH, str(path), 0) # Use multicell so test will be wrapped

    pdf.cell(colwidth, colH, 'Time period of analysis (days):', 0, 0, 'L')
    pdf.cell(0, colH, str(main.periodReport), 0, 2, 'L') # New line via table row
    pdf.cell(0, colH*2, '', 0, 2, 'L') # Table space
    pdf.cell(-colwidth) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Start inserted graphs
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)
    #pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/global_session_time_setup_histogram.png', x = 5, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/global_session_time_collection_histogram.png', x = 105, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/global_session_counts_per_year_week.png', x = 5, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/global_session_mic_count_histogram.png', x = 105, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/global_session_area_rate_histogram.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/global_session_area_total_histogram.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_runtime.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_mics.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_grids.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_rate.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')

    pdf.cell(0, 1, '', 0, 2, 'L') # Table space

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Global Session Report", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "Visit & session analysis", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    ## Comprehensive data header
    # Visits and sessions
    pdf.cell(col2width, col2H, 'Visit directories found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.visitNoAll), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Sessions found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.totalEPU), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Filters
    pdf.cell(col2width, col2H, 'Exclude search term(s):', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(args.exclude), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Filtered sessions:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.excludeFilterNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Include search term:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(args.include), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Filtered sessions:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.includeFilterNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Plotted report
    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, '', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'EPU sessions with processed:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.processedNo), 1, 2, 'L') # New line via table row
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Start inserted graphs
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)
    #pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/visit_type_counts.png', x = 5, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/BAG_session_counts.png', x = 105, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/visit_type_runtime.png', x = 5, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/BAG_session_runtime.png', x = 105, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/BAG_instrument_allocation.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/BAG_session_type.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_runtime.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_mics.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_grids.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_rate.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')

    pdf.cell(0, 1, '', 0, 2, 'L') # Table space

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Global Session Report", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "Pipeline analysis", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    ## Comprehensive data header
    # Visits and sessions
    pdf.cell(col2width, col2H, 'Visit directories found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.visitNoAll), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Sessions found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.totalEPU), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Filters
    pdf.cell(col2width, col2H, 'Exclude search term(s):', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(args.exclude), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Filtered sessions:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.excludeFilterNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Include search term:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(args.include), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Filtered sessions:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.includeFilterNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Plotted report
    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, '', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'EPU sessions with processed:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.processedNo), 1, 2, 'L') # New line via table row
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Enter graphs
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)
    #pdf.cell(-colwidth) # Reset to 0 column position
    #pdf.image(str(plots)+'/visit_type_counts.png', x = 5, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/BAG_processed_counts.png', x = 105, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/visit_type_runtime.png', x = 5, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/BAG_session_runtime.png', x = 105, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/BAG_instrument_allocation.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/BAG_processed_counts.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Global Session Report", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "General box plot analysis", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    ## Comprehensive data header
    # Visits and sessions
    pdf.cell(col2width, col2H, 'Visit directories found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.visitNoAll), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Sessions found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.totalEPU), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Filters
    pdf.cell(col2width, col2H, 'Exclude search term(s):', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(args.exclude), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Filtered sessions:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.excludeFilterNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Include search term:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(args.include), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Filtered sessions:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.includeFilterNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Plotted report
    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, '', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'EPU sessions with processed:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.processedNo), 1, 2, 'L') # New line via table row
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Insert optics box plot graphs
    startx = pdf.get_x()
    starty = pdf.get_y()
    # Rate analysis
    figWidth = figHeight * plotAreaRateAnalysis.figRatio
    #pdf.cell(colwidth, colH, 'Collection rates analysis:', 0, 0, 'L')
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)
    #pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/micrograph_session_mic_rate_box_plot.png', x = 5, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/micrograph_session_area_rate_box_plot.png', x = 105, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/micrograph_session_mic_shotno_box_plot.png', x = 5, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/micrograph_session_mic_holesize_box_plot.png', x = 105, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/micrograph_session_mic_rate_shots_box_plot.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/micrograph_session_area_rate_shots_box_plot.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')

    pdf.cell(0, 1, '', 0, 2, 'L') # Table space

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Global Session Report", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "BAG box plot analysis", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    ## Comprehensive data header
    # Visits and sessions
    pdf.cell(col2width, col2H, 'Visit directories found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.visitNoAll), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Sessions found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.totalEPU), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Filters
    pdf.cell(col2width, col2H, 'Exclude search term(s):', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(args.exclude), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Filtered sessions:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.excludeFilterNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Include search term:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(args.include), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Filtered sessions:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.includeFilterNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Plotted report
    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, '', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'EPU sessions with processed:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.processedNo), 1, 2, 'L') # New line via table row
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Insert BAG analysis box plot graphs
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)
    #pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/micrograph_session_BAG_grids.png', x = 5, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/micrograph_session_BAG_mics.png', x = 105, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/micrograph_session_BAG_runtime.png', x = 5, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/micrograph_session_BAG_rate.png', x = 105, y = starty+(80*1), w = 100, h = 80, type = '', link = '')

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Global Session Report", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "Script and session advisories", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    ## Capture script execution command
    pdf.cell(colwidth, colH, 'Script:', 0, 2, 'L')
    pdf.multi_cell(0, colH, str(main.command), 1) # New line via table row
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    #pdf.cell(-colwidth) # Reset to 0 column position

    ## Gather advisories on sessions
    # Add advisories as table if any are present
    advisorycsv=(csv+'/'+str(main.outputName)+'_advisory.csv')

    if os.path.exists(advisorycsv):
       # All data main stats
       pdf.cell(colwidth, colH, 'All session advisories summary:', 0, 2, 'L')
       #pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row

       # Essential code, this is where the csv data is read in for populated the data table
       # Read in csv, keep NaN so you can handle these as strings, force all floats to read as string
       df = pd.read_csv(csv+'/'+str(main.outputName)+'_advisory.csv', parse_dates=["EPU session date"], keep_default_na=False, na_values=[''])#, dtype=str)
       row = df.shape[0]  # Gives number of rows
       col = df.shape[1]  # Gives number of columns
       df['EPU session date'] = pd.to_datetime(df['EPU session date'])
       df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d")
       # Sort by date in descending order
       df.sort_values(by='EPU session date', ascending = False, inplace=True)
       #print(df['EPU session date'])

       # This dataframe is used to create the global report table, you can add new lines to lookup data in the csv files
       reportTable = [{'name':'Proposal', 'unit':'', 'data':'Proposal', 'width':12, 'height':3, 'line':'1', 'return':'0', 'align':'R'},
                     {'name':'BAG', 'unit':'', 'data':'BAG', 'width':14, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                     {'name':'EPU session date', 'unit':'', 'data':'EPU session date', 'width':20, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                     {'name':'Visit', 'unit':'', 'data':'Visit', 'width':16, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                     {'name':'Scope', 'unit':'', 'data':'Instrument', 'width':10, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                     {'name':'Type', 'unit':'', 'data':'Advisory_type', 'width':10, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                     {'name':'Advisory', 'unit':'', 'data':'Advisory', 'width':55, 'height':3, 'line':'1', 'return':'2', 'align':'L'}]
       reportDf = pd.DataFrame(reportTable)

       reportDf = pd.DataFrame(reportTable)
       # Reset needs to be total cell width, less one cell width
       widthSum = reportDf['width'].sum()
       widthLast = reportDf['width'].iloc[-1]
       reset = (widthSum-widthLast)*-1

       # Populate all data list table
       pdf.set_font('arial','B', 6)

       #
       print('Gathering all data for global report advisories table')

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
            insert = str(data)

            pdf.cell(int(reportDf.iloc[j]['width']),int(reportDf.iloc[j]['height']),insert,int(reportDf.iloc[j]['line']),int(reportDf.iloc[j]['return']),str(reportDf.iloc[j]['align']))
            j = j+1
         pdf.cell(reset)
         h = h+1

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Global Session Report", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "All session data summary", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    ## Comprehensive data header with additional error reporting
    # Visits and sessions
    pdf.cell(col2width, col2H, 'Visit directories found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.visitNoAll), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Sessions found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.totalEPU), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Calculate percentages
    emptyPct = str(roundup(((float(main.sessionErrorEmpty) / main.sessionNo)*100),1))
    nosupPct = str(roundup(((float(main.sessionErrorSuper) / main.sessionNo)*100),1))
    noatlPct = str(roundup(((float(main.sessionErrorAtlas) / main.sessionNo)*100),1))
    noepuPct = str(roundup(((float(main.sessionErrorEpu) / main.sessionNo)*100),1))
    tomoPct = str(roundup(((float(main.sessionErrorTomo) / main.sessionNo)*100),1))
    noerrPct = str(roundup(((float(main.sessionErrorNone) / main.sessionNo)*100),1))
    procPct = str(roundup(((float(main.processedNo) / main.sessionNo)*100),1))

    spaPct = str(float(noatlPct)+float(noepuPct)+float(noerrPct))
    unkPct = str(nosupPct)

    # Errors
    pdf.cell(col2width, col2H, 'Session classification:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Metadata breakdown:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, 'SPA %', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(spaPct+' %'), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'No EPU metadata errors', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionErrorNone)+' ('+noerrPct+'%)', 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, 'Tomo %', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(tomoPct+' %'), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'No Atlas dm', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionErrorAtlas)+' ('+noatlPct+'%)', 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, 'Unknown %', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(unkPct+' %'), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'No EPU dm', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionErrorEpu)+' ('+noepuPct+'%)', 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, 'Empty session directory %', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(emptyPct+' %'), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'No Supervisor directories', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionErrorSuper)+' ('+nosupPct+'%)', 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, '-', 1, 0, 'L')
    pdf.cell(col2width, col2H, str('-'), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Is Tomo', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionErrorTomo)+' ('+tomoPct+'%)', 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    pdf.cell(col2width, col2H, '-', 1, 0, 'L')
    pdf.cell(col2width, col2H, str('-'), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Empty session directory', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionErrorEmpty)+' ('+emptyPct+'%)', 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, '', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Processed:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.processedNo)+' ('+procPct+'%)', 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3-4) # Reset to 0 column position

    # Essential code, this is where the csv data is read in for populated the data table
    # Read in csv, keep NaN so you can handle these as strings, force all floats to read as string
    df = pd.read_csv(csv+'/'+str(main.outputName)+'_session.csv', parse_dates=["EPU session date"], keep_default_na=False, na_values=[''])#, dtype=str)
    row = df.shape[0]  # Gives number of rows
    col = df.shape[1]  # Gives number of columns
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d")
    # Sort by date in descending order
    df.sort_values(by='EPU session date', ascending = False, inplace=True)
    #print(df['EPU session date'])

    # This dataframe is used to create the global report table, you can add new lines to lookup data in the csv files
    reportTable = [{'name':'Visit', 'unit':'', 'data':'Visit', 'width':14, 'height':3, 'line':'1', 'return':'0', 'align':'R'},
                  {'name':'BAG', 'unit':'', 'data':'BAG', 'width':11, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Date', 'unit':'(end)', 'data':'EPU session date', 'width':13, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Scope', 'unit':'', 'data':'Instrument', 'width':8, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Type', 'unit':'', 'data':'Type', 'width':7, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Error', 'unit':'', 'data':'Error', 'width':10, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Grids', 'unit':'#', 'data':'Grids', 'width':7, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Setup', 'unit':'hrs', 'data':'Setup time (hrs)', 'width':7, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Run', 'unit':'hrs', 'data':'Collection time (hrs)', 'width':7, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Sqrs', 'unit':'#', 'data':'Collected squares', 'width':7, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'speed', 'unit':'mic/hr', 'data':'Rate (mic/hr)', 'width':8, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  #{'name':'mics', 'unit':'#', 'data':'Total_EPU_mics', 'width':9, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'movies', 'unit':'#', 'data':'Total_movies', 'width':9, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'rate', 'unit':'um2/hr', 'data':'Rate (um^2/hr)', 'width':8, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  #{'name':'~rate', 'unit':'~um2/hr', 'data':'~rate', 'width':9, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'total', 'unit':'um2', 'data':'Total area (um^2)', 'width':8, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'pixel', 'unit':'a/pix', 'data':'apix (A/px)', 'width':8, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Hole', 'unit':'um', 'data':'Hole (um)', 'width':6, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Shots', 'unit':'#', 'data':'Shots per hole', 'width':7, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Beam', 'unit':'um', 'data':'Beam (um)', 'width':9, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Time', 'unit':'sec', 'data':'Exposure time (s)', 'width':8, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Fdose', 'unit':'e/A^2', 'data':'Fraction dose (e/A2)', 'width':8, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Dose', 'unit':'', 'data':'Dose', 'width':11, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Target', 'unit':'', 'data':'AFIS/Accurate', 'width':9, 'height':3, 'line':'1', 'return':'0', 'align':'C'},
                  {'name':'Proc', 'unit':'', 'data':'Processed', 'width':7, 'height':3, 'line':'1', 'return':'2', 'align':'C'}]
    reportDf = pd.DataFrame(reportTable)
    # Reset needs to be total cell width, less one cell width
    widthSum = reportDf['width'].sum()
    widthLast = reportDf['width'].iloc[-1]
    reset = (widthSum-widthLast)*-1

    # Populate all data list table
    pdf.set_font('arial','B', 6)

    print('Gathering all data for global report data table')

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
         if data == 'nan' or data == 'NaN':
            insert = '-'
         # DEV DEV the detection of numbers and rounding is failing for all but mics column...
         elif data.isnumeric():
            if colName == 'mics':
               insert = str(round(float(data)))
            else:
               insert = str(roundup(float(data),2))
         else:
            insert = str(data)

         if data == 'Unknown':
            insert = '-'

         pdf.cell(int(reportDf.iloc[j]['width']),int(reportDf.iloc[j]['height']),insert,int(reportDf.iloc[j]['line']),int(reportDf.iloc[j]['return']),str(reportDf.iloc[j]['align']))
         j = j+1
      pdf.cell(reset)
      h = h+1

    ## Write out report
    pdf.output(report+'/'+str(main.outputName)+'_session.pdf', 'F')

    # Report
    print('Generated PDF report in '+report+'/'+str(main.outputName)+'_session.pdf')
    print('')

def plotProcessedAnalysis(f, output):
    print('plot processed analyses')

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    df = pd.read_csv(f)
    # Use datetime YY-MM-DD HH:MM to seperate sessions but plot with session name
    df['EPU session date'] = pd.to_datetime(df['EPU session date'])
    df['EPU session date'] = df['EPU session date'].dt.strftime("%Y-%m-%d %H:%M")

    # Make sure no other plots interfere with this plotting
    plt.clf
    plt.close('all')

    # Piking and extraction params
    # Marker styles: https://matplotlib.org/stable/api/markers_api.html
    fig14 = plt.figure(14)
    ax14 = df.plot.scatter(x="EPU session date" , y="Box_ang", s=25, figsize=(8,4), color='black', label="Box_ang", marker='s')
    ax15 = df.plot.scatter(x="EPU session date" , y="Mask_ang", s=20, figsize=(8,4), color='white', edgecolors='red', label="Mask_ang", marker='o', ax=ax14)
    ax16 = df.plot.scatter(x="EPU session date" , y="Pick_max", s=20, figsize=(8,4), color='orange', label="Pick_max", ax=ax14)
    ax17 = df.plot.scatter(x="EPU session date" , y="Pick_min", s=10, figsize=(8,4), color='green', label="Pick_min", ax=ax14)
    ax14.axes.xaxis.set_visible(False)
    ax14.set_ylabel("Angstrom")
    plt.tight_layout()
    fig14 = ax14.get_figure()
    fig14.savefig(str(output)+'/processed_picking.png', dpi=300)
    plt.figure(14).clear()
    plt.close(14)

    # Micrograph, motion and ctf counts
    fig18,ax = plt.subplots()
    ax.axes.xaxis.set_visible(False)
    df.plot.scatter(x="EPU session date" , y="Total_EPU_mics", s=50, figsize=(8,4), color='black', label="Total micrographs", marker='o', ax=ax)
    df.plot.scatter(x="EPU session date" , y="MotionCorr", s=30, figsize=(8,4), color='blue', label="MotionCorr", marker='o', ax=ax)
    df.plot.scatter(x="EPU session date" , y="CtfFind", s=15, figsize=(8,4), color='none', edgecolors='cyan', label="CtfFind", marker='o', ax=ax)
    ax.set_ylabel("Total number", color='green')
    plt.tight_layout()
    fig18.savefig(str(output)+'/processed_preprocessing.png', dpi=300)
    plt.figure(18).clear()
    plt.close(18)

    #fig18 = plt.figure(18)
    #ax18 = df.plot.scatter(x="EPU session date" , y="MotionCorr", s=30, figsize=(8,4), color='blue', label="MotionCorr", marker='o')
    #ax19 = df.plot.scatter(x="EPU session date" , y="CtfFind", s=15, figsize=(8,4), color='none', edgecolors='cyan', label="CtfFind", marker='o', ax=ax18)
    #ax20 = df.plot.scatter(x="EPU session date" , y="Pick_max", s=20, figsize=(8,4), color='orange', label="Pick_max", ax=ax18)
    #ax21 = df.plot.scatter(x="EPU session date" , y="Pick_min", s=10, figsize=(8,4), color='green', label="Pick_min", ax=ax18)
    #ax18.axes.xaxis.set_visible(False)
    #ax18.set_ylabel("Total number")
    #plt.tight_layout()
    #fig18 = ax18.get_figure()
    #fig18.savefig(str(output)+'/processed_preprocessing.png', dpi=300)
    #plt.figure(18).clear()
    #plt.close(18)

    # Motion data - DEV DEV THIS IS A MUCH NEATER WAY TO PLOT
    fig22,ax = plt.subplots()
    ax.axes.xaxis.set_visible(False)
    df.plot.scatter(x="EPU session date" , y="Total_motion_max", s=30, figsize=(8,4), color='green', label="Total motion (best)", marker='o', ax=ax)
    ax.set_ylabel("Total motion best (Angstroms)", color='green')
    ax2=ax.twinx()
    df.plot.scatter(x="EPU session date" , y="Total_motion_min", s=50, figsize=(8,4), color='orange', label="Total motion (worst)", marker='o', ax=ax2)
    ax2.set_ylabel("Total motion worst (Angstroms)", color='orange')
    plt.tight_layout()
    fig22.savefig(str(output)+'/processed_total_motion.png', dpi=300)
    plt.figure(22).clear()
    plt.close(22)

    #fig22 = plt.figure(22)
    #ax22 = df.plot.scatter(x="EPU session date" , y="Total_motion_max", s=30, figsize=(8,4), color='orange', label="Total motion (max)", marker='o')
    #ax23 = ax22.twinx()
    #ax23 = df.plot.scatter(x="EPU session date" , y="Total_motion_min", s=30, figsize=(8,4), color='green', label="Total motion (min)", marker='o', ax=ax22)
    #ax22.axes.xaxis.set_visible(False)
    #ax22.set_ylabel("Total motion (Angstroms)")
    #
    #fig22 = ax22.get_figure()
    #fig22.savefig(str(output)+'/processed_total_motion.png', dpi=300)
    #plt.figure(22).clear()
    #plt.close(22)

    # Ctf data
    fig24,ax = plt.subplots()
    ax.axes.xaxis.set_visible(False)
    df.plot.scatter(x="EPU session date" , y="Ctf_min_A", s=30, figsize=(8,4), color='blue', label="Max CTF res (best)", marker='o', ax=ax)
    ax.set_ylabel("CTF resolution best (Angstroms)", color='blue')
    ax2=ax.twinx()
    df.plot.scatter(x="EPU session date" , y="Ctf_max_A", s=50, figsize=(8,4), color='red', label="Max CTF res (worst)", marker='o', ax=ax2)
    ax2.set_ylabel("CTF resolution worst (Angstroms)", color='red')
    plt.tight_layout()
    fig24.savefig(str(output)+'/processed_ctf_res.png', dpi=300)
    plt.figure(24).clear()
    plt.close(24)

    #fig24 = plt.figure(24)
    #ax24 = df.plot.scatter(x="EPU session date" , y="Ctf_max_A", s=30, figsize=(8,4), color='red', label="Max CTF res (worst)", marker='o')
    #ax25 = df.plot.scatter(x="EPU session date" , y="Ctf_min_A", s=30, figsize=(8,4), color='blue', label="Max CTF res (best)", marker='o', ax=ax24)
    #ax24.axes.xaxis.set_visible(False)
    #ax24.set_ylabel("Resolution (Angstroms)")
    #plt.tight_layout()
    #fig24 = ax24.get_figure()
    #fig24.savefig(str(output)+'/processed_ctf_res.png', dpi=300)
    #plt.figure(24).clear()
    #plt.close(24)

def construct_global_process_report(plots, csv, report):
    print('Generating global processed PDF report...')

    # Column widths and heights
    colwidth = 50
    col0width = 65
    col1width = 50
    col2width = 44
    colH = 4
    col2H = 4.5
    col0H = 2 # This is a table spacer cell

    # https://towardsdatascience.com/how-to-create-pdf-reports-with-python-the-essential-guide-c08dd3ebf2ee
    # https://stackoverflow.com/questions/51864730/python-what-is-the-process-to-create-pdf-reports-with-charts-from-a-db
    pdf = FPDF()

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Processed Reporting and Analysis", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "All sessions", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

## Comprehensive data header with additional error reporting
    # Visits and sessions
    pdf.cell(col2width, col2H, 'Visit directories found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.visitNoAll), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Sessions found:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.totalEPU), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position
    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    pdf.cell(col2width, col2H, 'Tabulated unique visits:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(combineCsv.BAGvisits), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Tabulated sessions', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.sessionNo), 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Calculate percentages
    emptyPct = str(roundup(((float(main.sessionErrorEmpty) / main.sessionNo)*100),1))
    nosupPct = str(roundup(((float(main.sessionErrorSuper) / main.sessionNo)*100),1))
    noatlPct = str(roundup(((float(main.sessionErrorAtlas) / main.sessionNo)*100),1))
    noepuPct = str(roundup(((float(main.sessionErrorEpu) / main.sessionNo)*100),1))
    tomoPct = str(roundup(((float(main.sessionErrorTomo) / main.sessionNo)*100),1))
    noerrPct = str(roundup(((float(main.sessionErrorNone) / main.sessionNo)*100),1))
    procPct = str(roundup(((float(main.processedNo) / main.sessionNo)*100),1))

    spaPct = str(float(noatlPct)+float(noepuPct)+float(noerrPct))
    unkPct = str(nosupPct)

    # Processed statistics
    pdf.cell(col2width, col2H, '', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    pdf.cell(col2width, col2H, 'Processed:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(main.processedNo)+' ('+procPct+'%)', 1, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3) # Reset to 0 column position

    # column 1
    pdf.cell(col2width, col2H, 'Processed characterisation:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    # column 2
    pdf.cell(col2width, col2H, '', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 2, 'L')
    # Reset to 0 column position
    pdf.cell(-col2width*3)
    # column 1
    pdf.cell(col2width, col2H, 'Requested 2D:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    # column 2
    pdf.cell(col2width, col2H, 'Class 2D ran:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 2, 'L')
    # Reset to 0 column position
    pdf.cell(-col2width*3)
    # column 1
    pdf.cell(col2width, col2H, 'Requested 3D:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 0, 'L')
    # column 2
    pdf.cell(col2width, col2H, 'Class 3D ran:', 1, 0, 'L')
    pdf.cell(col2width, col2H, str(''), 1, 2, 'L')
    # Reset to 0 column position
    pdf.cell(-col2width*3)

    # Table break space
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 0, 'L')
    pdf.cell(col2width, col0H, '', 0, 2, 'L')
    pdf.cell(-col2width*3-4) # Reset to 0 column position
    
    # Start inserted graphs
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)
    #pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/session_size_supervisor_histogram.png', x = 5, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/global_session_time_collection_histogram.png', x = 105, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/session_size_mics_histogram.png', x = 5, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/global_session_mic_count_histogram.png', x = 105, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/session_size_processed_histogram.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/global_session_area_total_histogram.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_runtime.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_mics.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_grids.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_rate.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')

    #pdf.cell(col0width, colH, 'Number of CTFFIND star lines:', 0, 0, 'L')
    #pdf.cell(0, colH, str(analyseProcessed.ctfFindNo), 0, 2, 'L')
    #pdf.cell(-col0width) # Reset column position
    #pdf.cell(col0width, colH, 'Min CTF resolution fit (A):', 0, 0, 'L')
    #pdf.cell(0, colH, str(analyseProcessed.ctfMin), 0, 2, 'L')
    #pdf.cell(-col0width) # Reset column position
    #pdf.cell(col0width, colH, 'Max CTF resolution fit (A):', 0, 0, 'L')
    #pdf.cell(0, colH, str(analyseProcessed.ctfMax), 0, 2, 'L')
    #pdf.cell(-col0width, colH, '', 0, 2, 'L') # Table space
    #pdf.cell(-col0width) # Reset column position

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Processed Reporting and Analysis", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "All sessions", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    # Start inserted graphs
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)
    #pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/processed_particle_size_histogram.png', x = 5, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/processed_particle_dist_histogram.png', x = 105, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/processed_particle_total_histogram.png', x = 5, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/processed_particle_clustered_histogram.png', x = 105, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/session_size_processed_histogram.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/global_session_area_total_histogram.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_runtime.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_mics.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_grids.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_rate.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Processed Reporting and Analysis", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "All sessions", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)

    # Start inserted graphs
    startx = pdf.get_x()
    starty = pdf.get_y()
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-20)
    #pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/processed_ctf_best_histogram.png', x = 5, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/processed_particle_size_histogram.png', x = 105, y = starty+(80*0), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/processed_total_motion_max_histogram.png', x = 5, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    pdf.image(str(plots)+'/processed_particle_total_histogram.png', x = 105, y = starty+(80*1), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/session_size_processed_histogram.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/global_session_area_total_histogram.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_runtime.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_mics.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_grids.png', x = 5, y = starty+(80*2), w = 100, h = 80, type = '', link = '')
    #pdf.image(str(plots)+'/micrograph_session_BAG_rate.png', x = 105, y = starty+(80*2), w = 100, h = 80, type = '', link = '')

    # New page for template
    pdf.add_page()
    pdf.set_xy(10, 0)
    pdf.set_font('arial','B', 14)
    pdf.cell(190, 5, "", 0, 2)
    pdf.cell(190, 8, "eBIC Processed Reporting and Analysis", 0, 2, 'R')
    pdf.set_font('arial','I', 8)
    pdf.cell(190, 6, "All sessions", 0, 2, 'R')
    # Establish initial position for text
    pdf.set_xy(0, 22)
    pdf.set_font('arial','', 8)
    pdf.cell(10,0)
    # Basic preprocessing stats
    pdf.cell(colwidth, colH, 'Preprocessing: motion correction and ctf estimation:', 0, 0, 'L')
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/processed_preprocessing.png', x = 30, y = None, w = 160, h = 80, type = '', link = '')
    pdf.cell(0, 1, '', 0, 2, 'L') # Table space

    # Basic preprocessing stats
    pdf.cell(colwidth, colH, 'Total motion (worst and best):', 0, 0, 'L')
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/processed_total_motion.png', x = 30, y = None, w = 160, h = 80, type = '', link = '')
    pdf.cell(0, 1, '', 0, 2, 'L') # Table space

    # Basic preprocessing stats
    pdf.cell(colwidth, colH, 'CTF max resolutions (worst and best):', 0, 0, 'L')
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/processed_ctf_res.png', x = 30, y = None, w = 160, h = 80, type = '', link = '')
    pdf.cell(0, 1, '', 0, 2, 'L') # Table space

    # Basic picking stats
    pdf.cell(colwidth, colH, 'Picking and extraction:', 0, 0, 'L')
    pdf.cell(0, colH, '', 0, 2, 'L') # New line via table row
    pdf.cell(-colwidth) # Reset to 0 column position
    pdf.image(str(plots)+'/processed_picking.png', x = 30, y = None, w = 160, h = 80, type = '', link = '')
    pdf.cell(0, 1, '', 0, 2, 'L') # Table space

    # Write out report
    pdf.output(report+'/'+str(main.outputName)+'_processed.pdf', 'F')

    # Report
    print('Generated PDF report in '+report+'/'+str(main.outputName)+'_processed.pdf')
    print()

def dev():
    print('')

# Run main
main()
