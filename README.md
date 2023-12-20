## EMinsight

# Background

These scripts will parse the xml and mrc files in an EPU session directory and produce a report on the sessions set up parameters as well as performance analysis of the data collection sessions. It can be run anytime during a data collection, but should ideally be run at the end of a session to accuractely report on the whole session.

# Installation

Clone the repository to your local machine

In a terminal, change directory into the cloned repository directory

Run the following command to set up a vritual environment for EMinsight

$ python -m EMinsight venv

$ source EMinsight/bin/activate

$ pip install .

Each time you want to use EMinsight, remember to activate that virtual environment

$ source EMinsight/bin/activate

# Usage instructions

Navigate to the directory you want to send outputs to (optional)
$ cd /path/to/working/output/area

Run a report on a single visit diretory
$ epu.xml_summary.py --help
$ epu.xml_summary.py --i /dls/m07/data/2023/bi23047-106 --o bi23047-106

Run a global report across all visits from all visits on m07 in 2023
$ epu.reporting.py --o global --input /dls/m07/data/2023

Use a GUI to achieve the same as above
$ eminsight.py

# Data structure

EMinsight currently expects data structured in the following way

|-dls
  |-m02 # eBIC/DLS instrument identifier
    |-data
      |-2023 # year of collection
        |-bi23047-106 # unique eBIC/DLS data collection identifer
          |-Atlas
            |-Supervisor_[date]_[time]-[session-identifier]_Atlas # Atlas directory used by EPU session
          |-processed # Preprocessing pipeline results
          |-processing # User processing area
          |-raw
            |-GridSquare_[Square-ID] # Raw movie data from EPU session
              |-Data # Raw data from this GridSquare
            |-metadata
              |-Supervisor_[date]_[time]_[unique-session-identifier]_EPU # EPU session files, xml metadata, compressed images
                |-Images-Disc1
                  |-GridSquare_[Square-ID]
                    |-Data
                    |-FoilHoles

Further development is required to make this robust to scenarios where other facilities use:
    Different data collection identifiers (i.e. bi23047-106 versus YYYYMMDD)
    Different instrument identifiers (i.e. m02 versus Krios2)
    Different directory structures altogether

    