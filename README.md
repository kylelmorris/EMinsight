## EMinsight

# Background

These scripts will parse the xml and mrc files in an EPU session directory and produce a report on the sessions set up parameters as well as performance analysis of the data collection sessions. It can be run anytime during a data collection, but should ideally be run at the end of a session to accuractely report on the whole session.

# Installation

Clone the repository to your local machine

In a terminal, change directory into the cloned repository directory

Run the following command to set up a vritual environment for EMinsight

$ python -m venv python

$ source python/bin/activate

$ pip install .

Each time you want to use EMinsight, remember to activate that virtual environment

$ source python/bin/activate

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

dls <br />
├── m02 # eBIC/DLS instrument identifier <br />
│   ├── data <br />
│   │   ├── 2023 # year of collection <br />
|   |   |   └── bi23047-106 # unique eBIC/DLS data collection identifer <br />
|   |   |   |   ├── Atlas <br />
|   |   |   |   |   └── Supervisor\_[date]_[time]-[session-identifier]_Atlas # Atlas directory used by EPU session <br />
|   |   |   |   ├── processed <br />
|   |   |   |   |   └── raw
|   |   |   |   |   |   └── Relion <br />
|   |   |   |   ├── processing <br />
|   |   |   |   ├── raw <br />
|   |   |   |   |   ├── GridSquare_[Square-ID] # Raw movie data from EPU session <br />
|   |   |   |   |   |   └── Data # Raw data from this GridSquare <br />
|   |   |   |   |   ├── metadata <br />
|   |   |   |   |   |   └── Supervisor_[date]_[time]_[unique-session-identifier]_EPU # EPU session files, xml metadata, compressed images <br />
|   |   |   |   |   |   |   └── Images-Disc1 <br />
|   |   |   |   |   |   |   |   ├── Data <br />
|   |   |   |   |   |   |   |   ├── FoilHoles <br />

Further development is required to make this robust to scenarios where other facilities use: <br />
    - Different data collection identifiers (i.e. bi23047-106 versus YYYYMMDD) <br />
    - Different instrument identifiers (i.e. m02 versus Krios2) <br />
    - Different directory structures altogether <br />

    