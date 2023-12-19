# EMinsight

# Installation

Clone the repository to your local machine
The following instructions assume that you downlaod to your home directory.

In a terminal, change directory into the cloned repository directory

Run the following command to set up a vritual environment for EMinsight

$ cd ~/EMinsight

$ python -m venv python

$ source python/bin/activate

$ pip install .

Each time you want to use EMinsight, remember to activate that virtual environment

$ source ~/EMinsight/python/bin/activate

You may want to add EMinsight to your ~/.bash_profile, ~/.zshrc or equivalent

export PATH="~/EMinsight/bin:$PATH"