README
Guide for getting started with NREL's design optimization and analysis tools for wind turbines and plants

TOOLS
--------------------
Tools you will need for conceptual design of monopiles and jacket:
- Python (preferably via Anaconda)
- WISDEM (NREL's low-fidelity conceptual design optimization tool for wind turbines and plant)

Tools you will need for automated fatigue analysis (monopiles only) using OpenFAST (time-domain aero-servo-hydro-elastic multi-body dynamics tool)
- Python (preferably via Anaconda)
- Compilers (local to your machine or via Anaconda)
- WEIS (Essentially a stack of WISDEM, OpenFAST, ROSCO controller, and multiple pre- and post-processing utilities)


INSTALLATION
--------------------
WISDEM:
- Simplest option is: "conda install wisdem"
- To interact with the code and the example files, follow the instructions at: https://github.com/WISDEM/WISDEM
WISDEM will work on all operating systems.  Documentation is located at:
https://wisdem.readthedocs.io/en/master/

WEIS:
- Have to install from source following instructions at: https://github.com/WISDEM/WEIS
WEIS only works on Mac and Linux for now.  Windows users can use the Windows Subsystem Linux (preferred) or Cygwin (not well tested).  Documentation is still limited and located at:
https://weis.readthedocs.io/en/latest/


USAGE
--------------------
This archive contains two folders, one for monopile design and one for jacket design.  Both have top level drivers, "monopile_driver.py" and "jacket_driver.py".  Towards the top of these files are the user options that define what type of analysis is being run and key design options.  You can activate all options to run everything all at once, or run through things one at a time.  These would be run simply by doing:
"python monopile_driver.py" at the Anaconda prompt.
