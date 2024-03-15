# binSLT
Contains: 
  Description of the code
  Future updates
  How to install
  
Description:
Code to extract spin-orbit induced predissociation lifetimes from stabilization method Duo *.states files.

Currently only supports doublet sigma predissociation states as this was developed explicitly to extract predissociation lifetimes for the A2Sigma+ states of OH, SH, and SD. 

Future Updates:
State labels currently constrained to J, v, and e/f parity but will be expanded in future to generalise for all state symmetries. 

How to install:
This method works for linux machines and should work for Mac OS. If you are running Windows, Windows Subsystem for Linux is recomended here. 

In terminal create a directory named 'binslt' and download this repository. Enter this directory, running:
  **$ls**
Should return:
  **binslt  dist  examples**  pyproject.toml  README.md
and run in terminal:
  **$pip install .**

Don't skip the "**.**"

If you wish to make changes to the code while using it:

  **pip install -e .**

Any problems can be reported to **georgi.mitev.16@ucl.ac.uk**
