APBS 1.4.1-binary README :tada::clap:
=====================================

Contained herein are the steps for extracting and running the APBS 1.4.1
binary release.

First, a few notes are in order.

- Only 64-bit operating systems are supported with this release.
- [FETK](http://www.fetk.org) support is included for the OS X and Linux builds, but not for the Windows build.  This is mostly to do with the inability of Windows and Autotools to play nicely together.  _NB: FETK has been modified somewhat to better inter-operate with APBS._
- Intel C and C++ compilers (version 14.0.3.20140415) and libraries were used for this release.


OS-Specific Steps
-----------------
Below are the specific steps used to unpack/install the APBS release.  Following the steps below for your specific OS should produce the expected results. 


OS X & Linux tar files:
----------------------

1. tar -zxvf apbs-1.4.1-binary.tar.gz

OSX App Bundle:
--------------

1.  Drag the APBS bundle to your Applications directory
2.  From your Applications folder, double click on the APBS icon.  This
starts a new shell and adds apbs to your path.

Windows:
-------

1. double click on the .exe and follow the instructions


To run the APBS binary, see the instructions here:
http://www.poissonboltzmann.org/apbs-pdb2pqr/docs/calculating/
