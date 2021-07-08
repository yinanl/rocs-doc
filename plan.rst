Work In Progress
================


LTL-to-BA translator embedding
------------------------------
The fact that the user has to manually transcribe the DBA converted from an LTL specification into a text file in a specific form is not convenient. To improve the usability, we are working on integrating the LTL-to-BA translation function in Spot to ROCS.



Python bindings
---------------
Having a Python interface to the C++ source code in ROCS can facilitate research in formal control synthesis. The user does not need to struggle with writing and compiling a C++ program. We are working on developing Python methods that conveniently

- interpret the DBA translated by Spot Python interface
- define system dynamics
- define labeling functions
- call the two solver engines implemented in C++
- save and display control synthesis results

In this way, writing a C++ `main` program will be replaced by writing a Python script.


Command-line tool development
-----------------------------
A command-line tool is convenient for running control synthesis for systems with simple dynamics and simple control specifications.
