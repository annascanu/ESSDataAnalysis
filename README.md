# ESSDataAnalysis

This repository contains the code and scripts developed for the feasibility study of a monitored neutrino beam at the European Spallation Source (ESS). Currently, this project focuses on analyzing data derived from Geant4 simulations, specifically on the simulation of the hadron dump.

Key components of the project include:
- **skim_tree.cpp**: This script processes the .root output from the simulation, extracting relevant information and creating a new ROOT tree for further analysis.
- **analysis_muon.cpp**: This script takes the .root file produced by skim_tree.cpp as input, performing data analysis, generating plots, histograms, mapping particle positions in the x-y plane, etc.

This is _very_ much a work in progress.