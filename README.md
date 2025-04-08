# ESSDataAnalysis

This repository contains code and scripts developed as part of the feasibility study for a monitored neutrino beam at the European Spallation Source (ESS). The main goal of this project is to analyze data coming from Geant4 simulations — with a current focus on simulating the hadron dump.

## What’s in here:
- **skim_tree.cpp** — Processes the .root output from the simulation, picking out the relevant info and creating a new ROOT tree for further analysis.
- **data_analysis.cpp** — Takes the .root file produced by skim_tree.cpp and runs the actual data analysis: making plots, histograms, mapping particle positions in the x-y plane, etc.
- anything else that might come to mind (I will try to update this README whenever I add new code).

This is _very much_ a work in progress — things are likely to change and evolve.

---

## Cloning the repo

To get a local copy of this repository, simply run in your terminal:
```bash
git clone https://github.com/annascanu/ESSDataAnalysis
```

## Requirements
- [ROOT](https://root.cern/) (CERN data analysis framework): make sure ROOT is installed and properly configured on your machine.

If ROOT is available on your system, you should be able to run:
```bash
root-config --version
```

## Running the code

To compile the scripts:
```bash
c++ skim_tree.cpp `&#96;`root-config --cflags --libs`&#96;` -o skim_tree
c++ data_analysis.cpp `&#96;`root-config --cflags --libs`&#96;` -o data_analysis
```

and to run them:
```bash
./skim_tree input.root 
./data_analysis filtered_input.root
```

Replace **input.root** with the .root file produced by your Geant4 simulation.

## Notes
- This is still very early-stage work — feedback, suggestions, and contributions are welcome!
- If you run into any issues, please contact me: I have much to learn :)

## License
This project is licensed under the MIT License — see the [LICENSE](https://github.com/annascanu/ESSDataAnalysis/blob/main/LICENSE) file for details.