#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>

static const int MAXCELLS = 10000;

void filter_tree(const char* inputFile) 
{
    // Open input file
    TFile* inFile = TFile::Open(inputFile);
    if (!inFile || inFile->IsZombie()) 
    {
        std::cout << "Can't open the input file!" << std::endl;
        return;
    }

    // Retrieve the tree from the input file
    TTree* inTree = (TTree*)inFile->Get("t");
    if (!inTree) 
    {
        std::cout << "No tree found in the input file!" << std::endl;
        inFile->Close();
        return;
    }

    std::string inputFileName = inputFile;
    // std::cout << "Splitting: " << inputFileName << '\n';
    std::size_t found_slash = inputFileName.find_last_of("/\\");
    std::string baseFileName = inputFileName.substr(found_slash + 1);

    std::string outputFileName = "filtered_" + baseFileName;

    // Create the output file
    TFile* outFile = new TFile(outputFileName.c_str(), "RECREATE");
    
    // Check if the file is created successfully
    if (!outFile || outFile->IsZombie()) 
    {
        std::cout << "Failed to create output file " << outputFileName << std::endl;
        inFile->Close();
        return;
    }

    // Declare variables to hold branch data (arrays)
    int NParticles, NPart_dep, HD_PDG[MAXCELLS], HD_zLayer[MAXCELLS], HD_zLayer_dep[MAXCELLS];
    double HD_energy_kin[MAXCELLS], HD_x[MAXCELLS], HD_y[MAXCELLS], HD_time[MAXCELLS];
    double HD_energy_dep[MAXCELLS], HD_x_dep[MAXCELLS], HD_y_dep[MAXCELLS];

    // Set branch addresses for the input tree
    inTree->SetBranchAddress("NParticles", &NParticles);
    inTree->SetBranchAddress("HD_PDG", &HD_PDG);
    inTree->SetBranchAddress("HD_x", &HD_x);
    inTree->SetBranchAddress("HD_y", &HD_y);
    inTree->SetBranchAddress("HD_zLayer", &HD_zLayer);
    inTree->SetBranchAddress("HD_time", &HD_time);
    inTree->SetBranchAddress("HD_energy_kin", &HD_energy_kin);
    // ... same thing, but for the deposited energy
    inTree->SetBranchAddress("NPart_dep", &NPart_dep);
    inTree->SetBranchAddress("HD_energy_dep", &HD_energy_dep);
    inTree->SetBranchAddress("HD_zLayer_dep", &HD_zLayer_dep);
    inTree->SetBranchAddress("HD_x_dep", &HD_x_dep);
    inTree->SetBranchAddress("HD_y_dep", &HD_y_dep);

    // Create branches in the new tree for the filtered data
    TTree* outTree = new TTree("filteredTree", "Filtered tree");
    outTree->Branch("NParticles", &NParticles, "NParticles/I");
    outTree->Branch("HD_PDG", &HD_PDG, "HD_PDG[NParticles]/I");
    outTree->Branch("HD_x", &HD_x, "HD_x[NParticles]/D");
    outTree->Branch("HD_y", &HD_y, "HD_y[NParticles]/D");
    outTree->Branch("HD_zLayer", &HD_zLayer, "HD_zLayer[NParticles]/I");
    outTree->Branch("HD_time", &HD_time, "HD_time[NParticles]/D");
    outTree->Branch("HD_energy_kin", &HD_energy_kin, "HD_energy_kin[NParticles]/D");
    // ... same thing, but for the deposited energy
    outTree->Branch("NPart_dep", &NPart_dep, "NPart_dep/I");
    outTree->Branch("HD_energy_dep", &HD_energy_dep, "HD_energy_dep[NPart_dep]/D");
    outTree->Branch("HD_zLayer_dep", &HD_zLayer_dep, "HD_zLayer_dep[NPart_dep]/D");
    outTree->Branch("HD_x_dep", &HD_x_dep, "HD_x_dep[NPart_dep]/D");
    outTree->Branch("HD_y_dep", &HD_y_dep, "HD_y_dep[NPart_dep]/D");

    Long64_t nEntries = inTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++)
    {
        inTree->GetEntry(i);
        outTree->Fill(); 
    }

    // Write the output tree to the file
    outFile->Write();
    std::cout << "Filtered tree saved in " << outputFileName << std::endl;

    // Close the files
    inFile->Close();
    outFile->Close();
}

int main(int argc, char** argv) 
{
    if (argc < 2) 
    {
        std::cout << "\n Attention: you need to insert the name of a .root file containing the tree!" << std::endl;
        return 1;
    }

    const char* input_name = argv[1];
    filter_tree(input_name);

    return 0;
}