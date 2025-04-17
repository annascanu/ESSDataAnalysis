#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <cmath> // Needed for the absolute value of the PDG code

/* 
This version of data_analysis.cpp is old and I have created a new code to go back to this version
or to test some other stuff while I'm trying to code the new version.
It is fully functional, but it does not compute the rate of particles in the tree, and there is
no distinction for the zLayer variable. 
*/

static const int MAXCELLS = 10000;

int main(int argc, char** argv) 
{
    if (argc < 2)
    {
        std::cout << "\n You need to insert the name of a .root file containing the tree!" << std::endl;
        return 1;
    }

    TApplication* myApp = new TApplication("myApp", nullptr, nullptr);

    TFile* f = TFile::Open(argv[1]);
    if (!f || f->IsZombie()) 
    {
        std::cout << "Can't open root file!\n";
        return 1;
    }

    TTree* t = (TTree*)f->Get("filteredTree");
    if (!t) 
    {
        std::cout << "Did not find a root tree!\n";
        return 1;
    }

    // t->Print();
    // Declare variables for branches
    int NParticles, NPart_dep, HD_PDG[MAXCELLS], HD_zLayer[MAXCELLS], HD_zLayer_dep[MAXCELLS];
    double HD_energy_kin[MAXCELLS], HD_x[MAXCELLS], HD_y[MAXCELLS], HD_time[MAXCELLS];
    double HD_energy_dep[MAXCELLS], HD_x_dep[MAXCELLS], HD_y_dep[MAXCELLS];
    t->SetBranchAddress("NParticles", &NParticles);
    t->SetBranchAddress("HD_PDG", &HD_PDG);
    t->SetBranchAddress("HD_x", &HD_x);
    t->SetBranchAddress("HD_y", &HD_y);
    t->SetBranchAddress("HD_zLayer", &HD_zLayer);
    t->SetBranchAddress("HD_time", &HD_time);
    t->SetBranchAddress("HD_energy_kin", &HD_energy_kin);
    t->SetBranchAddress("NPart_dep", &NPart_dep);
    t->SetBranchAddress("HD_energy_dep", &HD_energy_dep);
    t->SetBranchAddress("HD_zLayer_dep", &HD_zLayer_dep);
    t->SetBranchAddress("HD_x_dep", &HD_x_dep);
    t->SetBranchAddress("HD_y_dep", &HD_y_dep);

    // Mapping the distribution of particles on the x-y plane and x/y directions
    TH2F* hd_xy = new TH2F("hd_xy", "HD_x vs HD_y", 80, -2500, 2500, 80, -2000, 2000);
    TH1F* hd_x  = new TH1F("hd_x", "HD_x", 80, -2500, 2500);
    TH1F* hd_y  = new TH1F("hd_y", "HD_y", 80, -2500, 2500);

    TH1F* h_mu       = new TH1F("mu-", "mu-", 100, 0, 1000);
    TH1F* h_mu_plus  = new TH1F("mu+", "mu+", 100, 0, 1000);
    TH1F* h_pi       = new TH1F("pi", "pi", 100, 0, 1000);
    TH1F* h_K        = new TH1F("K", "K", 100, 0, 1000);
    TH1F* h_e        = new TH1F("e", "e", 100, 0, 1000);
    TH1F* h_n        = new TH1F("n", "n", 100, 0, 1000);
    TH1F* h_p        = new TH1F("p", "p", 100, 0, 1000);
    TH1F* h_other    = new TH1F("other", "other", 100, 0, 1000);

    // Deposited energy maps
    TH2F* hd_xy_dep = new TH2F("hd_xy_dep", "HD_x_dep vs HD_y_dep", 80, -2500, 2500, 80, -2000, 2000);
    TH1F* hd_x_dep  = new TH1F("hd_x_dep", "HD_x_dep", 80, -2500, 2500);
    TH1F* hd_y_dep  = new TH1F("hd_y_dep", "HD_y_dep", 80, -2500, 2500);

    // Fill histograms by using the TTree branches directly
    // PLEASE NOTE: This method doesn’t work (it causes a bus error and crashes the TApplication), 
    // but it would definitely be a more straightforward approach. I need to understand what is 
    // causing this issue.

    /*
    t->Draw("HD_x:HD_y >> hd_xy(80,-2000,2000,80,-2000,2000)", "", "COLZ");
    t->Draw("HD_x >> hd_x(80,-2000,2000)", "", "HIST");
    t->Draw("HD_y >> hd_y(80,-2000,2000)", "", "HIST");
    t->Draw("HD_energy_kin >> h_mu(100,0,100)", "HD_PDG == 13", "HIST");
    t->Draw("HD_energy_kin >> h_mu_plus(100,0,100)", "HD_PDG == -13", "HIST");
    t->Draw("HD_energy_kin >> h_pi(100,0,100)", "abs(HD_PDG) == 211", "HIST");
    t->Draw("HD_energy_kin >> h_K(100,0,100)", "abs(HD_PDG) == 321", "HIST");
    t->Draw("HD_energy_kin >> h_e(100,0,100)", "abs(HD_PDG) == 11", "HIST");
    t->Draw("HD_energy_kin >> h_n(100,0,100)", "HD_PDG == 2112", "HIST");
    t->Draw("HD_energy_kin >> h_p(100,0,100)", "abs(HD_PDG) == 2212", "HIST");
    t->Draw("HD_energy_kin >> h_other(100,0,100)", "abs(HD_PDG) != 2212 && abs(HD_PDG) != 13 && abs(HD_PDG) != 211 && abs(HD_PDG) != 321 && abs(HD_PDG) != 11 && HD_PDG != 2112", "HIST");
    */
    
    int count_mu, count_mu_plus, count_pi, count_K, count_e, count_n, count_p, count_other;
    count_mu = count_mu_plus = count_pi = count_K = count_e = count_n = count_p = count_other = 0;

    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; i++)  // You need to loop over the events...
    {
        t->GetEntry(i); 

        // ... and then over the number of particles in each event.
        for (int j = 0; j < NParticles; j++) 
        {
            hd_xy->Fill(HD_x[j], HD_y[j]);
            hd_x->Fill(HD_x[j]);
            hd_y->Fill(HD_y[j]);

            if(HD_zLayer[j] == 2) 
            {
                if (HD_PDG[j] == 13) 
                {
                    h_mu->Fill(HD_energy_kin[j]);
                    count_mu++;
                }
                else if (HD_PDG[j] == -13)
                {
                    h_mu_plus->Fill(HD_energy_kin[j]);
                    count_mu_plus++;
                }
                else if (abs(HD_PDG[j]) == 211) 
                {
                    h_pi->Fill(HD_energy_kin[j]);
                    count_pi++;
                }
                else if (abs(HD_PDG[j]) == 321) 
                {
                    h_K->Fill(HD_energy_kin[j]);
                    count_K++;
                }
                else if (abs(HD_PDG[j]) == 11) 
                {
                    h_e->Fill(HD_energy_kin[j]);
                    count_e++;
                }
                else if (HD_PDG[j] == 2112) 
                {
                    h_n->Fill(HD_energy_kin[j]);
                    count_n++;
                }    
                else if (abs(HD_PDG[j]) == 2212) 
                {
                    h_p->Fill(HD_energy_kin[j]);
                    count_p++;
                }
                else 
                {
                    h_other->Fill(HD_energy_kin[j]);
                    count_other++;
                }
            }   

            // Fill the deposited energy histograms
            hd_xy_dep->Fill(HD_x_dep[j], HD_y_dep[j], HD_energy_dep[j]);
            hd_x_dep->Fill(HD_x_dep[j], HD_energy_dep[j]);
            hd_y_dep->Fill(HD_y_dep[j], HD_energy_dep[j]);
        }
    }

    // Canvas for x and y distributions for all events
    TCanvas* c1 = new TCanvas("c1", "Hadron dump x, y distributions", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); hd_xy->Draw("COLZ");
    c1->cd(2); hd_x->Draw();
    c1->cd(3); hd_y->Draw();

    // Stack the kinetic energies of all the types of particles
    TCanvas* c2 = new TCanvas("c2", "Kinetic energy stack", 800, 600);
    c2->SetLogy();
    THStack* stack = new THStack("stack", "Kinetic energy by particle;Kinetic energy [MeV];Number of entries");
    h_mu      -> SetFillColorAlpha(2, 0.6);
    h_mu_plus -> SetFillColorAlpha(3, 0.6);
    h_pi      -> SetFillColorAlpha(4, 0.6);
    h_K       -> SetFillColorAlpha(5, 0.6);
    h_e       -> SetFillColorAlpha(6, 0.6);
    h_n       -> SetFillColorAlpha(7, 0.6);
    h_p       -> SetFillColorAlpha(8, 0.6);
    h_other   -> SetFillColorAlpha(9, 0.6);
    stack -> Add(h_n);
    stack -> Add(h_p);
    stack -> Add(h_K);
    stack -> Add(h_e);
    stack -> Add(h_other);
    stack -> Add(h_pi);
    stack -> Add(h_mu);
    stack -> Add(h_mu_plus);
    stack -> Draw("pads");
    // Adding a legend to the stack plot
    TLegend* legend = new TLegend(0.7, 0.6, 0.9, 0.88);
    legend->AddEntry(h_mu, "mu-", "f");
    legend->AddEntry(h_mu_plus, "mu+", "f");
    legend->AddEntry(h_pi, "pi", "f");
    legend->AddEntry(h_K, "K", "f");
    legend->AddEntry(h_e, "e", "f");
    legend->AddEntry(h_n, "n", "f");
    legend->AddEntry(h_p, "p", "f");
    legend->AddEntry(h_other, "other particles", "f");
    legend->Draw();

    // Canvas for deposited energy distributions
    TCanvas* c3 = new TCanvas("c3", "Hadron dump x, y distributions", 1200, 800);
    c3->Divide(2,2);
    c3->cd(1); hd_xy_dep->Draw("COLZ");
    c3->cd(2); hd_x_dep->Draw("hist");
    c3->cd(3); hd_y_dep->Draw("hist");

    myApp->Run();
    return 0;
}