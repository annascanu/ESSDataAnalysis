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

    // Histograms for mapping the distribution of particles on the x-y plane and x/y directions
    TH2F* hd_xy = new TH2F("hd_xy", "HD_x vs HD_y", 100, -250, 250, 100, -250, 250);
    TH1F* hd_x  = new TH1F("hd_x", "HD_x", 100, -250, 250);
    TH1F* hd_y  = new TH1F("hd_y", "HD_y", 100, -250, 250);

    TH1F* h_mu       = new TH1F("mu-", "mu-", 100, 0, 1000);
    TH1F* h_mu_plus  = new TH1F("mu+", "mu+", 100, 0, 1000);
    TH1F* h_pi       = new TH1F("pi", "pi", 100, 0, 1000);
    TH1F* h_K        = new TH1F("K", "K", 100, 0, 1000);
    TH1F* h_e        = new TH1F("e", "e", 100, 0, 1000);
    TH1F* h_n        = new TH1F("n", "n", 100, 0, 1000);
    TH1F* h_p        = new TH1F("p", "p", 100, 0, 1000);
    TH1F* h_other    = new TH1F("other", "other", 100, 0, 1000);

    Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; i++) 
    {
        t->GetEntry(i);

        hd_xy->Fill(HD_x[i], HD_y[i]);
        hd_x->Fill(HD_x[i]);
        hd_y->Fill(HD_y[i]);

        if (HD_PDG[i] == -13)            
            h_mu_plus->Fill(HD_energy_kin[i]);
        else if (HD_PDG[i] == 13)      
            h_mu->Fill(HD_energy_kin[i]);
        else if (abs(HD_PDG[i]) == 211) 
            h_pi->Fill(HD_energy_kin[i]);
        else if (abs(HD_PDG[i]) == 321) 
            h_K->Fill(HD_energy_kin[i]);
        else if (abs(HD_PDG[i]) == 11)  
            h_e->Fill(HD_energy_kin[i]);
        else if (HD_PDG[i] == 2112)     
            h_n->Fill(HD_energy_kin[i]);
        else if (abs(HD_PDG[i]) == 2212)
            h_p->Fill(HD_energy_kin[i]);
        else                         
            h_other->Fill(HD_energy_kin[i]);
    }

    TCanvas* c1 = new TCanvas("c1", "Hadron dump x, y distributions", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); hd_xy->Draw("COLZ");
    c1->cd(2); hd_x->Draw();
    c1->cd(3); hd_y->Draw();

    // Stack the kinetic energies of all the type of particles
    TCanvas* c2 = new TCanvas("c2", "Kinetic energy stack", 800, 600);
    THStack* s = new THStack("s", "Kinetic energy by particle;Kinetic energy [MeV];Number of entries");

    h_mu      ->SetFillColorAlpha(2, 0.6);
    h_mu_plus ->SetFillColorAlpha(3, 0.6);
    h_pi      ->SetFillColorAlpha(4, 0.6);
    h_K       ->SetFillColorAlpha(5, 0.6);
    h_e       ->SetFillColorAlpha(6, 0.6);
    h_n       ->SetFillColorAlpha(7, 0.6);
    h_p       ->SetFillColorAlpha(8, 0.6);
    h_other   ->SetFillColorAlpha(9, 0.6);

    s->Add(h_mu);
    s->Add(h_mu_plus);
    s->Add(h_pi);
    s->Add(h_K);
    s->Add(h_e);
    s->Add(h_n);
    s->Add(h_p);
    s->Add(h_other);
    s->Draw("hist");

    TLegend* legend = new TLegend(0.7, 0.6, 0.9, 0.88);
    legend->AddEntry(h_mu, "mu-", "f");
    legend->AddEntry(h_mu_plus, "mu+", "f");
    legend->AddEntry(h_pi, "pi", "f");
    legend->AddEntry(h_K, "K", "f");
    legend->AddEntry(h_e, "e", "f");
    legend->AddEntry(h_n, "n", "f");
    legend->AddEntry(h_p, "p", "f");
    legend->AddEntry(h_other, "other", "f");
    legend->Draw();

    myApp->Run();
    return 0;
}
