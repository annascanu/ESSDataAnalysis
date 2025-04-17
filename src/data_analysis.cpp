#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <cmath>

static const int MAXCELLS = 10000;
const int nLayers = 4;

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
    double vtxT;

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
    t->SetBranchAddress("vtxT", &vtxT); // Time for a preliminary rate estimation

    // Histograms for the x, y distributions; also doing a check on the distribution
    // on the x and y directions of muons, pions, other particles except n
    TH2F* hd_xy[nLayers];
    TH1F* hd_x[nLayers];
    TH1F* hd_y[nLayers];
    TH1F* hd_mu_plus_x[nLayers];
    TH1F* hd_mu_plus_y[nLayers];
    TH1F* hd_pi_x[nLayers];
    TH1F* hd_pi_y[nLayers];
    TH1F* hd_other_noNeutrons_x[nLayers];
    TH1F* hd_other_noNeutrons_y[nLayers];

    // Histograms for the rate of mu+ and other charged particles (pi, p, K, e)
    TH2F* h_mu_plus_rate[nLayers]; // note: This is going to be the number of particles divided by DeltaT and the bin area
    TH2F* h_otherChargedParticles_rate[nLayers];

    // Histograms for the kinetic energy: stacking the kinetic energy of each "type" of particle
    TH1F* h_mu[nLayers];
    TH1F* h_mu_plus[nLayers];
    TH1F* h_pi[nLayers];
    TH1F* h_K[nLayers];
    TH1F* h_e[nLayers];
    TH1F* h_n[nLayers];
    TH1F* h_p[nLayers];
    TH1F* h_other[nLayers];

    // Histograms for the deposited energy: finish later
    TH2F* hd_xy_dep[nLayers];
    // TH2F* hd_xy_dep = new TH2F("hd_xy_dep", "HD_x_dep vs HD_y_dep", 80, -2500, 2500, 80, -2000, 2000);
    // TH1F* hd_x_dep  = new TH1F("hd_x_dep", "HD_x_dep", 80, -2500, 2500);
    // TH1F* hd_y_dep  = new TH1F("hd_y_dep", "HD_y_dep", 80, -2500, 2500);

    for (int i = 0; i < nLayers; i++) 
    {
        hd_xy[i] = new TH2F(Form("hd_xy_z%d", i+1), Form("HD_x vs HD_y with zLayer = %d", i+1), 80, -2500, 2500, 80, -2000, 2000);
        hd_x[i]  = new TH1F(Form("hd_x_z%d", i+1), Form("HD_x with zLayer = %d", i+1), 80, -2500, 2500);
        hd_y[i]  = new TH1F(Form("hd_y_z%d", i+1), Form("HD_y with zLayer = %d", i+1), 80, -2500, 2500);
        hd_mu_plus_x[i] = new TH1F(Form("hd_mu_plus_x_z%d", i+1), Form("mu+ x distribution with zLayer = %d", i+1), 80, -2500, 2500);
        hd_mu_plus_y[i] = new TH1F(Form("hd_mu_plus_y_z%d", i+1), Form("mu+ y distribution with zLayer = %d", i+1), 80, -2500, 2500);
        hd_pi_x[i]       = new TH1F(Form("hd_pi_x_z%d", i+1), Form("pi x distribution with zLayer = %d", i+1), 80, -2500, 2500);
        hd_pi_y[i]       = new TH1F(Form("hd_pi_y_z%d", i+1), Form("pi y distribution with zLayer = %d", i+1), 80, -2500, 2500);
        hd_other_noNeutrons_x[i] = new TH1F(Form("hd_other_noNeutrons_x_z%d", i+1), Form("other x distribution with zLayer = %d", i+1), 80, -2500, 2500);
        hd_other_noNeutrons_y[i] = new TH1F(Form("hd_other_noNeutrons_y_z%d", i+1), Form("other y distribution with zLayer = %d", i+1), 80, -2500, 2500);

        h_mu_plus_rate[i] = new TH2F(Form("h_mu_plus_rate_z%d", i+1), Form("mu+ rate with zLayer = %d", i+1), 80, -2500, 2500, 80, -2000, 2000);
        h_otherChargedParticles_rate[i] = new TH2F(Form("h_otherChargedParticles_rate_z%d", i+1), Form("other charged particles rate with zLayer = %d", i+1), 80, -2500, 2500, 80, -2000, 2000);

        h_mu[i]       = new TH1F(Form("h_mu_z%d", i+1), Form("mu- with zLayer = %d", i+1), 100, 0, 1000);
        h_mu_plus[i]  = new TH1F(Form("h_mu_plus_z%d", i+1), Form("mu+ with zLayer = %d", i+1), 100, 0, 1000);
        h_pi[i]       = new TH1F(Form("h_pi_z%d", i+1), Form("pi with zLayer = %d", i+1), 100, 0, 1000);
        h_K[i]        = new TH1F(Form("h_K_z%d", i+1), Form("K with zLayer = %d", i+1), 100, 0, 1000);
        h_e[i]        = new TH1F(Form("h_e_z%d", i+1), Form("e with zLayer = %d", i+1), 100, 0, 1000);
        h_n[i]        = new TH1F(Form("h_n_z%d", i+1), Form("n with zLayer = %d", i+1), 100, 0, 1000);
        h_p[i]        = new TH1F(Form("h_p_z%d", i+1), Form("p with zLayer = %d", i+1), 100, 0, 1000);
        h_other[i]    = new TH1F(Form("h_other_z%d", i+1), Form("other with zLayer = %d", i+1), 100, 0, 1000);

        hd_xy_dep[i] = new TH2F(Form("hd_xy_dep_z%d", i+1), Form("HD_x_dep vs HD_y_dep with zLayer = %d", i+1), 80, -2500, 2500, 80, -2000, 2000);
    }

    int count_mu = 0, count_mu_plus = 0, count_pi = 0, count_K = 0, count_e = 0, count_n = 0, count_p = 0, count_other = 0;
    double max_vtxT = 0; // Needed for a preliminary rate estimation. This time basically depends on the "duration"
    // of the job; our preliminary estimate will be that the DeltaT needed for the rate calculation is equal to 10
    // times the max_vtxT, since we have 10 jobs. We will then divide the number of entries in each histogam by this
    // DeltaT and the area of each bin. (NOTE: The reason for this preliminary estimate is in my personal notes)

    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; i++)  
    {
        t->GetEntry(i);

        if (vtxT > max_vtxT)
            max_vtxT = vtxT;

        for (int j = 0; j < NParticles; j++) 
        {
            int z = HD_zLayer[j];

            if (z == 1) // First layer
            {
                int idx = z - 1;

                // x, y distributions
                hd_xy[idx]->Fill(HD_x[j], HD_y[j]);
                hd_x[idx]->Fill(HD_x[j]);
                hd_y[idx]->Fill(HD_y[j]);

                // x, y distributions for mu+, pi, and other particles
                if (HD_PDG[j] == -13)
                {
                    hd_mu_plus_x[idx]->Fill(HD_x[j]);
                    hd_mu_plus_y[idx]->Fill(HD_y[j]);
                }
                else if (abs(HD_PDG[j]) == 211)
                {
                    hd_pi_x[idx]->Fill(HD_x[j]);
                    hd_pi_y[idx]->Fill(HD_y[j]);
                }
                else if (abs(HD_PDG[j]) != 2112)
                {
                    hd_other_noNeutrons_x[idx]->Fill(HD_x[j]);
                    hd_other_noNeutrons_y[idx]->Fill(HD_y[j]);
                }

                // Fill histograms for the rate of mu+ and other charged particles
                if (HD_PDG[j] == -13)
                    h_mu_plus_rate[idx]->Fill(HD_x[j], HD_y[j]);   
                else if (abs(HD_PDG[j]) == 211 || abs(HD_PDG[j]) == 2212 || abs(HD_PDG[j]) == 321 || abs(HD_PDG[j]) == 11 || HD_PDG[j] == 13) // pi, p, K, e, mu-
                    h_otherChargedParticles_rate[idx]->Fill(HD_x[j], HD_y[j]);

                // Fill histograms for the kinetic energy, for each particle to stack later
                if (HD_PDG[j] == 13 && HD_energy_kin[j] > 10) 
                {
                    h_mu[idx]->Fill(HD_energy_kin[j]);
                    count_mu++;
                }
                else if (HD_PDG[j] == -13  && HD_energy_kin[j] > 10)
                {
                    h_mu_plus[idx]->Fill(HD_energy_kin[j]);
                    count_mu_plus++;
                }
                else if (abs(HD_PDG[j]) == 211  && HD_energy_kin[j] > 10) 
                {
                    h_pi[idx]->Fill(HD_energy_kin[j]);
                    count_pi++;
                }
                else if (abs(HD_PDG[j]) == 321  && HD_energy_kin[j] > 10) 
                {
                    h_K[idx]->Fill(HD_energy_kin[j]);
                    count_K++;
                }
                else if (abs(HD_PDG[j]) == 11  && HD_energy_kin[j] > 10) 
                {
                    h_e[idx]->Fill(HD_energy_kin[j]);
                    count_e++;
                }
                else if (HD_PDG[j] == 2112  && HD_energy_kin[j] > 10) 
                {
                    h_n[idx]->Fill(HD_energy_kin[j]);
                    count_n++;
                }    
                else if (abs(HD_PDG[j]) == 2212 && HD_energy_kin[j] > 10) 
                {
                    h_p[idx]->Fill(HD_energy_kin[j]);
                    count_p++;
                }
                else if (HD_energy_kin[j] > 10)
                {
                    h_other[idx]->Fill(HD_energy_kin[j]);
                    count_other++;
                }

                // Fill the deposited energy histograms
                hd_xy_dep[idx]->Fill(HD_x_dep[j], HD_y_dep[j], HD_energy_dep[j]);
            }
            else if (z == 2) // Second layer
            {
                int idx = z - 1;

                // x, y distributions
                hd_xy[idx]->Fill(HD_x[j], HD_y[j]);
                hd_x[idx]->Fill(HD_x[j]);
                hd_y[idx]->Fill(HD_y[j]);

                // x, y distributions for mu+, pi, and other particles
                if (HD_PDG[j] == -13)
                {
                    hd_mu_plus_x[idx]->Fill(HD_x[j]);
                    hd_mu_plus_y[idx]->Fill(HD_y[j]);
                }
                else if (abs(HD_PDG[j]) == 211)
                {
                    hd_pi_x[idx]->Fill(HD_x[j]);
                    hd_pi_y[idx]->Fill(HD_y[j]);
                }
                else if (abs(HD_PDG[j]) != 2112)
                {
                    hd_other_noNeutrons_x[idx]->Fill(HD_x[j]);
                    hd_other_noNeutrons_y[idx]->Fill(HD_y[j]);
                }

                // Fill histograms for the rate of mu+ and other charged particles
                if (HD_PDG[j] == -13)
                    h_mu_plus_rate[idx]->Fill(HD_x[j], HD_y[j]);   
                else if (abs(HD_PDG[j]) == 211 || abs(HD_PDG[j]) == 2212 || abs(HD_PDG[j]) == 321 || abs(HD_PDG[j]) == 11 || HD_PDG[j] == 13) // pi, p, K, e, mu-
                    h_otherChargedParticles_rate[idx]->Fill(HD_x[j], HD_y[j]);

                // Fill histograms for the kinetic energy, for each particle to stack later
                if (HD_PDG[j] == 13 && HD_energy_kin[j] > 10) 
                {
                    h_mu[idx]->Fill(HD_energy_kin[j]);
                    count_mu++;
                }
                else if (HD_PDG[j] == -13  && HD_energy_kin[j] > 10)
                {
                    h_mu_plus[idx]->Fill(HD_energy_kin[j]);
                    count_mu_plus++;
                }
                else if (abs(HD_PDG[j]) == 211  && HD_energy_kin[j] > 10) 
                {
                    h_pi[idx]->Fill(HD_energy_kin[j]);
                    count_pi++;
                }
                else if (abs(HD_PDG[j]) == 321  && HD_energy_kin[j] > 10) 
                {
                    h_K[idx]->Fill(HD_energy_kin[j]);
                    count_K++;
                }
                else if (abs(HD_PDG[j]) == 11  && HD_energy_kin[j] > 10) 
                {
                    h_e[idx]->Fill(HD_energy_kin[j]);
                    count_e++;
                }
                else if (HD_PDG[j] == 2112  && HD_energy_kin[j] > 10) 
                {
                    h_n[idx]->Fill(HD_energy_kin[j]);
                    count_n++;
                }    
                else if (abs(HD_PDG[j]) == 2212 && HD_energy_kin[j] > 10) 
                {
                    h_p[idx]->Fill(HD_energy_kin[j]);
                    count_p++;
                }
                else if (HD_energy_kin[j] > 10)
                {
                    h_other[idx]->Fill(HD_energy_kin[j]);
                    count_other++;
                }

                // Fill the deposited energy histograms
                hd_xy_dep[idx]->Fill(HD_x_dep[j], HD_y_dep[j], HD_energy_dep[j]);
            }
            else if (z == 3) // Third layer
            {
                int idx = z - 1;

                // x, y distributions
                hd_xy[idx]->Fill(HD_x[j], HD_y[j]);
                hd_x[idx]->Fill(HD_x[j]);
                hd_y[idx]->Fill(HD_y[j]);

                // x, y distributions for mu+, pi, and other particles
                if (HD_PDG[j] == -13)
                {
                    hd_mu_plus_x[idx]->Fill(HD_x[j]);
                    hd_mu_plus_y[idx]->Fill(HD_y[j]);
                }
                else if (abs(HD_PDG[j]) == 211)
                {
                    hd_pi_x[idx]->Fill(HD_x[j]);
                    hd_pi_y[idx]->Fill(HD_y[j]);
                }
                else if (abs(HD_PDG[j]) != 2112)
                {
                    hd_other_noNeutrons_x[idx]->Fill(HD_x[j]);
                    hd_other_noNeutrons_y[idx]->Fill(HD_y[j]);
                }

                // Fill histograms for the rate of mu+ and other charged particles
                if (HD_PDG[j] == -13)
                    h_mu_plus_rate[idx]->Fill(HD_x[j], HD_y[j]);   
                else if (abs(HD_PDG[j]) == 211 || abs(HD_PDG[j]) == 2212 || abs(HD_PDG[j]) == 321 || abs(HD_PDG[j]) == 11 || HD_PDG[j] == 13) // pi, p, K, e, mu-
                    h_otherChargedParticles_rate[idx]->Fill(HD_x[j], HD_y[j]);

                // Fill histograms for the kinetic energy, for each particle to stack later
                if (HD_PDG[j] == 13 && HD_energy_kin[j] > 10) 
                {
                    h_mu[idx]->Fill(HD_energy_kin[j]);
                    count_mu++;
                }
                else if (HD_PDG[j] == -13  && HD_energy_kin[j] > 10)
                {
                    h_mu_plus[idx]->Fill(HD_energy_kin[j]);
                    count_mu_plus++;
                }
                else if (abs(HD_PDG[j]) == 211  && HD_energy_kin[j] > 10) 
                {
                    h_pi[idx]->Fill(HD_energy_kin[j]);
                    count_pi++;
                }
                else if (abs(HD_PDG[j]) == 321  && HD_energy_kin[j] > 10) 
                {
                    h_K[idx]->Fill(HD_energy_kin[j]);
                    count_K++;
                }
                else if (abs(HD_PDG[j]) == 11  && HD_energy_kin[j] > 10) 
                {
                    h_e[idx]->Fill(HD_energy_kin[j]);
                    count_e++;
                }
                else if (HD_PDG[j] == 2112  && HD_energy_kin[j] > 10) 
                {
                    h_n[idx]->Fill(HD_energy_kin[j]);
                    count_n++;
                }    
                else if (abs(HD_PDG[j]) == 2212 && HD_energy_kin[j] > 10) 
                {
                    h_p[idx]->Fill(HD_energy_kin[j]);
                    count_p++;
                }
                else if (HD_energy_kin[j] > 10)
                {
                    h_other[idx]->Fill(HD_energy_kin[j]);
                    count_other++;
                }

                // Fill the deposited energy histograms
                hd_xy_dep[idx]->Fill(HD_x_dep[j], HD_y_dep[j], HD_energy_dep[j]);
            }
            else if (z == 4) // Fourth layer
            {
                int idx = z - 1;

                // x, y distributions
                hd_xy[idx]->Fill(HD_x[j], HD_y[j]);
                hd_x[idx]->Fill(HD_x[j]);
                hd_y[idx]->Fill(HD_y[j]);

                // x, y distributions for mu+, pi, and other particles
                if (HD_PDG[j] == -13)
                {
                    hd_mu_plus_x[idx]->Fill(HD_x[j]);
                    hd_mu_plus_y[idx]->Fill(HD_y[j]);
                }
                else if (abs(HD_PDG[j]) == 211)
                {
                    hd_pi_x[idx]->Fill(HD_x[j]);
                    hd_pi_y[idx]->Fill(HD_y[j]);
                }
                else if (abs(HD_PDG[j]) != 2112)
                {
                    hd_other_noNeutrons_x[idx]->Fill(HD_x[j]);
                    hd_other_noNeutrons_y[idx]->Fill(HD_y[j]);
                }

                // Fill histograms for the rate of mu+ and other charged particles
                if (HD_PDG[j] == -13)
                    h_mu_plus_rate[idx]->Fill(HD_x[j], HD_y[j]);   
                else if (abs(HD_PDG[j]) == 211 || abs(HD_PDG[j]) == 2212 || abs(HD_PDG[j]) == 321 || abs(HD_PDG[j]) == 11 || HD_PDG[j] == 13) // pi, p, K, e, mu-
                    h_otherChargedParticles_rate[idx]->Fill(HD_x[j], HD_y[j]);

                // Fill histograms for the kinetic energy, for each particle to stack later
                if (HD_PDG[j] == 13 && HD_energy_kin[j] > 10) 
                {
                    h_mu[idx]->Fill(HD_energy_kin[j]);
                    count_mu++;
                }
                else if (HD_PDG[j] == -13  && HD_energy_kin[j] > 10)
                {
                    h_mu_plus[idx]->Fill(HD_energy_kin[j]);
                    count_mu_plus++;
                }
                else if (abs(HD_PDG[j]) == 211  && HD_energy_kin[j] > 10) 
                {
                    h_pi[idx]->Fill(HD_energy_kin[j]);
                    count_pi++;
                }
                else if (abs(HD_PDG[j]) == 321  && HD_energy_kin[j] > 10) 
                {
                    h_K[idx]->Fill(HD_energy_kin[j]);
                    count_K++;
                }
                else if (abs(HD_PDG[j]) == 11  && HD_energy_kin[j] > 10) 
                {
                    h_e[idx]->Fill(HD_energy_kin[j]);
                    count_e++;
                }
                else if (HD_PDG[j] == 2112  && HD_energy_kin[j] > 10) 
                {
                    h_n[idx]->Fill(HD_energy_kin[j]);
                    count_n++;
                }    
                else if (abs(HD_PDG[j]) == 2212 && HD_energy_kin[j] > 10) 
                {
                    h_p[idx]->Fill(HD_energy_kin[j]);
                    count_p++;
                }
                else if (HD_energy_kin[j] > 10)
                {
                    h_other[idx]->Fill(HD_energy_kin[j]);
                    count_other++;
                }

                // Fill the deposited energy histograms
                hd_xy_dep[idx]->Fill(HD_x_dep[j], HD_y_dep[j], HD_energy_dep[j]);
            }   
            
            // hd_xy_dep->Fill(HD_x_dep[j], HD_y_dep[j], HD_energy_dep[j]);
            // hd_x_dep->Fill(HD_x_dep[j], HD_energy_dep[j]);
            // hd_y_dep->Fill(HD_y_dep[j], HD_energy_dep[j]);
        }
    }

    double Delta_T = max_vtxT * 10; // Preliminary estimate for the DeltaT needed for the rate calculation: it is 
    // equal to 10 times the duration of each job, assuming that the number of jobs is equal to 10.
    // I might need to make the number of jobs a variable in the future.
    double area = 25.0; // Area of each bin in the x-y plane (in cm^2)
    double rate_conversion_factor = 1 / (Delta_T * area); // Conversion factor for the rate calculation

    // ------------------------------------------------------------------------------------
    // x, y, x vs. y distributions; also dividing by number of layer and type of particle
    //     just for the case of muons, pions and other particles (except neutrons)
    // ------------------------------------------------------------------------------------

    TCanvas* c1 = new TCanvas("c1", "Hadron dump x, y distributions", 1200, 800);
    c1 -> Print("Plots_ESS_HDAnalysis.pdf[");

    TCanvas* canvas_x_y_distrib[nLayers];
    TCanvas* canvas_mu_distrib_x[nLayers];
    TLegend* legend_mu_distrib_x[nLayers]; 
    TCanvas* canvas_mu_distrib_y[nLayers];
    TLegend* legend_mu_distrib_y[nLayers];
    TCanvas* canvas_rate[nLayers];
    /// TLegend* legend_rate[nLayers];

    for (int i = 0; i < nLayers; i++)
    {
        legend_mu_distrib_x[i] = new TLegend(0.7, 0.6, 0.9, 0.88);
        legend_mu_distrib_y[i] = new TLegend(0.7, 0.6, 0.9, 0.88);

        canvas_x_y_distrib[i] = new TCanvas(Form("canvas_x_y_distrib%d", i+1), Form("Hadron dump x, y distributions for zLayer = %d", i+1), 1200, 800);
        canvas_x_y_distrib[i] -> cd();
        canvas_x_y_distrib[i] -> Divide(2, 2);
        canvas_x_y_distrib[i] -> cd(1);
        hd_xy[i] -> Draw("COLZ");
        canvas_x_y_distrib[i] -> cd(2);
        hd_x[i] -> Draw("hist");
        canvas_x_y_distrib[i] -> cd(3);
        hd_y[i] -> Draw("hist");
        canvas_x_y_distrib[i] -> Print("Plots_ESS_HDAnalysis.pdf");

        // x distribution for mu+, pi, and other particles
        hd_mu_plus_x[i]->SetLineColor(kRed);
        hd_mu_plus_x[i]->SetLineWidth(2);
        hd_pi_x[i]->SetLineColor(kBlue);
        hd_pi_x[i]->SetLineWidth(2);
        hd_other_noNeutrons_x[i]->SetLineColor(kGreen+2);
        hd_other_noNeutrons_x[i]->SetLineWidth(2);

        canvas_mu_distrib_x[i] = new TCanvas(Form("canvas_mu_distrib_x%d", i+1), Form("Hadron dump x distribution for mu+ in zLayer = %d", i+1), 1200, 800);
        canvas_mu_distrib_x[i]->cd();
        hd_mu_plus_x[i]->Draw("hist");
        hd_pi_x[i]->Draw("hist same");
        hd_other_noNeutrons_x[i]->Draw("hist same");
        legend_mu_distrib_x[i]->AddEntry(hd_mu_plus_x[i], "mu+", "f");
        legend_mu_distrib_x[i]->AddEntry(hd_pi_x[i], "pi", "f");
        legend_mu_distrib_x[i]->AddEntry(hd_other_noNeutrons_x[i], "other particles (no neutrons)", "f");
        legend_mu_distrib_x[i]->Draw();
        canvas_mu_distrib_x[i]->Print("Plots_ESS_HDAnalysis.pdf");

        // y distribution for mu+, pi, and other particles
        hd_mu_plus_y[i]->SetLineColor(kRed);
        hd_mu_plus_y[i]->SetLineWidth(2);
        hd_pi_y[i]->SetLineColor(kBlue);
        hd_pi_y[i]->SetLineWidth(2);
        hd_other_noNeutrons_y[i]->SetLineColor(kGreen+2);
        hd_other_noNeutrons_y[i]->SetLineWidth(2);

        canvas_mu_distrib_y[i] = new TCanvas(Form("canvas_mu_distrib_y%d", i+1), Form("Hadron dump y distribution for mu+ in zLayer = %d", i+1), 1200, 800);
        canvas_mu_distrib_y[i]->cd();
        hd_mu_plus_y[i]->Draw("hist");
        hd_pi_y[i]->Draw("hist same");
        hd_other_noNeutrons_y[i]->Draw("hist same");
        legend_mu_distrib_y[i]->AddEntry(hd_mu_plus_y[i], "mu+", "f");
        legend_mu_distrib_y[i]->AddEntry(hd_pi_y[i], "pi", "f");
        legend_mu_distrib_y[i]->AddEntry(hd_other_noNeutrons_y[i], "other particles (no neutrons)", "f");
        legend_mu_distrib_y[i]->Draw();
        canvas_mu_distrib_y[i]->Print("Plots_ESS_HDAnalysis.pdf");

        // Another canvas for the rate
        canvas_rate[i] = new TCanvas(Form("canvas_rate%d", i+1), Form("Rate of mu+ and other charged particles in zLayer = %d", i+1), 1200, 800);
        canvas_rate[i]->cd();
        canvas_rate[i] -> Divide(2, 1);
        h_mu_plus_rate[i]->Scale(rate_conversion_factor);
        h_mu_plus_rate[i]->SetLineColor(kRed);
        h_mu_plus_rate[i]->SetLineWidth(2);
        h_otherChargedParticles_rate[i]->Scale(rate_conversion_factor);
        h_otherChargedParticles_rate[i]->SetLineColor(kBlue);
        h_otherChargedParticles_rate[i]->SetLineWidth(2);

        canvas_rate[i] -> cd(1);
        h_mu_plus_rate[i]->SetTitle(Form("Rate of mu+ in zLayer = %d", i+1));
        h_mu_plus_rate[i]->GetXaxis()->SetTitle("x (cm)");
        h_mu_plus_rate[i]->GetYaxis()->SetTitle("y (cm)");
        h_mu_plus_rate[i]->GetZaxis()->SetTitle("Rate (Hz/cm^2)");
        h_mu_plus_rate[i]->Draw("COLZ");
        canvas_rate[i] -> cd(2);
        h_otherChargedParticles_rate[i]->SetTitle(Form("Rate of other charged particles in zLayer = %d", i+1));
        h_otherChargedParticles_rate[i]->GetXaxis()->SetTitle("x (mm)");
        h_otherChargedParticles_rate[i]->GetYaxis()->SetTitle("y (mm)");
        h_otherChargedParticles_rate[i]->GetZaxis()->SetTitle("Rate (Hz/cm^2)");
        h_otherChargedParticles_rate[i]->Draw("COLZ");
        canvas_rate[i]->Print("Plots_ESS_HDAnalysis.pdf");
    }

    // ------------------------------------------------------------------------------------
    //        Kinetic energy histograms: divided by type of particle and zLayer
    // ------------------------------------------------------------------------------------

    TCanvas* canvas_stack[nLayers];
    THStack* stack[nLayers];
    TLegend* legend[nLayers];

    // Creating the number of canvases, stacks and legends for each layer
    for (int i = 0; i < nLayers; i++) 
    {
        canvas_stack[i] = new TCanvas(Form("canvas_kinEnergy%d", i+1), Form("Kinetic energy by particle in zLayer = %d", i+1), 800, 600);
        stack[i] = new THStack(Form("stack_z%d", i+1), Form("Kinetic energy by particle in zLayer %d", i+1));
        legend[i] = new TLegend(0.7, 0.6, 0.9, 0.88);
    }

    for (int i = 0; i < nLayers; i++)
    {
        canvas_stack[i]->cd();
        canvas_stack[i]->SetLogy();

        h_mu[i] -> SetFillColorAlpha(2, 0.6);
        h_mu_plus[i] -> SetFillColorAlpha(3, 0.6);
        h_pi[i] -> SetFillColorAlpha(4, 0.6);
        h_K[i] -> SetFillColorAlpha(5, 0.6);
        h_e[i] -> SetFillColorAlpha(6, 0.6);
        h_n[i] -> SetFillColorAlpha(7, 0.6);
        h_p[i] -> SetFillColorAlpha(8, 0.6);
        h_other[i] -> SetFillColorAlpha(9, 0.6);
        stack[i] -> Add(h_n[i]);
        stack[i] -> Add(h_p[i]);
        stack[i] -> Add(h_K[i]);
        stack[i] -> Add(h_e[i]);
        stack[i] -> Add(h_other[i]);
        stack[i] -> Add(h_pi[i]);
        stack[i] -> Add(h_mu[i]);
        stack[i] -> Add(h_mu_plus[i]);
        stack[i] ->Draw();

        legend[i]->AddEntry(h_mu[i], "mu-", "f");
        legend[i]->AddEntry(h_mu_plus[i], "mu+", "f");
        legend[i]->AddEntry(h_pi[i], "pi", "f");
        legend[i]->AddEntry(h_K[i], "K", "f");
        legend[i]->AddEntry(h_e[i], "e", "f");
        legend[i]->AddEntry(h_n[i], "n", "f");
        legend[i]->AddEntry(h_p[i], "p", "f");
        legend[i]->AddEntry(h_other[i], "Other particles", "f");
        legend[i]->Draw();

        canvas_stack[i]->Print("Plots_ESS_HDAnalysis.pdf");
    }

    /*
    TCanvas* c3 = new TCanvas("c3", "Distribution of the dose on the x-y plane", 1200, 800);
    c3->Divide(2,2);
    //hd_xy_dep->Scale(3e1022/1e9); 
    c3->cd(1); 
    hd_xy_dep->Draw("COLZ");
    c3->cd(2); 
    hd_x_dep->Draw("hist");
    c3->cd(3); 
    hd_y_dep->Draw("hist");
    */

    TCanvas* c4 = new TCanvas("c4", "Hadron dump x, y distributions", 1200, 800);
    c4->Print("Plots_ESS_HDAnalysis.pdf]");

    myApp->Run();
    return 0;
}