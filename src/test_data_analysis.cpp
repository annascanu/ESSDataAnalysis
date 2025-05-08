#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
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

    // Variables for the inizialization of histograms
    int x_max = 2500;      // In mm 
    int x_min = -2500;
    int y_max = 2500;
    int y_min = -2500;
    int n_bins = (x_max-x_min) / 50;    // 5 cm wide bins

    for (int i = 0; i < nLayers; i++) 
    {
        hd_xy[i] = new TH2F(Form("hd_xy_z%d", i+1), Form("HD_x vs HD_y with zLayer = %d", i+1), n_bins, -2500, 2500, n_bins, -2500, 2500);
        hd_x[i]  = new TH1F(Form("hd_x_z%d", i+1), Form("HD_x with zLayer = %d", i+1), n_bins, -2500, 2500);
        hd_y[i]  = new TH1F(Form("hd_y_z%d", i+1), Form("HD_y with zLayer = %d", i+1), n_bins, -2500, 2500);
        hd_mu_plus_x[i] = new TH1F(Form("hd_mu_plus_x_z%d", i+1), Form("mu+ x distribution with zLayer = %d", i+1), n_bins, -2500, 2500);
        hd_mu_plus_y[i] = new TH1F(Form("hd_mu_plus_y_z%d", i+1), Form("mu+ y distribution with zLayer = %d", i+1), n_bins, -2500, 2500);
        hd_pi_x[i]       = new TH1F(Form("hd_pi_x_z%d", i+1), Form("pi x distribution with zLayer = %d", i+1), n_bins, -2500, 2500);
        hd_pi_y[i]       = new TH1F(Form("hd_pi_y_z%d", i+1), Form("pi y distribution with zLayer = %d", i+1), n_bins, -2500, 2500);
        hd_other_noNeutrons_x[i] = new TH1F(Form("hd_other_noNeutrons_x_z%d", i+1), Form("other x distribution with zLayer = %d", i+1), n_bins, -2500, 2500);
        hd_other_noNeutrons_y[i] = new TH1F(Form("hd_other_noNeutrons_y_z%d", i+1), Form("other y distribution with zLayer = %d", i+1), n_bins, -2500, 2500);

        h_mu_plus_rate[i] = new TH2F(Form("h_mu_plus_rate_z%d", i+1), Form("mu+ rate with zLayer = %d", i+1), n_bins, -2500, 2500, n_bins, -2500, 2500);
        h_otherChargedParticles_rate[i] = new TH2F(Form("h_otherChargedParticles_rate_z%d", i+1), Form("other charged particles rate with zLayer = %d", i+1), n_bins, -2500, 2500, n_bins, -2500, 2500);

        h_mu[i]       = new TH1F(Form("h_mu_z%d", i+1), Form("mu- with zLayer = %d", i+1), 100, 0, 1000);
        h_mu_plus[i]  = new TH1F(Form("h_mu_plus_z%d", i+1), Form("mu+ with zLayer = %d", i+1), 100, 0, 1000);
        h_pi[i]       = new TH1F(Form("h_pi_z%d", i+1), Form("pi with zLayer = %d", i+1), 100, 0, 1000);
        h_K[i]        = new TH1F(Form("h_K_z%d", i+1), Form("K with zLayer = %d", i+1), 100, 0, 1000);
        h_e[i]        = new TH1F(Form("h_e_z%d", i+1), Form("e with zLayer = %d", i+1), 100, 0, 1000);
        h_n[i]        = new TH1F(Form("h_n_z%d", i+1), Form("n with zLayer = %d", i+1), 100, 0, 1000);
        h_p[i]        = new TH1F(Form("h_p_z%d", i+1), Form("p with zLayer = %d", i+1), 100, 0, 1000);
        h_other[i]    = new TH1F(Form("h_other_z%d", i+1), Form("other with zLayer = %d", i+1), 100, 0, 1000);

        hd_xy_dep[i] = new TH2F(Form("hd_xy_dep_z%d", i+1), Form("HD_x_dep vs HD_y_dep with zLayer = %d", i+1), n_bins, -2500, 2500, n_bins, -2500, 2500);
    }

    int count_mu = 0, count_mu_plus = 0, count_pi = 0, count_K = 0, count_e = 0, count_n = 0, count_p = 0, count_other = 0;
    double max_vtxT = 0; // Needed for a preliminary rate estimation. This time basically depends on the "duration"
    // of the job; our preliminary estimate will be that the DeltaT needed for the rate calculation is equal to 10
    // times the max_vtxT, since we have 10 jobs. We will then divide the number of entries in each histogam by this
    // DeltaT and the area of each bin. (NOTE: The reason for this preliminary estimate is in my personal notes)

    Long64_t nentries = t->GetEntries();
for (Long64_t i = 0; i < nentries; i++) {
    t->GetEntry(i);

    if (vtxT > max_vtxT)
        max_vtxT = vtxT;

    for (int j = 0; j < NParticles; j++) {
        int z = HD_zLayer[j];
        if (z < 1 || z > 4) continue; // skip invalid layers
        int idx = z - 1;

        // x, y distributions
        hd_xy[idx]->Fill(HD_x[j], HD_y[j]);
        hd_x[idx]->Fill(HD_x[j]);
        hd_y[idx]->Fill(HD_y[j]);

        // Particle-specific x/y histograms
        if (HD_PDG[j] == -13) {
            hd_mu_plus_x[idx]->Fill(HD_x[j]);
            hd_mu_plus_y[idx]->Fill(HD_y[j]);
        } else if (abs(HD_PDG[j]) == 211) {
            hd_pi_x[idx]->Fill(HD_x[j]);
            hd_pi_y[idx]->Fill(HD_y[j]);
        } else if (abs(HD_PDG[j]) != 2112) {
            hd_other_noNeutrons_x[idx]->Fill(HD_x[j]);
            hd_other_noNeutrons_y[idx]->Fill(HD_y[j]);
        }

        // Rate histograms
        if (HD_PDG[j] == -13) {
            h_mu_plus_rate[idx]->Fill(HD_x[j], HD_y[j]);   
        } else if (abs(HD_PDG[j]) == 211 || abs(HD_PDG[j]) == 2212 ||
                   abs(HD_PDG[j]) == 321 || abs(HD_PDG[j]) == 11 || HD_PDG[j] == 13) {
            h_otherChargedParticles_rate[idx]->Fill(HD_x[j], HD_y[j]);
        }

        // Kinetic energy histograms (energy threshold > 10)
        double E_kin = HD_energy_kin[j];
        if (E_kin > 10) {
            int pdg = HD_PDG[j];
            if (pdg == 13) {
                h_mu[idx]->Fill(E_kin);
                count_mu++;
            } else if (pdg == -13) {
                h_mu_plus[idx]->Fill(E_kin);
                count_mu_plus++;
            } else if (abs(pdg) == 211) {
                h_pi[idx]->Fill(E_kin);
                count_pi++;
            } else if (abs(pdg) == 321) {
                h_K[idx]->Fill(E_kin);
                count_K++;
            } else if (abs(pdg) == 11) {
                h_e[idx]->Fill(E_kin);
                count_e++;
            } else if (pdg == 2112) {
                h_n[idx]->Fill(E_kin);
                count_n++;
            } else if (abs(pdg) == 2212) {
                h_p[idx]->Fill(E_kin);
                count_p++;
            } else {
                h_other[idx]->Fill(E_kin);
                count_other++;
            }
        }

        // Deposited energy
        hd_xy_dep[idx]->Fill(HD_x_dep[j], HD_y_dep[j], HD_energy_dep[j]);
    }
}

    double Delta_T = max_vtxT * 10; // Preliminary estimate for the DeltaT needed for the rate calculation: it is 
    // equal to 10 times the duration of each job, assuming that the number of jobs is equal to 10.
    // I might need to make the number of jobs a variable in the future.
    double area = 50 * 50; // Area of each bin in the x-y plane (in mm^2)
    double rate_conversion_factor = 1 / (Delta_T * area); // Conversion factor for the rate calculation

    // ------------------------------------------------------------------------------------
    // x, y, x vs. y distributions; also dividing by number of layer and type of particle
    //     just for the case of muons, pions and other particles (except neutrons)
    // ------------------------------------------------------------------------------------

    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "Hadron dump x, y distributions", 1200, 800);
    c1 -> Print("Plots_ESS_HDAnalysis.pdf[");

    TCanvas* canvas_x_y_distrib[nLayers];
    TCanvas* canvas_mu_distrib_x[nLayers];
    TLegend* legend_mu_distrib_x[nLayers]; 
    TCanvas* canvas_mu_distrib_y[nLayers];
    TLegend* legend_mu_distrib_y[nLayers];
    TCanvas* canvas_rate[nLayers];
    THStack* stack_x_distr[nLayers];
    THStack* stack_y_distr[nLayers];
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
        //hd_mu_plus_x[i]->Draw("hist");
        //hd_pi_x[i]->Draw("hist same");
        //hd_other_noNeutrons_x[i]->Draw("hist same");
        stack_x_distr[i] = new THStack(Form("stack_x_distr_z%d", i+1), Form("x distribution for mu+, pi, and other particles in zLayer = %d", i+1));
        stack_x_distr[i]->Add(hd_mu_plus_x[i]);
        stack_x_distr[i]->Add(hd_pi_x[i]);
        stack_x_distr[i]->Add(hd_other_noNeutrons_x[i]);
        legend_mu_distrib_x[i]->AddEntry(hd_mu_plus_x[i], "mu+", "f");
        legend_mu_distrib_x[i]->AddEntry(hd_pi_x[i], "pi", "f");
        legend_mu_distrib_x[i]->AddEntry(hd_other_noNeutrons_x[i], "other particles (no neutrons)", "f");
        stack_x_distr[i]->Draw("nostack");
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
        //hd_mu_plus_y[i]->Draw("hist");
        //hd_pi_y[i]->Draw("hist same");
        //hd_other_noNeutrons_y[i]->Draw("hist same");
        stack_y_distr[i] = new THStack(Form("stack_y_distr_z%d", i+1), Form("y distribution for mu+, pi, and other particles in zLayer = %d", i+1));
        stack_y_distr[i]->Add(hd_mu_plus_y[i]);
        stack_y_distr[i]->Add(hd_pi_y[i]);
        stack_y_distr[i]->Add(hd_other_noNeutrons_y[i]);
        legend_mu_distrib_y[i]->AddEntry(hd_mu_plus_y[i], "mu+", "f");
        legend_mu_distrib_y[i]->AddEntry(hd_pi_y[i], "pi", "f");
        legend_mu_distrib_y[i]->AddEntry(hd_other_noNeutrons_y[i], "other particles (no neutrons)", "f");
        stack_y_distr[i]->Draw("nostack");
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
        h_mu_plus_rate[i]->GetXaxis()->SetTitle("x (mm)");
        h_mu_plus_rate[i]->GetYaxis()->SetTitle("y (mm)");
        h_mu_plus_rate[i]->GetZaxis()->SetTitle("Rate (Hz/mm^2)");
        h_mu_plus_rate[i]->Draw("COLZ");
        canvas_rate[i] -> cd(2);
        h_otherChargedParticles_rate[i]->SetTitle(Form("Rate of other charged particles in zLayer = %d", i+1));
        h_otherChargedParticles_rate[i]->GetXaxis()->SetTitle("x (mm)");
        h_otherChargedParticles_rate[i]->GetYaxis()->SetTitle("y (mm)");
        h_otherChargedParticles_rate[i]->GetZaxis()->SetTitle("Rate (Hz/mm^2)");
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
        canvas_stack[i] = new TCanvas(Form("canvas_kinEnergy%d", i+1), Form("Kinetic energy by particle in zLayer = %d", i+1), 1000, 800);
        stack[i] = new THStack(Form("stack_z%d", i+1), Form("Kinetic energy by particle in zLayer %d", i+1));
        legend[i] = new TLegend(0.7, 0.7, 0.8, 0.89);
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

    // ------------------------------------------------------------------------------------
    //                Deposited energy histograms: divided by zLayer
    // ------------------------------------------------------------------------------------

    double bin_volume = 1 * 5 * 5; // in cm^3 
    double bin_mass = 0;

    double density_Ne = 0.839e-3; // g/cm^3
    double density_Eth = 1.263e-3; // g/cm^3
    double density_CF4 = 1.25e-3; // g/cm^3
    double mix_Ne = 0.8; // 80% of Ne
    double mix_Eth = 0.1; // 10% of Eth
    double mix_CF4 = 0.1; // 10% of CF4

    bin_mass = (bin_volume*(density_Ne * mix_Ne + density_Eth * mix_Eth + density_CF4 * mix_CF4))/1e3; // In kg
    // Conversion factor: multiply by conversion from MeV to J and divide by the mass of the bin
    // All of this is also normalized by the total number of POTs in a year, which is 3e22.

    // Later on: might make this conversion factor not hard-coded so that it depends on the actual number
    // of jobs and primaries per job. For now: 10 jobs of 10000 pi+ each = 10^5 primaries
    double number_POT_year = 3e22;
    // int number_jobs = 10;
    // int number_piPlus = 10000;
    // int number_total = number_jobs * number_primaries;
    // int number_POTs_jobs = number_total * 1e4; 
    // WHY THIS CONVERSION FACTOR? Since for each POT we have about 1e-4 pi+ entering
    // the tunnel, to obtain the number of POTs needed for the number of primaries that we set we need to
    // multiply the total number of primaries simulated at the entrance of the tunnel by a factor 1e4.
    // int scale = number_POT_year / number_POT_jobs;

    TCanvas* canvas_dep_energy[nLayers];

    for (int i = 0; i < nLayers; i++)
    {
        canvas_dep_energy[i] = new TCanvas(Form("canvas_dep_energy%d", i+1), Form("Dose estimation in zLayer = %d", i+1), 800, 600);
        canvas_dep_energy[i]->cd();

        hd_xy_dep[i]->Scale( (number_POT_year/1e9) * (1e-13/bin_mass) ); // To obtain the dose in Gy for 3e22 POT/year
        hd_xy_dep[i]->SetTitle(Form("Estimated dose in zLayer = %d", i+1));
        hd_xy_dep[i]->GetXaxis()->SetTitle("x (mm)");
        hd_xy_dep[i]->GetYaxis()->SetTitle("y (mm)");
        hd_xy_dep[i]->GetZaxis()->SetTitle("Dose/year (Gy)");
        hd_xy_dep[i]->Draw("COLZ");

        canvas_dep_energy[i]->Print("Plots_ESS_HDAnalysis.pdf");
    }

    // This is just to close the .pdf file with the plots. I have no idea how to do this otherwise
    // TCanvas* c4 = new TCanvas("c4", "Hadron dump x, y distributions", 1200, 800);
    // c4->Print("Plots_ESS_HDAnalysis.pdf]");

    c1 -> Print("Plots_ESS_HDAnalysis.pdf]");

    myApp->Run();
    return 0;
}