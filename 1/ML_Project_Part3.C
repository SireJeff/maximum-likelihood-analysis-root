// ML_Project_Part3.C //
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <random>

#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

void ML_Project_Part3() {
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleOffset(1.2, "Y");

    // --- Provided Data ---
    std::vector<double> times = {
        4.99, 4.87, 2.59, 3.04, 3.39, 6.20, 10.61, 7.64, 3.92, 5.33, 4.85, 2.39,
        4.16, 6.74, 3.53, 5.86, 5.41, 26.25, 4.40, 10.79, 7.08, 2.86, 33.92, 3.03, 0.98,
        5.63, 4.89, 2.26, 10.49, 6.51, 7.36, 2.13, 6.45, 2.29, 21.15, 4.07, 4.34, 5.38,
        7.69, 4.93
    };
    const int N_data = times.size();

    std::cout << "--- Part 1: Binned Data Analysis with Poisson Likelihood ---" << std::endl;
    std::cout << "\n--- 1.1: Binned Fit (t_max=15, nBins=3, Poisson method) ---" << std::endl;

    double t_min_hist1 = 0.0;
    double t_max_hist1 = 15.0;
    int n_bins_hist1 = 3;
    double bin_width1 = (t_max_hist1 - t_min_hist1) / n_bins_hist1;

    TH1F* h1 = new TH1F("h1", Form("Binned Decay (N_bins=%d, Poisson);Time (t);Entries / (%.1f time units)", n_bins_hist1, bin_width1),
                        n_bins_hist1, t_min_hist1, t_max_hist1);
    for (double t_val : times) {
        if (t_val >= t_min_hist1 && t_val < t_max_hist1) {
            h1->Fill(t_val);
        }
    }
    std::cout << "Number of entries in histogram h1: " << h1->GetEntries() << std::endl;

    TF1* fitFunc1 = new TF1("fitFunc1", "[0]*exp(-[1]*x)", t_min_hist1, t_max_hist1);
    fitFunc1->SetParName(0, "Amplitude");
    fitFunc1->SetParName(1, "#lambda");
    fitFunc1->SetParameter(1, 0.15);
    fitFunc1->SetParameter(0, h1->GetEntries());
    TFitResultPtr r1 = h1->Fit(fitFunc1, "LME S");;

    double lambda_binned1 = -1, err_lambda_binned1 = -1;
    if (r1->IsValid()) {
        lambda_binned1 = r1->Parameter(1);
        err_lambda_binned1 = r1->ParError(1);
        std::cout << "  Poisson LML fit result for lambda (3 bins) = " << lambda_binned1 << " +/- " << err_lambda_binned1 << std::endl;
    } else {
        std::cout << "  Fit 1 (3 bins) failed or was invalid." << std::endl;
    }

    TCanvas* c1 = new TCanvas("c1", "Binned Fit (3 bins, Poisson LML)", 800, 600);
    h1->Draw("E1");
    if (r1->IsValid()) {
        fitFunc1->Draw("SAME");
    }
    c1->Update();
    c1->SaveAs("BinnedFit_3Bins_Poisson.png");

    std::cout << "\n--- 1.2: Binned Fit (t_max=15, nBins=5, Poisson method) ---" << std::endl;
    int n_bins_hist2 = 5;
    double bin_width2 = (t_max_hist1 - t_min_hist1) / n_bins_hist2;

    TH1F* h2 = new TH1F("h2", Form("Binned Decay (N_bins=%d, Poisson);Time (t);Entries / (%.1f time units)", n_bins_hist2, bin_width2),
                        n_bins_hist2, t_min_hist1, t_max_hist1);
    for (double t_val : times) {
        if (t_val >= t_min_hist1 && t_val < t_max_hist1) {
            h2->Fill(t_val);
        }
    }
    std::cout << "Number of entries in histogram h2: " << h2->GetEntries() << std::endl;

    TF1* fitFunc2 = new TF1("fitFunc2", "[0]*exp(-[1]*x)", t_min_hist1, t_max_hist1);
    fitFunc2->SetParName(0, "Amplitude");
    fitFunc2->SetParName(1, "#lambda");
    fitFunc2->SetParameter(1, lambda_binned1 > 0 ? lambda_binned1 : 0.15);
    fitFunc2->SetParameter(0, h2->GetEntries());
    TFitResultPtr r2 = h2->Fit(fitFunc2, "LME S");

    double lambda_binned2 = -1, err_lambda_binned2 = -1;
    if (r2->IsValid()) {
        lambda_binned2 = r2->Parameter(1);
        err_lambda_binned2 = r2->ParError(1);
        std::cout << "  Poisson LML fit result for lambda (5 bins) = " << lambda_binned2 << " +/- " << err_lambda_binned2 << std::endl;
    } else {
        std::cout << "  Fit 2 (5 bins) failed or was invalid." << std::endl;
    }

    TCanvas* c2 = new TCanvas("c2", "Binned Fit (5 bins, Poisson LML)", 800, 600);
    h2->Draw("E1");
    if (r2->IsValid()) {
        fitFunc2->Draw("SAME");
    }
    c2->Update();
    c2->SaveAs("BinnedFit_5Bins_Poisson.png");

    std::cout << "\n--- 1.3: Comparison with LS fit ---" << std::endl;
    TFitResultPtr r_LS = h2->Fit(fitFunc2, "LS Q S"); // LS fit on the 5-bin histogram

    double lambda_LS = -1, err_lambda_LS = -1;
    if (r_LS->IsValid()) {
        lambda_LS = r_LS->Parameter(1);
        err_lambda_LS = r_LS->ParError(1);
        std::cout << "  LS fit result for lambda (5 bins) = " << lambda_LS << " +/- " << err_lambda_LS << std::endl;
        std::cout << "  Chi2/NDF for LS fit: " << r_LS->Chi2() << "/" << r_LS->Ndf() << std::endl;
    } else {
        std::cout << "  LS fit failed or was invalid." << std::endl;
    }

    TCanvas* c_LS = new TCanvas("c_LS", "Binned Fit (5 bins, LS)", 800, 600);
    h2->Draw("E1");
    if (r_LS->IsValid()) {
        fitFunc2->Draw("SAME");
    }
    c_LS->Update();
    c_LS->SaveAs("BinnedFit_5Bins_LS.png");

    std::cout << "\n--- Part 2: Monte Carlo Simulation and Bias Study ---" << std::endl;
    std::cout << "\n--- 2.1: MC Simulation of lambda ---" << std::endl;

    const double true_lambda = 0.15;
    const int N_trials = 1000;
    const int N_events_per_trial = 40;
    std::vector<double> lambda_estimates;

    TRandom3 rng(0); // Seed RNG for reproducibility

    for (int i = 0; i < N_trials; ++i) {
        double sum_t = 0.0;
        for (int j = 0; j < N_events_per_trial; ++j) {
            sum_t += rng.Exp(1.0 / true_lambda); // Generate exponential variate
        }
        if (sum_t > 0) {
            lambda_estimates.push_back(static_cast<double>(N_events_per_trial) / sum_t);
        }
    }
    
    TH1F* h_mc = new TH1F("h_mc", "Distribution of #hat{#lambda}_{ML} Estimates; #hat{#lambda}_{ML}; Entries", 50, 0, 0.5);
    for (double lambda_est : lambda_estimates) {
        h_mc->Fill(lambda_est);
    }

    TCanvas* c_mc_dist = new TCanvas("c_mc_dist", "MC Distribution of lambda_hat", 800, 600);
    h_mc->Draw();
    h_mc->Fit("gaus");
    c_mc_dist->Update();
    c_mc_dist->SaveAs("MC_Lambda_Distribution.png");

    double mean_mc = h_mc->GetMean();
    double sigma_mc = h_mc->GetRMS();
    double bias = mean_mc - true_lambda;

    std::cout << "  True lambda used in MC: " << true_lambda << std::endl;
    std::cout << "  Mean of MC lambda estimates: " << mean_mc << std::endl;
    std::cout << "  Bias (Mean - True): " << bias << std::endl;
    std::cout << "  RMS of MC distribution (Statistical Error): " << sigma_mc << std::endl;
    std::cout << "  Predicted error from Part 1 analytical formula: " << true_lambda / TMath::Sqrt(static_cast<double>(N_events_per_trial)) << std::endl;
    std::cout << "  Ratio of Bias / RMS: " << TMath::Abs(bias / sigma_mc) << std::endl;
    std::cout << "  If Bias/RMS is small, the estimator is approximately unbiased." << std::endl;

    std::cout << "\n--- Part 3: Relationship between LS and LML errors ---" << std::endl;
    std::cout << "\n--- 3.1: Chi-squared vs N_obs ---" << std::endl;

    std::vector<double> n_obs_list;
    std::vector<double> chi2_diff_list;
    std::vector<double> n_obs_err_list;
    std::vector<double> chi2_diff_err_list;
    
    TCanvas* c_chi2_vs_n = new TCanvas("c_chi2_vs_n", "Chi-squared difference vs N_obs", 800, 600);

    const int n_bins_chi2 = 10;
    TH1F* h_chi2 = new TH1F("h_chi2", "Chi-squared vs N_obs;N_{obs};#chi^{2}_{diff}", n_bins_chi2, 0, 20);

    for (int i = 0; i < 1000; ++i) {
        int n_obs = TRandom3(0).Poisson(5.0);
        if (n_obs > 0) {
            double mu = n_obs * 1.1;
            double sigma_sq = n_obs;
            if (sigma_sq > 0) {
                double chi2 = TMath::Power(n_obs - mu, 2) / sigma_sq;
                double logL_diff = 2 * (mu - n_obs + n_obs * TMath::Log(n_obs/mu)); // 2*Delta_ln(L)
                h_chi2->Fill(n_obs, TMath::Abs(chi2 - logL_diff));
            }
        }
    }
    h_chi2->Draw("E1");
    
    // Fit with a linear function for trend visualization
    TF1* fit_chi2 = new TF1("fit_chi2", "[0] + [1]*x", 0, 20);
    h_chi2->Fit(fit_chi2, "Q");
    c_chi2_vs_n->Update();
    c_chi2_vs_n->SaveAs("Chi2_vs_N_obs.png");
    
    std::cout << "\n--- Project Part 3 Finished ---" << std::endl;
}