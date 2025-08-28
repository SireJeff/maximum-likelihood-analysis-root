// ML_Project_Part2.C // Ensure your filename matches this if you use this function name

#include <iostream>
#include <vector>
#include <numeric> // For std::accumulate
#include <cmath>   // For std::log, std::exp, std::sqrt

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
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TPaveStats.h>

// Main function for the project - RENAMED TO MATCH THE FILENAME
void ML_Project_Part2(){ // <--- RENAMED HERE
    // Improve default plot aesthetics
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);    // Turn off statistics box on histograms
    gStyle->SetOptFit(1111);  // Show fit parameters in stat box if it were on
    gStyle->SetOptTitle(1);   // Show title

    // --- Data Provided ---
    std::vector<double> times = {
        4.99, 4.87, 2.59, 3.04, 3.39, 6.20, 10.61, 7.64, 3.92, 5.33, 4.85, 2.39,
        4.16, 6.74, 3.53, 5.86, 5.41, 26.25, 4.40, 10.79, 7.08, 2.86, 33.92, 3.03, 0.98,
        5.63, 4.89, 2.26, 10.49, 6.51, 7.36, 2.13, 6.45, 2.29, 21.15, 4.07, 4.34, 5.38,
        7.69, 4.93
    };
    int N_data = times.size();
    std::cout << "Number of data points (N): " << N_data << std::endl;

    // --- First: with the unbinned data ---
    std::cout << "\n--- Part 1: Unbinned Data Analysis ---" << std::endl;

    // Sub-bullet 1.1: Using the analytical solutions of the ML estimates
    std::cout << "\n--- 1.1: Analytical Unbinned ML Estimate ---" << std::endl;
    double sum_t = 0.0;
    for (double t_val : times) {
        sum_t += t_val;
    }
    std::cout << "Sum of decay times (Sigma t_i): " << sum_t << std::endl;

    double lambda_ml_unbinned = 0.0;
    if (sum_t > 1e-9) { // Avoid division by zero
        lambda_ml_unbinned = static_cast<double>(N_data) / sum_t;
    }

    double err_lambda_ml_unbinned = 0.0;
    if (N_data > 0 && lambda_ml_unbinned > 1e-9) {
        err_lambda_ml_unbinned = lambda_ml_unbinned / TMath::Sqrt(static_cast<double>(N_data));
    }

    std::cout << "Unbinned ML estimate for lambda (lambda_hat_ML): " << lambda_ml_unbinned << std::endl;
    std::cout << "Analytical error on lambda_hat_ML: " << err_lambda_ml_unbinned << std::endl;

    // Sub-bullet 1.2: Plot the Logarithm of Likelihood function (logL)
    std::cout << "\n--- 1.2: Log-Likelihood Plot and Graphical Error ---" << std::endl;

    auto logLikelihoodFunc = [&](double* x_lambda_ptr, double* /*params*/) {
        double current_lambda = x_lambda_ptr[0];
        if (current_lambda <= 1e-9) return -1e10;
        return static_cast<double>(N_data) * TMath::Log(current_lambda) - current_lambda * sum_t;
    };

    double lambda_plot_min = lambda_ml_unbinned > 3 * err_lambda_ml_unbinned ? lambda_ml_unbinned - 3 * err_lambda_ml_unbinned : 0.01;
    if (lambda_plot_min <= 0) lambda_plot_min = 0.001;
    double lambda_plot_max = lambda_ml_unbinned + 3 * err_lambda_ml_unbinned;
    if (lambda_plot_max <= lambda_plot_min) lambda_plot_max = lambda_plot_min + 0.1;

    TF1* fLogL = new TF1("fLogL", logLikelihoodFunc, lambda_plot_min, lambda_plot_max, 0);
    fLogL->SetNpx(500);

    double logL_max = fLogL->Eval(lambda_ml_unbinned);
    double target_logL_for_sigma = logL_max - 0.5;

    std::cout << "LogL_max at lambda_hat_ML (" << lambda_ml_unbinned << ") = " << logL_max << std::endl;
    std::cout << "Target LogL for 1-sigma band (LogL_max - 0.5) = " << target_logL_for_sigma << std::endl;

    double lambda_sigma_minus = lambda_ml_unbinned;
    double lambda_sigma_plus = lambda_ml_unbinned;

    if (fLogL->Eval(lambda_plot_min) < target_logL_for_sigma && lambda_ml_unbinned > lambda_plot_min) {
         lambda_sigma_minus = fLogL->GetX(target_logL_for_sigma, lambda_plot_min, lambda_ml_unbinned - 1e-6);
    } else {
        std::cout << "Warning: Could not find lower sigma bound reliably with GetX. Check plot range or function behavior." << std::endl;
    }
    if (fLogL->Eval(lambda_plot_max) < target_logL_for_sigma && lambda_ml_unbinned < lambda_plot_max) {
        lambda_sigma_plus = fLogL->GetX(target_logL_for_sigma, lambda_ml_unbinned + 1e-6, lambda_plot_max);
    } else {
         std::cout << "Warning: Could not find upper sigma bound reliably with GetX. Check plot range or function behavior." << std::endl;
    }

    std::cout << "Graphically determined 1-sigma interval for lambda:" << std::endl;
    std::cout << "  lambda_sigma_minus = " << lambda_sigma_minus << std::endl;
    std::cout << "  lambda_sigma_plus = " << lambda_sigma_plus << std::endl;
    double err_minus_graphical = lambda_ml_unbinned - lambda_sigma_minus;
    double err_plus_graphical = lambda_sigma_plus - lambda_ml_unbinned;
    std::cout << "  Asymmetric errors: -" << err_minus_graphical << " / +" << err_plus_graphical << std::endl;
    double avg_err_graphical = (err_minus_graphical + err_plus_graphical) / 2.0;
    std::cout << "  Approx. symmetric error (average of +/-): " << avg_err_graphical << std::endl;

    std::cout << "Compatibility with analytical error (" << err_lambda_ml_unbinned << "): ";
    if (TMath::Abs(avg_err_graphical - err_lambda_ml_unbinned) < 0.1 * err_lambda_ml_unbinned) {
        std::cout << "Compatible." << std::endl;
    } else {
        std::cout << "Check compatibility (graphical error " << avg_err_graphical
                  << " vs analytical " << err_lambda_ml_unbinned << ")." << std::endl;
    }

    TCanvas* c_logL = new TCanvas("c_logL", "Log-Likelihood Function for Unbinned Data", 800, 600);
    fLogL->SetTitle(Form("Log-Likelihood for #lambda; #lambda; ln L(#lambda) (N=%d, #Sigma t=%.2f)", N_data, sum_t));
    fLogL->Draw();
    TLine* line_max_val = new TLine(lambda_ml_unbinned, fLogL->GetMinimum(), lambda_ml_unbinned, logL_max);
    line_max_val->SetLineStyle(2); line_max_val->SetLineColor(kGray+2); line_max_val->Draw();
    TLine* line_sigma_level = new TLine(fLogL->GetXmin(), target_logL_for_sigma, fLogL->GetXmax(), target_logL_for_sigma);
    line_sigma_level->SetLineStyle(2); line_sigma_level->SetLineColor(kRed); line_sigma_level->Draw();
    TMarker* marker_ml = new TMarker(lambda_ml_unbinned, logL_max, 20);
    marker_ml->SetMarkerColor(kBlue); marker_ml->SetMarkerSize(1.5); marker_ml->Draw();
    TMarker* marker_s1 = new TMarker(lambda_sigma_minus, target_logL_for_sigma, 20);
    marker_s1->SetMarkerColor(kRed); marker_s1->SetMarkerSize(1.5); marker_s1->Draw();
    TMarker* marker_s2 = new TMarker(lambda_sigma_plus, target_logL_for_sigma, 20);
    marker_s2->SetMarkerColor(kRed); marker_s2->SetMarkerSize(1.5); marker_s2->Draw();
    TLegend* leg_logL = new TLegend(0.15, 0.15, 0.5, 0.35);
    leg_logL->AddEntry(marker_ml, Form("#hat{#lambda}_{ML} = %.3f", lambda_ml_unbinned), "p");
    leg_logL->AddEntry(line_sigma_level, "ln L_{max} - 0.5 level", "l");
    leg_logL->AddEntry(marker_s1, Form("1#sigma band [%.3f, %.3f]", lambda_sigma_minus, lambda_sigma_plus), "p");
    leg_logL->Draw();
    c_logL->Update();
    c_logL->SaveAs("LogLikelihoodFunction.png");

    // --- Second: with the data stored in a histogram (binned fit) ---
    std::cout << "\n\n--- Part 2: Binned Data Analysis (Histogram Fits) ---" << std::endl;

    // Sub-bullet 2.1: Histogram (t_min=0, t_max=15, nBins=3), LogL fit
    std::cout << "\n--- 2.1: Binned Fit (t_max=15, nBins=3, LogL method) ---" << std::endl;
    double t_min_hist1 = 0.0;
    double t_max_hist1 = 15.0;
    int n_bins_hist1 = 3;
    TH1F* h1 = new TH1F("h1", Form("Binned Decay (N_bins=%d);Time (t);Entries / (%.1f time units)", n_bins_hist1, (t_max_hist1-t_min_hist1)/n_bins_hist1),
                        n_bins_hist1, t_min_hist1, t_max_hist1);
    for (double t_val : times) { h1->Fill(t_val); }
    h1->SetLineColor(kBlue);
    h1->SetMarkerColor(kBlue);
    h1->SetMarkerStyle(20);
    h1->SetLineWidth(2);
    TF1* fitFunc1 = new TF1("fitFunc1", "[0]*exp(-[1]*x)", t_min_hist1, t_max_hist1);
    fitFunc1->SetParName(0, "Amplitude"); fitFunc1->SetParName(1, "#lambda");
    fitFunc1->SetParameter(1, lambda_ml_unbinned);
    double approx_norm1 = (double)N_data * ( (t_max_hist1-t_min_hist1)/n_bins_hist1 ) * lambda_ml_unbinned;
    fitFunc1->SetParameter(0, approx_norm1 > 0 ? approx_norm1 : N_data);
    fitFunc1->SetLineColor(kRed);
    fitFunc1->SetLineWidth(2);
    TFitResultPtr r1 = h1->Fit(fitFunc1, "Q S"); // L for likelihood method, Q for quiet, S for save
    double lambda_binned1 = -1, err_lambda_binned1 = -1;
    if (r1->IsValid()) {
        lambda_binned1 = r1->Parameter(1); err_lambda_binned1 = r1->ParError(1);
        std::cout << "  lambda_binned (3 bins) = " << lambda_binned1 << " +/- " << err_lambda_binned1 << std::endl;
        std::cout << "  Fit Chi2/NDF = " << r1->Chi2() << "/" << r1->Ndf() << " (Note: Using Log-Likelihood method)" << std::endl;
    } else { std::cout << "  Fit 1 (3 bins) failed or was invalid." << std::endl; }
    TCanvas* c_hist1 = new TCanvas("c_hist1", "Binned Fit (3 bins, LogL)", 800, 600);
    h1->Draw("E1"); 
    if (r1->IsValid()) {
        fitFunc1->Draw("SAME");
    }
    // Add a legend
    TLegend* leg1 = new TLegend(0.78, 0.7, 0.7, 0.78);
    leg1->AddEntry(h1, "Data", "lep");
    leg1->AddEntry(fitFunc1, "Fit", "l");
    leg1->Draw();
    c_hist1->Update();
    c_hist1->SaveAs("BinnedFit_3Bins.png");

    // Sub-bullet 2.2: Evaluation of fit (3 bins)
    std::cout << "\n--- 2.2: Evaluation of Fit (3 bins) ---" << std::endl;
    if (r1->IsValid()) {
        double integral_fit1 = fitFunc1->Integral(t_min_hist1, t_max_hist1);
        double bin_width1 = h1->GetXaxis()->GetBinWidth(1);
        double Q1 = (bin_width1 > 1e-9) ? (integral_fit1 / bin_width1) : 0.0;
        double sum_bin_contents1 = h1->GetSumOfWeights();
        std::cout << "  Integral of fitted function (Area under curve) = " << integral_fit1 << std::endl;
        std::cout << "  Bin width = " << bin_width1 << std::endl;
        std::cout << "  Calculated Q (Integral / bin_width) = " << Q1 << std::endl;
        std::cout << "  Sum of histogram bin contents (N_obs_in_range) = " << sum_bin_contents1 << std::endl;
        std::cout << "  Total entries in histogram (N_filled_in_range) = " << h1->GetEntries() << std::endl;
    }

    // Sub-bullet 2.3: Increase nBins to 5
    std::cout << "\n--- 2.3: Binned Fit (t_max=15, nBins=5, LogL method) ---" << std::endl;
    int n_bins_hist2 = 5;
    TH1F* h2 = new TH1F("h2", Form("Binned Decay (N_bins=%d);Time (t);Entries / (%.1f time units)", n_bins_hist2, (t_max_hist1-t_min_hist1)/n_bins_hist2),
                        n_bins_hist2, t_min_hist1, t_max_hist1);
    for (double t_val : times) { h2->Fill(t_val); }
    h2->SetLineColor(kGreen+2);
    h2->SetMarkerColor(kGreen+2);
    h2->SetMarkerStyle(20);
    h2->SetLineWidth(2);
    TF1* fitFunc2 = new TF1("fitFunc2", "[0]*exp(-[1]*x)", t_min_hist1, t_max_hist1);
    fitFunc2->SetParName(0, "Amplitude"); fitFunc2->SetParName(1, "#lambda");
    fitFunc2->SetParameter(1, lambda_ml_unbinned);
    double approx_norm2 = (double)N_data * ( (t_max_hist1-t_min_hist1)/n_bins_hist2 ) * lambda_ml_unbinned;
    fitFunc2->SetParameter(0, approx_norm2 > 0 ? approx_norm2 : N_data);
    fitFunc2->SetLineColor(kOrange+7);
    fitFunc2->SetLineWidth(2);
    TFitResultPtr r2 = h2->Fit(fitFunc2, "Q S"); // L for likelihood method, Q for quiet, S for save
    double lambda_binned2 = -1, err_lambda_binned2 = -1;
    if (r2->IsValid()) {
        lambda_binned2 = r2->Parameter(1); err_lambda_binned2 = r2->ParError(1);
        std::cout << "  lambda_binned (5 bins) = " << lambda_binned2 << " +/- " << err_lambda_binned2 << std::endl;
    } else { std::cout << "  Fit 2 (5 bins) failed or was invalid." << std::endl; }
    std::cout << "Observation: Compare lambda_binned (5 bins) = " << lambda_binned2
              << " with lambda_binned (3 bins) = " << lambda_binned1
              << " and unbinned lambda_ML = " << lambda_ml_unbinned << "." << std::endl;
    std::cout << "  The error err_lambda_binned (5 bins) = " << err_lambda_binned2
              << " vs err_lambda_binned (3 bins) = " << err_lambda_binned1
              << " and unbinned error = " << err_lambda_ml_unbinned << "." << std::endl;
    std::cout << "  With more bins (5 vs 3), the binned estimate might get closer to the unbinned one if statistics per bin are still sufficient." << std::endl;
    TCanvas* c_hist2 = new TCanvas("c_hist2", "Binned Fit (5 bins, LogL)", 800, 600);
    h2->Draw("E1"); 
    if (r2->IsValid()) {
        fitFunc2->Draw("SAME");
    }
    // Add a legend
    TLegend* leg2 = new TLegend(0.78, 0.7, 0.7, 0.78);
    leg2->AddEntry(h2, "Data", "lep");
    leg2->AddEntry(fitFunc2, "Fit", "l");
    leg2->Draw();
    c_hist2->Update();
    c_hist2->SaveAs("BinnedFit_5Bins.png");

    // Sub-bullet 2.4: Extend x-axis (t_max=35), nBins=350
    std::cout << "\n--- 2.4: Binned Fit (t_max=35, nBins=350, LogL method) ---" << std::endl;
    double t_min_hist3 = 0.0; double t_max_hist3 = 35.0; int n_bins_hist3 = 350;
    TH1F* h3 = new TH1F("h3", "Binned Decay (350 bins, t_{max}=35);Time (t);Entries / bin",
                        n_bins_hist3, t_min_hist3, t_max_hist3);
    for (double t_val : times) { h3->Fill(t_val); }
    h3->SetLineColor(kViolet-3);
    h3->SetMarkerColor(kViolet-3);
    h3->SetMarkerStyle(20);
    h3->SetLineWidth(2);
    std::cout << "  Total entries filled in h3 (t_max=35): " << h3->GetEntries() << " (should be " << N_data << ")" << std::endl;
    TF1* fitFunc3 = new TF1("fitFunc3", "[0]*exp(-[1]*x)", t_min_hist3, t_max_hist3);
    fitFunc3->SetParName(0, "Amplitude"); fitFunc3->SetParName(1, "#lambda");
    fitFunc3->SetParameter(1, lambda_ml_unbinned);
    fitFunc3->SetParameter(0, N_data * lambda_ml_unbinned > 0 ? N_data * lambda_ml_unbinned : N_data);
    fitFunc3->SetLineColor(kCyan+2);
    fitFunc3->SetLineWidth(2);
    TFitResultPtr r3 = h3->Fit(fitFunc3, "Q S"); // L for likelihood method, Q for quiet, S for save
    double lambda_binned3 = -1, err_lambda_binned3 = -1;
    if (r3->IsValid()) {
        lambda_binned3 = r3->Parameter(1); err_lambda_binned3 = r3->ParError(1);
        std::cout << "  lambda_binned (350 bins, t_max=35) = " << lambda_binned3 << " +/- " << err_lambda_binned3 << std::endl;
    } else { std::cout << "  Fit 3 (350 bins) failed or was invalid." << std::endl; }
    std::cout << "How results change: Compare lambda_binned (350 bins) = " << lambda_binned3
              << " with unbinned lambda_ML = " << lambda_ml_unbinned << "." << std::endl;
    std::cout << "  The error err_lambda_binned (350 bins) = " << err_lambda_binned3
              << " vs unbinned error = " << err_lambda_ml_unbinned << "." << std::endl;
    std::cout << "  With very fine binning and all data included, the binned LogL fit should closely approximate the unbinned ML result." << std::endl;
    TCanvas* c_hist3 = new TCanvas("c_hist3", "Binned Fit (350 bins, t_max=35, LogL)", 800, 600);
    c_hist3->SetTopMargin(0.15);  // Add more space for title
    h3->SetTitleOffset(1.2, "y"); // Move y-axis title further from axis
    h3->SetTitleOffset(1.2, "x"); // Move x-axis title further from axis
    h3->Draw("E1"); 
    if (r3->IsValid()) {
        fitFunc3->Draw("SAME");
    }
    // Add a legend
    TLegend* leg3 = new TLegend(0.78, 0.7, 0.7, 0.78);
    leg3->AddEntry(h3, "Data", "lep");
    leg3->AddEntry(fitFunc3, "Fit", "l");
    leg3->Draw();
    c_hist3->Update();
    c_hist3->SaveAs("BinnedFit_350Bins.png");

    std::cout << "\n--- Project Part 2 Finished ---" << std::endl;
}