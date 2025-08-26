// ML_Project_Part1.C //
// Goal: Illustration of the general idea behind the Maximum Likelihood (ML) estimation method

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

void ML_Project_Part1() {
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(1111);  // Show statistics box with all info
    gStyle->SetOptFit(1111);   // Show fit parameters in stat box
    gStyle->SetOptTitle(1);    // Show title
    gStyle->SetTitleOffset(1.2, "Y");

    // Step 1: Generate a sample of 50 measurements according to a Gaussian p.d.f.
    std::cout << "--- Step 1: Generating 50 Gaussian measurements ---" << std::endl;
    const double true_mean = 0.2;   // True mean (μ)
    const double true_sigma = 0.1;  // True standard deviation (σ)
    const int N_measurements = 50;  // Number of measurements
    
    TRandom3 rng(12345);  // Seed RNG for reproducibility
    std::vector<double> measurements;
    
    for (int i = 0; i < N_measurements; ++i) {
        measurements.push_back(rng.Gaus(true_mean, true_sigma));
    }
    
    // Create a histogram to visualize the data
    double data_min = *std::min_element(measurements.begin(), measurements.end());
    double data_max = *std::max_element(measurements.begin(), measurements.end());
    double margin = 0.2 * (data_max - data_min);
    double hist_min = data_min - margin;
    double hist_max = data_max + margin;
    
    TH1F* h_data = new TH1F("h_data", "Sample of 50 Gaussian Measurements;x;Entries", 20, hist_min, hist_max);
    for (double x : measurements) {
        h_data->Fill(x);
    }
    
    // Enhance visual appearance of data histogram
    h_data->SetFillColor(38);  // Light blue
    h_data->SetLineColor(4);   // Blue
    h_data->SetLineWidth(2);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerColor(4); // Blue
    h_data->SetMarkerSize(0.8);
    
    TCanvas* c_data = new TCanvas("c_data", "Gaussian Measurements", 800, 600);
    c_data->SetFillColor(10);
    c_data->SetFrameFillColor(10);
    c_data->SetGridx();
    c_data->SetGridy();
    h_data->Draw("E1 HIST");
    c_data->Update();
    c_data->SaveAs("Gaussian_Measurements.png");
    
    // Step 2: Find the ML estimations of the mean and variance of the parent distribution
    std::cout << "\n--- Step 2: Finding ML estimates for Gaussian parameters ---" << std::endl;
    
    // For Gaussian distribution, the ML estimate of the mean is the sample mean
    double sum = std::accumulate(measurements.begin(), measurements.end(), 0.0);
    double ml_mean = sum / N_measurements;
    
    // ML estimate of variance is the biased estimator (dividing by N, not N-1)
    double sum_squared_diff = 0.0;
    for (double x : measurements) {
        sum_squared_diff += std::pow(x - ml_mean, 2);
    }
    double ml_variance = sum_squared_diff / N_measurements;
    double ml_sigma = std::sqrt(ml_variance);
    
    // Unbiased estimator of variance (for comparison)
    double unbiased_variance = sum_squared_diff / (N_measurements - 1);
    double unbiased_sigma = std::sqrt(unbiased_variance);
    
    std::cout << "  True mean (μ): " << true_mean << std::endl;
    std::cout << "  ML estimate of mean: " << ml_mean << std::endl;
    std::cout << "  True standard deviation (σ): " << true_sigma << std::endl;
    std::cout << "  ML estimate of standard deviation: " << ml_sigma << std::endl;
    std::cout << "  Unbiased estimate of standard deviation: " << unbiased_sigma << " (not used for ML)" << std::endl;
    
    // Step 3: Plot the Gaussian distribution using the ML parameters and true parameters
    std::cout << "\n--- Step 3: Plotting Gaussian PDF with ML and true parameters ---" << std::endl;
    
    // Create a canvas for the plot
    TCanvas* c_pdf = new TCanvas("c_pdf", "Gaussian PDF Comparison", 800, 600);
    
    // Function with ML parameters
    TF1* f_ml = new TF1("f_ml", "[0]*TMath::Gaus(x, [1], [2], 1)", hist_min, hist_max);
    f_ml->SetParameter(0, 1.0);  // Normalization
    f_ml->SetParameter(1, ml_mean);  // Mean
    f_ml->SetParameter(2, ml_sigma);  // Sigma
    f_ml->SetLineColor(4);  // Blue
    f_ml->SetLineWidth(2);
    
    // Function with true parameters
    TF1* f_true = new TF1("f_true", "[0]*TMath::Gaus(x, [1], [2], 1)", hist_min, hist_max);
    f_true->SetParameter(0, 1.0);  // Normalization
    f_true->SetParameter(1, true_mean);  // Mean
    f_true->SetParameter(2, true_sigma);  // Sigma
    f_true->SetLineColor(2);  // Red
    f_true->SetLineStyle(2);
    f_true->SetLineWidth(2);
    
    // Scale the functions to match the histogram
    double bin_width = (hist_max - hist_min) / 20.0;
    double scale_factor = h_data->Integral() * bin_width;
    f_ml->SetParameter(0, scale_factor);
    f_true->SetParameter(0, scale_factor);
    
    // Draw histogram and functions
    h_data->Draw();
    f_ml->Draw("SAME");
    f_true->Draw("SAME");
    
    // Add a legend
    TLegend* legend = new TLegend(0.14, 0.95, 0.35, 0.7);
    legend->AddEntry(h_data, "Data", "lep");
    legend->AddEntry(f_ml, Form("ML Fit: #mu=%.3f, #sigma=%.3f", ml_mean, ml_sigma), "l");
    legend->AddEntry(f_true, Form("True: #mu=%.3f, #sigma=%.3f", true_mean, true_sigma), "l");
    legend->Draw();
    
    c_pdf->Update();
    c_pdf->SaveAs("Gaussian_PDF_Comparison.png");
    
    // Step 4: Estimate the variance of the mean with the Monte Carlo method
    std::cout << "\n--- Step 4: Monte Carlo estimation of the variance of the mean ---" << std::endl;
    
    // Use the ML estimates as "true" parameters for MC simulation
    const double mc_true_mean = ml_mean;
    const double mc_true_sigma = ml_sigma;
    const int N_mc_trials = 1000;
    std::vector<double> mean_estimates;
    
    std::cout << "  Simulating " << N_mc_trials << " experiments with " << N_measurements << " measurements each..." << std::endl;
    
    // Generate and analyze many simulated experiments
    for (int trial = 0; trial < N_mc_trials; ++trial) {
        std::vector<double> mc_measurements;
        
        // Generate measurements for this trial
        for (int i = 0; i < N_measurements; ++i) {
            mc_measurements.push_back(rng.Gaus(mc_true_mean, mc_true_sigma));
        }
        
        // Calculate ML estimate of mean for this trial
        double mc_sum = std::accumulate(mc_measurements.begin(), mc_measurements.end(), 0.0);
        double mc_mean = mc_sum / N_measurements;
        mean_estimates.push_back(mc_mean);
    }
    
    // Create histogram of mean estimates
    double mc_min = *std::min_element(mean_estimates.begin(), mean_estimates.end());
    double mc_max = *std::max_element(mean_estimates.begin(), mean_estimates.end());
    double mc_margin = 0.2 * (mc_max - mc_min);
    
    TH1F* h_mc = new TH1F("h_mc", "Distribution of ML Mean Estimates;#hat{#mu}_{ML};Entries", 50, mc_min - mc_margin, mc_max + mc_margin);
    for (double mean_est : mean_estimates) {
        h_mc->Fill(mean_est);
    }
    
    // Enhance visual appearance of MC histogram
    h_mc->SetFillColor(46);  // Red/orange
    h_mc->SetLineColor(2);   // Red
    h_mc->SetLineWidth(2);
    h_mc->SetMarkerStyle(21);
    h_mc->SetMarkerColor(2); // Red
    h_mc->SetMarkerSize(0.7);
    
    // Calculate statistics of the distribution
    double mean_of_means = h_mc->GetMean();
    double rms_of_means = h_mc->GetRMS();
    
    // Calculate unbiased estimator of variance
    double sum_squared_diff_mc = 0.0;
    for (double mean_est : mean_estimates) {
        sum_squared_diff_mc += std::pow(mean_est - mean_of_means, 2);
    }
    double s_squared = sum_squared_diff_mc / (N_mc_trials - 1);
    double s = std::sqrt(s_squared);
    
    // Plot the histogram with a Gaussian fit
    TCanvas* c_mc = new TCanvas("c_mc", "Monte Carlo Distribution of Mean Estimates", 800, 600);
    c_mc->SetFillColor(10);
    c_mc->SetFrameFillColor(10);
    c_mc->SetGridx();
    c_mc->SetGridy();
    
    // First draw histogram
    h_mc->Draw("HIST");
    
    // Then perform fit and draw it prominently
    h_mc->Fit("gaus", "L");  // "L" for likelihood fit method
    
    // Get fit function and make it more visible
    TF1* gaussian_fit = h_mc->GetFunction("gaus");
    gaussian_fit->SetLineColor(1);  // Black
    gaussian_fit->SetLineWidth(3);  // Thicker line
    gaussian_fit->Draw("SAME");     // Draw it again on top
    
    // Add a line showing the mean and annotate it
    TLine* mean_line = new TLine(mean_of_means, 0, mean_of_means, h_mc->GetMaximum()*0.8);
    mean_line->SetLineColor(4);  // Blue
    mean_line->SetLineWidth(2);
    mean_line->SetLineStyle(2);  // Dashed
    mean_line->Draw();
    
    // Get fit parameters (already have the function above)
    double fit_mean = gaussian_fit->GetParameter(1);
    double fit_sigma = gaussian_fit->GetParameter(2);
    
    // Add text box with stats
    TPaveText* stats = new TPaveText(0.14, 0.95, 0.35, 0.7, "NDC");
    stats->SetFillColor(10);
    stats->SetBorderSize(1);
    stats->SetTextAlign(12);
    stats->AddText(Form("Mean = %.5f", mean_of_means));
    stats->AddText(Form("RMS = %.5f", rms_of_means));
    stats->AddText(Form("Fit #sigma = %.5f", fit_sigma));
    stats->AddText(Form("MC Error = %.5f", s));
    stats->Draw();
    
    c_mc->Update();
    c_mc->SaveAs("ML_Mean_Distribution.png");
    
    std::cout << "  Mean of ML mean estimates: " << mean_of_means << std::endl;
    std::cout << "  RMS of ML mean estimates: " << rms_of_means << std::endl;
    std::cout << "  Unbiased standard deviation (s): " << s << std::endl;
    std::cout << "  Gaussian fit sigma: " << fit_sigma << std::endl;
    
    // Step 5: Compare Monte Carlo error with analytical solution
    std::cout << "\n--- Step 5: Comparing Monte Carlo error with analytical solution ---" << std::endl;
    
    // Calculate analytical error: σ/√n
    double analytical_error = mc_true_sigma / std::sqrt(N_measurements);
    
    std::cout << "  Monte Carlo error (s): " << s << std::endl;
    std::cout << "  Analytical error (σ/√n): " << analytical_error << std::endl;
    std::cout << "  Ratio (MC/Analytical): " << s / analytical_error << std::endl;
    
    if (std::abs(s / analytical_error - 1.0) < 0.1) {
        std::cout << "  The Monte Carlo result is in good agreement with the analytical solution." << std::endl;
    } else if (s / analytical_error > 1.0) {
        std::cout << "  The Monte Carlo error is larger than the analytical solution." << std::endl;
    } else {
        std::cout << "  The Monte Carlo error is smaller than the analytical solution." << std::endl;
    }
    
    std::cout << "\n--- Project Part 1 Finished ---" << std::endl;
}