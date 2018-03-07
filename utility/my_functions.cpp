/*
 * =====================================================================================
 *
 *       Filename:  my_functions.cpp
 *
 *    Description:  some Utilities
 *
 *        Version:  1.0
 *        Created:  02/28/2018 02:01:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yaopeng Ma (), sdumyp@126.com
 *   Organization:  
 *
 * =====================================================================================
 */

#include <iostream>
#include <vector>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include "../windowsfft/spectrum_for_RRI.h"

using std::vector;

void draw_apnea_start_and_end(std::string file){
    std::string rrifile = "/Users/apple/workspace/Bar_Ilan/specturm_RRI/somno/" + file + ".rri";
    std::string csvfile = "/Users/apple/Desktop/data/" + file + " - events.csv";
    SpectrumForRRI myspectrum(rrifile, 4, 1024, 1020);
    myspectrum.cal_band();
    myspectrum.cal_around_apnea(csvfile,';');
    vector< vector<double> >apnea_start, apnea_end, breath_start, breath_end;
    apnea_start = myspectrum.get_apnea_start_vec();
    apnea_end = myspectrum.get_apnea_end_vec();
    breath_start = myspectrum.get_breathing_start_vec();
    breath_end = myspectrum.get_breathing_end_vec();

    double *t1 = new double[apnea_start.size()];
    double *t2 = new double[apnea_end.size()];
    double t1_start = -10;
    double t2_start = 10;
    double delta_t = 20 / (apnea_start.size() - 1);
    t1[0] = t1_start;
    t2[0] = t2_start;
    for(unsigned int i = 1; i < apnea_start.size(); i++){
        t1[i] = t1[i - 1] + delta_t;
        t2[i] = t2[i - 1] + delta_t;
    }

    double *ap_start = new double[apnea_start.size()];
    double *ap_end   = new double[apnea_end.size()];
    double *br_start = new double[breath_start.size()];
    double *br_end   = new double[breath_end.size()];


    double *ap_start_er = new double[apnea_start.size()];
    double *ap_end_er   = new double[apnea_end.size()];
    double *br_start_er = new double[breath_start.size()];
    double *br_end_er   = new double[breath_end.size()];

    for(unsigned int i  = 0; i < apnea_end.size(); ++i){
        ap_start[i]     = mean_value(apnea_start[i]);
        ap_start_er[i]  = std_err(apnea_start[i]);
        ap_end[i]       = mean_value(apnea_end[i]);
        ap_end_er[i]    = std_err(apnea_end[i]);
        br_start[i]     = mean_value(breath_start[i]);
        br_start_er[i]  = std_err(apnea_start[i]);
        br_end[i]       = mean_value(breath_end[i]);
        br_end_er[i]    = std_err(apnea_end[i]);
    }

    TCanvas *mycanva = new TCanvas("canva","around apnea",800,600);
    mycanva->cd();
    TLegend *leg = new TLegend(0.8,0.8,0.95,0.9);
    double x[6] = {-10, 0, 0, 10, 10, 30};
    double y[6] = {0, 0, 1, 1, 0, 0};
    double x2[4] = {10, 10, 20, 20};
    double y2[4] = {0, 1, 1, 0};

    int N = apnea_start.size();
    TGraphErrors *gr_ap_start = new TGraphErrors(N, t1, ap_start,0,ap_start_er);
    TGraphErrors *gr_ap_end   = new TGraphErrors(N, t2, ap_end, 0, ap_end_er);
    TGraphErrors *gr_br_start = new TGraphErrors(N, t1, br_start,0 , br_start_er);
    TGraphErrors *gr_br_end   = new TGraphErrors(N, t2, br_end, 0, br_end_er);

    leg->AddEntry(gr_ap_start,"apnea band","pl");
    leg->AddEntry(gr_br_start, "breathing band","pl");
    leg->AddEntry(gr_ap_end, "apnea band","pl");
    leg->AddEntry(gr_br_end, "breathing band","pl");

    gr_ap_start->SetMarkerStyle(26);
    gr_br_start->SetMarkerStyle(22);
    gr_ap_end->SetMarkerStyle(32);
    gr_br_end->SetMarkerStyle(23);

    TGraph *shadow_start = new TGraph(6, x, y);
    TGraph *shadow_end = new TGraph(4, x2, y2);

    shadow_start->SetTitle(file.c_str());

    shadow_start->SetFillColorAlpha(kRed, 0.3);
    shadow_start->GetXaxis()->SetTitle("time/s");
    shadow_start->GetXaxis()->CenterTitle();
    shadow_start->GetXaxis()->SetRangeUser(-11,31);
    shadow_start->GetYaxis()->SetRangeUser(0,1);
    shadow_end->SetFillColorAlpha(kBlue, 0.3);
    shadow_start->Draw("AF");
    shadow_end->Draw("F");
    gr_ap_start->Draw("PL");
    gr_ap_end->Draw("PL");
    gr_br_start->Draw("PL");
    gr_br_end->Draw("PL");
    leg->Draw();
}

void draw_apnea_one_start_and_end(std::string file){
    std::string rrifile = "/Users/apple/workspace/Bar_Ilan/specturm_RRI/somno/" + file + ".rri";
    std::string csvfile = "/Users/apple/Desktop/data/" + file + " - events.csv";
    SpectrumForRRI myspectrum(rrifile, 4, 1024, 1020);
    myspectrum.cal_band();
    myspectrum.cal_around_apnea(csvfile,';');
    vector< vector<double> >apnea_start, apnea_end, breath_start, breath_end;
    apnea_start = myspectrum.get_apnea_start_one_vec();
    apnea_end = myspectrum.get_apnea_end_one_vec();
    breath_start = myspectrum.get_breathing_start_one_vec();
    breath_end = myspectrum.get_breathing_end_one_vec();

    double *t1 = new double[apnea_start.size()];
    double *t2 = new double[apnea_end.size()];
    double t1_start = -10;
    double t2_start = 10;
    double delta_t = 20 / (apnea_start.size() - 1);
    t1[0] = t1_start;
    t2[0] = t2_start;
    for(unsigned int i = 1; i < apnea_start.size(); i++){
        t1[i] = t1[i - 1] + delta_t;
        t2[i] = t2[i - 1] + delta_t;
    }

    double *ap_start = new double[apnea_start.size()];
    double *ap_end   = new double[apnea_end.size()];
    double *br_start = new double[breath_start.size()];
    double *br_end   = new double[breath_end.size()];


    double *ap_start_er = new double[apnea_start.size()];
    double *ap_end_er   = new double[apnea_end.size()];
    double *br_start_er = new double[breath_start.size()];
    double *br_end_er   = new double[breath_end.size()];

    for(unsigned int i  = 0; i < apnea_end.size(); ++i){
        ap_start[i]     = mean_value(apnea_start[i]);
        ap_start_er[i]  = std_err(apnea_start[i]);
        ap_end[i]       = mean_value(apnea_end[i]);
        ap_end_er[i]    = std_err(apnea_end[i]);
        br_start[i]     = mean_value(breath_start[i]);
        br_start_er[i]  = std_err(apnea_start[i]);
        br_end[i]       = mean_value(breath_end[i]);
        br_end_er[i]    = std_err(apnea_end[i]);
    }

    TCanvas *mycanva = new TCanvas("canva","around apnea",800,600);
    mycanva->cd();
    double x[6] = {-10, 0, 0, 10, 10, 30};
    double y[6] = {0, 0, 1, 1, 0, 0};
    double x2[4] = {10, 10, 20, 20};
    double y2[4] = {0, 1, 1, 0};

    int N = apnea_start.size();
    TGraphErrors *gr_ap_start = new TGraphErrors(N, t1, ap_start,0,ap_start_er);
    TGraphErrors *gr_ap_end   = new TGraphErrors(N, t2, ap_end, 0, ap_end_er);
    TGraphErrors *gr_br_start = new TGraphErrors(N, t1, br_start,0 , br_start_er);
    TGraphErrors *gr_br_end   = new TGraphErrors(N, t2, br_end, 0, br_end_er);
    TLegend *leg = new TLegend(0.8,0.8,0.95,0.9);

    gr_ap_start->SetMarkerStyle(26);
    gr_br_start->SetMarkerStyle(22);
    gr_ap_end->SetMarkerStyle(32);
    gr_br_end->SetMarkerStyle(23);


    leg->AddEntry(gr_ap_start,"apnea band","pl");
    leg->AddEntry(gr_br_start, "breathing band","pl");
    leg->AddEntry(gr_ap_end, "apnea band","pl");
    leg->AddEntry(gr_br_end, "breathing band","pl");

    TGraph *shadow_start = new TGraph(6, x, y);
    TGraph *shadow_end = new TGraph(4, x2, y2);
    shadow_start->SetTitle((file + " one_event").c_str());
    shadow_start->SetFillColorAlpha(kRed, 0.3);
    shadow_start->GetXaxis()->SetTitle("time/s");
    shadow_start->GetXaxis()->CenterTitle();
    shadow_start->GetXaxis()->SetRangeUser(-11,31);
    shadow_start->GetYaxis()->SetRangeUser(0,1);
    shadow_end->SetFillColorAlpha(kBlue, 0.3);
    shadow_start->Draw("AF");
    shadow_end->Draw("F");
    gr_ap_start->Draw("PL");
    gr_ap_end->Draw("PL");
    gr_br_start->Draw("PL");
    gr_br_end->Draw("PL");
    leg->Draw();
}
