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
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <TLegend.h>
#include "../windowsfft/spectrum_for_RRI.h"
#include "../csvfiles/Read_CSV.h"

using std::vector;

TCanvas* draw_apnea_start_and_end(std::string file, bool savefile = kFALSE){
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
    mycanva->Divide(0, 2, 0, 0);
    mycanva->cd(1);
    double x[6] = {-10, 0, 0, 10, 10, 30};
    double y[6] = {0, 0, 1, 1, 0, 0};
    double x2[4] = {10, 10, 20, 20};
    double y2[4] = {0, 1, 1, 0};

    int N = apnea_start.size();
    TGraphErrors *gr_ap_start = new TGraphErrors(N, t1, ap_start,0,ap_start_er);
    TGraphErrors *gr_ap_end   = new TGraphErrors(N, t2, ap_end, 0, ap_end_er);
    TGraphErrors *gr_br_start = new TGraphErrors(N, t1, br_start,0 , br_start_er);
    TGraphErrors *gr_br_end   = new TGraphErrors(N, t2, br_end, 0, br_end_er);
    TMultiGraph *gr_ap = new TMultiGraph();
    TMultiGraph *gr_br = new TMultiGraph();


    gr_ap_start->SetMarkerStyle(26);
    gr_br_start->SetMarkerStyle(22);
    gr_ap_end->SetMarkerStyle(32);
    gr_br_end->SetMarkerStyle(23);

    gr_ap_start->SetMarkerColor(kRed);
    gr_ap_end->SetMarkerColor(kBlue);
    gr_ap_start->SetLineColor(kRed);
    gr_ap_end->SetLineColor(kBlue);

    gr_ap->SetTitle("apnea; time/s");

    gr_br_start->SetMarkerColor(kRed);
    gr_br_end->SetMarkerColor(kBlue);
    gr_br_start->SetLineColor(kRed);
    gr_br_end->SetLineColor(kBlue);

    gr_ap->Add(gr_ap_start);
    gr_ap->Add(gr_ap_end);
    gr_br->Add(gr_br_end);
    gr_br->Add(gr_br_start);
    gr_br->SetTitle("breathing; time/s");

    TGraph *shadow_start = new TGraph(6, x, y);
    TGraph *shadow_end = new TGraph(4, x2, y2);

    gr_ap_start->SetTitle(file.c_str());

    shadow_start->SetFillColorAlpha(kRed, 0.3);
    shadow_start->GetXaxis()->SetTitle("time/s");
    shadow_start->GetXaxis()->CenterTitle();
    shadow_end->SetFillColorAlpha(kBlue, 0.3);
    gr_ap->Draw("APL");
    shadow_start->Draw("F");
    shadow_end->Draw("F");

    mycanva->cd(2);
    gr_br->Draw("ALP");
    shadow_start->Draw("F");
    shadow_end->Draw("F");
    if(savefile)
        mycanva->SaveAs( (file + ".pdf").c_str() );

    return mycanva;
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
    TLegend *leg = new TLegend(0.7, 0.7, 1, 1);

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
void print_info(TGraph *apnea_band, int starttime){
    int N = apnea_band->GetN();
    double *y = apnea_band->GetY();
    double *x = apnea_band->GetX();
    int i;
    int time_start = 0;
    for (i = 0; i < N; ++i) {
        if(y[i] == 1 && y[i - 1] == 0){
            time_start = x[i];
            second_to_time(time_start + starttime).print();
        }
        if(y[i] == 1 && y[i + 1] == 0){
            printf(" %f \n",x[i] - time_start);
        }
    }
}
