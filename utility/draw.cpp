#include <TH2D.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "../windowsfft/spectrum_for_RRI.h"



void draw_pad(TH2D *hist2d,bool isLog){
    double contours[50];
    int n = 50;
    int i;
    for (i = 0; i < n; ++i) {
        contours[i] = -20 + i*0.5;
    }

    gPad->SetLogy(isLog);
    hist2d->SetContour(n,contours);
    hist2d->GetXaxis()->SetTitle("time/s");
    hist2d->GetYaxis()->SetTitle("freq/HZ");
    hist2d->GetYaxis()->CenterTitle();
    hist2d->GetXaxis()->CenterTitle();
    //    gStyle->SetPalette(53);
    hist2d->SetStats(kFALSE);
    hist2d->DrawCopy("cont4");
}
void draw_apnea(string filename,double tstart = 0, double tend = 0){
    TCanvas *c = new TCanvas("c1","c1", 800, 600);
    c->Divide(1, 3, 0, 0);
    string datafile = "/Users/apple/Desktop/data/" + filename + " - events.csv";
    string rrifile = "/Users/apple/workspace/Bar_Ilan/specturm_RRI/somno/" + filename + ".rri";
    SpectrumForRRI *myrrihand = new SpectrumForRRI(rrifile, 4, 1024, 1020);
    std::cout << "1" << std::endl;
    myrrihand->cal_apnea(datafile);
    myrrihand->cal_band();
    std::cout << "2" << std::endl;
    if(tstart != tend){
        myrrihand->get_apnea_band()->GetXaxis()->SetRangeUser(tstart,tend);
        myrrihand->get_breathing_band()->GetXaxis()->SetRangeUser(tstart,tend);
        myrrihand->get_apnea()->GetXaxis()->SetRangeUser(tstart,tend);
    }
    c->cd(2);
    myrrihand->get_apnea_band()->Draw("AL");
    c->cd(1);
    myrrihand->get_breathing_band()->Draw("AL");
    c->cd(3);
    myrrihand->get_apnea()->Draw("AL");
    gPad->Update();
    c->SaveAs((filename + "_apnea.pdf").c_str());
}

void draw_same(string filename,bool ave = false, double tstart = 0, double tend = 0){
    TCanvas *c = new TCanvas("c1","c1", 800, 600);
    TLegend *leg = new TLegend(0.7,0.8,0.9,0.9);
    leg->SetFillColorAlpha(kWhite, 0);
    leg->SetBorderSize(0);
    string datafile = "/Users/apple/Desktop/data/" + filename + " - events.csv";
    string rrifile = "/Users/apple/workspace/Bar_Ilan/specturm_RRI/somno/" + filename + ".rri";
    SpectrumForRRI *myrrihand = new SpectrumForRRI(rrifile, 4, 1024, 1020);
    myrrihand->cal_apnea(datafile);
    myrrihand->cal_band();
    TGraph *apn = NULL;
    if(ave){
        apn = myrrihand->get_ave_apnea();
    }
    else
        apn = myrrihand->get_apnea();
    if(tstart != tend){
        myrrihand->get_apnea_band()->GetXaxis()->SetRangeUser(tstart,tend);
        myrrihand->get_breathing_band()->GetXaxis()->SetRangeUser(tstart,tend);
        apn->GetXaxis()->SetRangeUser(tstart,tend);
    }
    c->cd();
    c->SetGridx();
    c->SetGridy();
    leg->AddEntry(apn,"apnea event", "f");
    leg->AddEntry(myrrihand->get_breathing_band(), "breathing band", "l");
    leg->AddEntry(myrrihand->get_apnea_band(),"apnea band", "l");
    apn->SetFillColorAlpha(kOrange, 0.3);
    apn->GetYaxis()->SetRangeUser(0, 1);
    apn->Draw("AB");
    myrrihand->get_breathing_band()->Draw("L");
    myrrihand->get_apnea_band()->Draw("L");
    leg->Draw();
}
