/*
 * =====================================================================================
 *
 *       Filename:  spectrum_for_RRI.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  02/27/2018 15:33:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yaopeng Ma (), sdumyp@126.com
 *   Organization:
 *
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TVirtualFFT.h>

#include "spectrum_for_RRI.h"
#include "../csvfiles/Read_CSV.h"


using std::string;
using std::cout;
using std::endl;
using std::ifstream;

SpectrumForRRI::SpectrumForRRI(string filepath,int freq,int window_size, int overflap)
    :frequency(freq), windowsize(window_size),over(overflap),filename(filepath)
{
    spectrum        = NULL;
    power_spectrum  = NULL;
    log_of_spectrum = NULL;
    frequency_fft   = NULL;
    R_R_position    = NULL;
    R_R_interval    = NULL;
    resample_t      = NULL;
    newRRI          = NULL;
    breathing       = NULL;
    apnea           = NULL;
    show_apnea      = NULL;
    ave_apnea       = NULL;
    read_data();
    cal_resample_time();
    cal_resample();
    windows_num = total_num_after_resample / (windowsize - overflap);
    br_fr_start = 0.1;
    br_fr_end   = 0.25;
    ap_fr_start = 0.0039;
    ap_fr_end   = 0.1;
    apnea_event = new double[total_num_after_resample];
    for(int i = 0; i < total_num_after_resample; i++)
        apnea_event[i] = 0;
}

void SpectrumForRRI::read_data(){
    std::ifstream in(filename.c_str());
    string line;
    in >> line;
    // ___________________________
    in >> line;
    int position_of_equal;
    // find the position of "="
    position_of_equal = line.find("=");
    line.erase(0,position_of_equal+1);
    num_of_RR_peaks = atoi(line.c_str());
    in >> line;
    in >> line;
    position_of_equal = line.find("=");
    mytime start = str2time(line.erase(0,position_of_equal+1));
    cout << "time start at "<< line << endl;
    starttime = start.get_total_second();
    for (int i = 4;  i < 7; ++i) {
        in >> line;
    }
    R_R_interval = new double[num_of_RR_peaks + 1];
    R_R_position = new double[num_of_RR_peaks + 1];

    double delta_time;
    in >> delta_time;
    in >> line;
    double factor = 1.0/1000.0;
    R_R_interval[0] = R_R_position[0] = 0;
    R_R_interval[1] = R_R_position[1] = delta_time*factor;
    for(int i = 2; i < num_of_RR_peaks + 1; i++){
        in >> delta_time;
        in >> line;
        R_R_interval[i] = delta_time*factor;
        R_R_position[i] = R_R_interval[i] + R_R_position[i-1];
    }
    in.close();
}

void SpectrumForRRI::cal_resample_time(){
    double endtime = R_R_position[num_of_RR_peaks];
    int nn = endtime * frequency;
    time_distance = 1./(double)frequency;
    std::cout << "time distance is "<<time_distance << std::endl;
    total_num_after_resample = nn;
    int i ;
    resample_t = new double[nn];
    for(i = 0; i < nn; i++){
        resample_t[i] = (double)i * time_distance;
    }
}

void SpectrumForRRI::cal_resample(){
    newRRI = new double[total_num_after_resample];
    int inow = 0;
    int i ;
    for (i = 0; i < total_num_after_resample; i++) {
        newRRI[i] = cal_value_at_t(resample_t[i],inow);
        if(resample_t[i] > R_R_position[inow+1]){
            inow++;
        }
    }
}

double SpectrumForRRI::cal_value_at_t(double t, int istart){
    if (t > R_R_position[istart + 1]){
        while(t > R_R_position[istart + 1]){
            istart ++ ;
        }
    }
    else if (t < R_R_position[istart]){
        while( t < R_R_position[istart] ){
            istart --;
        }
    }
    if(t == R_R_position[istart]){
        return R_R_interval[istart];
    }
    if(t == R_R_position[istart+1]){
        return R_R_interval[istart + 1];
    }
    /*
       line function : y = k(t - t0) + y0;
       k = (y1 - y0) / (x1 - x0)
       */
    if(istart == 0){
        return t;
    }
    double slope = (R_R_interval[istart + 1] - R_R_interval[istart]);
    slope = slope / (R_R_interval[istart + 1]);
    double cal_value = slope * (t - R_R_position[istart]) + R_R_interval[istart];
    return cal_value;
}

void SpectrumForRRI::apply_FFT(int i, bool dofft, bool save){
    fft_now = i;
    int istart, iend;
    int n;
    istart = (windowsize - over) * i;
    iend   = istart +windowsize -1;
    if( i == windows_num - 1 )
        iend = total_num_after_resample - 1;
    n = iend - istart + 1;
    freq_num = n/2;
    num_of_frequency = freq_num;
    if(save){
        std::ofstream file;
        file.open("rri_now.dat");

        int j;
        for (j = 0; j < n; ++j) {
            file << newRRI[istart + j] << std::endl;
        }
        file.close();
    }
    if(frequency_fft == NULL){
        frequency_fft = new double[n/2];
        for(int i = 0; i < n/2 ; i++){
            frequency_fft[i] = ((double)frequency *(double)i)/((double)n);
        }
    }
    if(dofft){
        TVirtualFFT *myfft = TVirtualFFT::FFT(1,&n,"R2C");
        if(log_of_spectrum == NULL)
            log_of_spectrum = new double[freq_num];
        if(power_spectrum == NULL)
            power_spectrum  = new double[freq_num];
        myfft->SetPoints(&newRRI[istart]);
        myfft->Transform();
        double im,re;
        for(int i = 0; i < n/2; i++){
            myfft->GetPointComplex(i,re,im);
            power_spectrum[i] = (re*re+im*im)/n;
            //            std::cout << "n = "<<n <<" i =" <<i <<" freq = "<<frequency_fft[i] <<"power_spectrum = " << power_spectrum[i]<< std::endl;
            log_of_spectrum[i] = TMath::Log(power_spectrum[i]);
        }
        delete myfft;
    }
}

void SpectrumForRRI::draw_fft_now(){
    TGraph *gr = new TGraph(freq_num-1, &frequency_fft[1],&power_spectrum[1]);
    gr->GetXaxis()->SetTitle("frequency/HZ");
    gr->SetTitle("fft result at windows");
    gr->Draw("AL");
}

void SpectrumForRRI::cal_spectrum(){
    int i, j;
    apply_FFT(0);
    spectrum = new TH2D("spectrogram","spectrogram",windows_num,0,windows_num*(windowsize-over)*time_distance,num_of_frequency,0,frequency/2);
    double time_now = (double)(windowsize-over)*time_distance / 2.;
    for(i = 0; i < windows_num; i++){
        if(i != 0){
            apply_FFT(i);
        }
        for(j = 0; j < num_of_frequency-1; j++){
            if(j == 0){
                spectrum->Fill(time_now,(frequency_fft[j+1]+frequency_fft[j])/2., 0);
            }
            else{
                spectrum->Fill(time_now,(frequency_fft[j+1]+frequency_fft[j])/2., power_spectrum[j]);
            }
        }
        spectrum->Fill(time_now,
                (frequency_fft[num_of_frequency-1]+2.)/2.,
                log_of_spectrum[num_of_frequency-1]);
        time_now = time_now + (windowsize-over)*time_distance;
    }
}

void SpectrumForRRI::cal_band(){
    double *ti = new double[windows_num - 1];
    double *br   = new double[windows_num - 1];
    double *ap   = new double[windows_num - 1];
    apply_FFT(0,kFALSE);
    int b_start, b_end, a_start, a_end;
    b_start = br_fr_start / frequency_fft[1] + 1;
    b_end   = br_fr_end   / frequency_fft[1] + 1;
    a_start = ap_fr_start / frequency_fft[1] + 1;
    a_end   = ap_fr_end   / frequency_fft[1] + 1;
    std::cout << frequency_fft[1] << std::endl;
    std::cout << "breathing : " <<frequency_fft[b_start] <<"HZ to" <<frequency_fft[b_end] << " HZ" << std::endl;
    std::cout << "apnea : "<<frequency_fft[a_start] <<"HZ to "<<frequency_fft[a_end] << "HZ" << std::endl;
    double frac = 1 / (double)frequency;
    double sum_br = 0;
    double sum_ap = 0;
    double sum_all = 0;
    int j;
    for(int i = 0; i < windows_num-1; i++){
        ti[i] = (double)(i * (windowsize - over) + windowsize/2);
        ti[i] *= frac;
        apply_FFT(i,kTRUE);
        /* for(int j = 0; j < num_of_frequency; j++){ */
        /*     if(j >= b_start && j <= b_end){ */
        /*         sum_br += power_spectrum[j]; */
        /*     } */
        /*     else if(j >= a_start && j <= a_end){ */
        /*         sum_ap += power_spectrum[j]; */
        /*     } */
        /* } */
        for (j = b_start;  j<= b_end; j++) {
            sum_br += power_spectrum[j];
        }
        for (j = a_start; j <= a_end; j++) {
            sum_ap += power_spectrum[j];
        }
        for (j = 1; j < num_of_frequency; ++j) {
            sum_all += power_spectrum[j];
        }
        /* std::cout << sum_br << " " << sum_ap << " " << sum_all << std::endl; */
        br[i]   = sum_br / sum_all;
        ap[i]   = sum_ap / sum_all;
        /* std::cout << br[i] + ap[i] << std::endl; */
        sum_ap  = 0;
        sum_br  = 0;
        sum_all = 0;
    }
    breathing = new TGraph(windows_num-1, ti, br);
    apnea     = new TGraph(windows_num-1, ti, ap);
    apnea->SetTitle("apnea band");
    apnea->GetXaxis()->SetTitle("Time/S");
    apnea->GetXaxis()->CenterTitle();
    apnea->SetLineColor(kBlue);
    apnea->SetLineWidth(1);
    apnea->GetYaxis()->SetTitle("ratio");
    apnea->GetYaxis()->CenterTitle();
    apnea->GetYaxis()->SetTitleSize(0.06);
    breathing->SetTitle("breathing band");
    breathing->GetXaxis()->SetTitle("Time/S");
    breathing->GetXaxis()->CenterTitle();
    breathing->GetYaxis()->SetTitle("ratio");
    breathing->GetYaxis()->CenterTitle();
    breathing->GetYaxis()->SetTitleSize(0.06);
    breathing->SetLineColor(kRed);
    breathing->SetLineWidth(1);
    delete[] ti;
    delete[] ap;
    delete[] br;
}

void SpectrumForRRI::cal_apnea(std::string filename, char c){
    Event Ap_event(filename, c);
    int i = 0;
    int num_of_apnea;
    for(i = 0; i < Ap_event.get_number_of_apnea(); i++){
        int apnea_start = Ap_event.get_apnea_starttime(i);
        apnea_start     = apnea_start - starttime ;
        double duration = Ap_event.get_apnea_duration(i);
        int idur        = (int)((double)duration*(double)frequency);
        int istart      = apnea_start*frequency;
        for(int j = 0; j < idur; j++)
            apnea_event[istart + j] = 1;
    }
    for(i = 0; i < Ap_event.get_number_of_hypopnea(); i++){
        int hyp_start   = Ap_event.get_hypopnea_starttime(i);
        hyp_start       = hyp_start - starttime ;
        double duration = Ap_event.get_hypopnea_duration(i);
        int idur        = (int)((double)duration*(double)frequency);
        int istart      = hyp_start*frequency;
        for (int j = 0; j < idur; j++) {
            apnea_event[istart + j] = 0.5;
        }
    }
    show_apnea = new TGraph(total_num_after_resample, resample_t, apnea_event);
    show_apnea->SetTitle("apnea events");
    show_apnea->GetXaxis()->SetTitle("time/s");
    show_apnea->GetXaxis()->CenterTitle();
    double *time = new double[windows_num - 1];
    double *average = new double[windows_num - 1];
    double frac = 1 / (double)frequency;
    for(int i = 0; i < windows_num-1; i++){
        num_of_apnea = 0;
        int windows_start = i * (windowsize - over);
        average[i] = 0;
        time[i] = (double)(i * (windowsize - over) + windowsize/2);
        time[i] *= frac;
        for(int j = 0; j < windowsize; j++){
            average[i] += apnea_event[windows_start + j];
            if (apnea_event[windows_start + j] == 1 && apnea_event[windows_start+j+1] == 0) {
                num_of_apnea ++;
            }
        }
        //        std::cout << "windows " << i+1 << " has " << num_of_apnea <<" apneas" << std::endl;
        average[i] = average[i] / windowsize;
    }
    ave_apnea = new TGraph(windows_num-1, time, average);
    ave_apnea->SetTitle("ave apnea");
    ave_apnea->GetXaxis()->SetTitle("time/s");
    ave_apnea->GetXaxis()->CenterTitle();
}

void SpectrumForRRI::cal_around_apnea(std::string filename, char c){
    Event ap_event(filename, c);
    double time_one;
    double time_delta;
    int num_of_apnea = 0;
    vector<double> bre_start;
    vector<double> bre_end;
    vector<double> apn_start;
    vector<double> apn_end;
    if(apnea != NULL){
        int num = apnea->GetN();
        double *y1 = apnea->GetY();
        double *y2 = breathing->GetY();
        time_one = apnea->GetX()[0];
        time_delta = apnea->GetX()[1] - apnea->GetX()[0];
        for(int i = 0 ; i < ap_event.get_number_of_apnea(); i++){
            double apnea_start = ap_event.get_apnea_starttime(i);
            apnea_start -= starttime;
            double apnea_end = apnea_start + ap_event.get_apnea_duration(i);
            if(apnea_end + 10 <= apnea->GetX()[num - 1]){
                num_of_apnea ++ ;
                int istart = ((apnea_start - 10) - time_one) / time_delta;
                int iend   = ((apnea_start + 10) - time_one) / time_delta;
                int jstart = ((apnea_end - 10) - time_one) / time_delta;
                int jend   = ((apnea_end + 10) - time_one) / time_delta;
                if(i == 0){
                    apneastart.resize(iend - istart + 1);
                    apneaend.resize(jend - jstart + 1);
                    breathstart.resize(iend - istart + 1);
                    breathend.resize(jend - jstart + 1);
                    std::cout << "resize finish" << std::endl;
                    for(int j = istart; j <= iend; j++){
                        apneastart[j - istart].resize(ap_event.get_number_of_apnea());
                        breathstart[j - istart].resize(ap_event.get_number_of_apnea());
                    }
                    for(int j = jstart; j <= jend; j++){
                        apneaend[j - jstart].resize(ap_event.get_number_of_apnea());
                        breathend[j - jstart].resize(ap_event.get_number_of_apnea());
                    }
                    std::cout << "resize finish" << std::endl;
                }
                for(int j = istart; j <= iend; j++){
                    apneastart[j - istart][i] = y1[j];
                    breathstart[j - istart][i] = y2[j];
                }
                for(int j = jstart; j <= jend; j++){
                    apneaend[j - jstart][i] = y1[j];
                    breathend[j - jstart][i] = y2[j];
                }
            }
        }
        for(unsigned int i = 0; i < apneastart.size(); ++i){
            apneastart[i].resize(num_of_apnea);
            apneaend[i].resize(num_of_apnea);
            breathstart[i].resize(num_of_apnea);
            breathend[i].resize(num_of_apnea);
        }
    }
}

SpectrumForRRI::~SpectrumForRRI(){
    delete[] R_R_position;
    delete[] R_R_interval;
    delete[] resample_t;
    delete[] newRRI;
    delete[] apnea_event;
    if(log_of_spectrum != NULL) delete[] log_of_spectrum;
    if(frequency_fft   != NULL) delete[] frequency_fft;
    if(breathing       != NULL) delete breathing;
    if(apnea           != NULL) delete apnea;
    if(power_spectrum  != NULL) delete[] power_spectrum;
    if(show_apnea      != NULL) delete show_apnea;
    if(ave_apnea       != NULL) delete ave_apnea;
}

double mean_value(vector<double> values){
    unsigned int i;
    double sum = 0;
    for (i = 0; i < values.size(); ++i) {
        sum += values[i];
    }
    return sum/values.size();
}
double mean_value(double *values, int size){
    int i;
    double sum = 0;
    for (i = 0; i < size; ++i) {
        sum += values[i];
    }
    return sum / size;
}
double std_dev(vector<double> values){
    unsigned int i;
    double sum = 0;
    for (i = 0; i < values.size(); ++i) {
        sum += values[i] * values[i];
    }
    double ave = mean_value(values);
    return TMath::Sqrt(sum / values.size() - ave * ave);
}
double std_dev(double *values, int size){
    int i;
    double sum = 0;
    for(i = 0; i < size; ++i){
        sum += values[i] * values[i];
    }
    double ave = mean_value(values, size);
    return TMath::Sqrt(sum / size +  ave*ave);
}

