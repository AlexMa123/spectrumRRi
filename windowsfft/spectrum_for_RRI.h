/*
 * =====================================================================================
 *
 *       Filename:  spectrum_for_RRI.h
 *
 *    Description:  calculate windows fft of RRI file
 *
 *        Version:  1.0
 *        Created:  02/27/2018 15:33:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yaopeng Ma (), sdumyp@126.com
 *   Organization:
 *
 * =====================================================================================
 */
#include <iostream>
#include <TH2D.h>
#include <TGraph.h>
#include <TVirtualFFT.h>

#include "../csvfiles/Read_CSV.h"

#ifndef SPECTRUM_FOR_RRI_H
#define SPECTRUM_FOR_RRI_H

class SpectrumForRRI
{
    private:
        /*
         *  about RRI
         */
        const int frequency;
        const int windowsize;
        const int over;
        const std::string filename;
        int num_of_RR_peaks;
        int total_num_after_resample;
        double time_distance;
        double *R_R_position;
        double *R_R_interval;
        double *resample_t;
        double *newRRI;
        /*
         *  end RRI
         */
        /*
         *  about window FFT
         *  Using ROOT    (fftw3)
         */
        TH2D *spectrum;
        int fft_now;
        int num_of_frequency;
        int windows_num;
        int freq_num;
        double *log_of_spectrum;
        double *frequency_fft;
        double *power_spectrum;
        /*
         *  end FFT
         */
        /*
         * breathing and apnea band
         */
        TGraph *breathing;
        TGraph *apnea;
        double br_fr_start;
        double br_fr_end;
        double ap_fr_start;
        double ap_fr_end;
        /*
         *  end band
         */
        int starttime;
        double *apnea_event;
        TGraph *show_apnea;
        TGraph *ave_apnea;
        //
        vector< vector<double> > apneastart;
        vector< vector<double> > apneaend;
        vector< vector<double> > breathstart;
        vector< vector<double> > breathend;

        vector< vector<double> > apneastart_one;
        vector< vector<double> > apneaend_one;
        vector< vector<double> > breathstart_one;
        vector< vector<double> > breathend_one;
    public:
        SpectrumForRRI(std::string filepath,int freq,int window_size, int overflap);
        /*
         * filepath : the path of rri file
         * freq     : resample frequency
         * windows_size : num of samples included in one window
         * overlap  : num of samples overlap between two windows
         */
        /*
         *  methods for cal new RRI
         */
        void read_data(); // read data from rri file
        void cal_resample(); // resample RRI to the frequency
        void cal_resample_time(); // cal resample time
        int get_num_of_RR_peaks(){
            return num_of_RR_peaks;
        }
        int get_total_num_resample(){return total_num_after_resample;}
        int get_num_of_window(){
            return windows_num;
        }
        double cal_value_at_t(double t, int istart=0);
        double *get_RR_interval(){return R_R_interval;}
        double *get_RR_position(){return R_R_position;}
        double *get_resample_t(){return resample_t;}
        double *get_new_RRI(){return newRRI;}
        /*
         * methods for windows FFT
         */
        void apply_FFT(int i,bool dofft = kTRUE,bool save = kFALSE);
        void cal_spectrum();
        void draw_fft_now();
        void cal_band();
        TGraph * get_breathing_band(){return breathing;}
        TGraph * get_apnea_band(){return apnea;}
        TH2D *get_spectrum(){return spectrum;}
        int get_total_number(){return num_of_frequency;}
        void cal_apnea(std::string filename, char c = ';');
        void cal_around_apnea(std::string filename, char c = ';');
        vector< vector <double> > get_apnea_start_vec(){return apneastart;}
        vector< vector<double> > get_apnea_end_vec(){return apneaend;}
        vector< vector<double> > get_breathing_start_vec(){return breathstart;}
        vector< vector<double> > get_breathing_end_vec(){return breathend;}

        vector< vector <double> > get_apnea_start_one_vec(){return apneastart_one;}
        vector< vector<double> > get_apnea_end_one_vec(){return apneaend_one;}
        vector< vector<double> > get_breathing_start_one_vec(){return breathstart_one;}
        vector< vector<double> > get_breathing_end_one_vec(){return breathend_one;}
        TGraph *get_apnea(){
            return show_apnea;
        }
        TGraph *get_ave_apnea(){
            return ave_apnea;
        }
        virtual ~SpectrumForRRI();
};


double mean_value(vector<double> values);
double mean_value(double *values, int size);
double std_dev(vector<double> values);
double std_dev(double *values, int size);
double std_err(vector<double> values);
double std_dev(double *values, int size);

#endif /* SPECTRUM_FOR_RRI_H */
