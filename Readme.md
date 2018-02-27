### Files:

1. csvfiles/Read_CSV : include classed for reading csv files and get apnea events from csv files:

   ```cpp
   struct mytime(int h,int m,int s); // h: hour, m: minuit, s: second
   int mytime::get_total_second(); // return the time in seconds
   double str2num(string str); // string to double
   mytime str2time(string str);// string to mytime struct
   ```

   ```cpp
   class ReadCSV(std::string filename, char c=';'); // filename : path of csvfiles ， c: delimiters of csv files
   string ReadCSV::get_value(int i, int j); // return the data in ith line, jth column
   int ReadCSV::get_colume(); // return the num of colume
   int ReadCSV::get_line();   // return the num of line
   ```

   ```cpp
   class Event(string filename, char c=';'): private ReadCSV;
   double Event::get_*_starttime(int i); // *(apnea or hypopnea), i start from 0
   double Event::get_*_duration(int i);
   int get_number_of_*();
   vector<double> get_*_starttime();
   vector<double> get_*_duration(); // get the vector of starttime or duration.
   ```

2. windowsfft/spectrum_for_RRI calculate windows fft of RRI file

   ```cpp
   class SpectrumForRRI(string filepath, int freq, int windows_size, int overlap);
   // filepath	    : the path of RRI file, 
   // freq    	    : resample frequency
   // windows_size : num of samples in one window
   // overlap      : num of samples overlap between two windows

   void SpectrumForRRI::apply_FFT(int i, bool dofft = kTRUE, bool save=kFALSE);
   // do fft in ith window, 
   // if dofft = true, fft will apply to the windows, or it won't
   // if save = true, the value of resample rri in the ith windows will save into a .dat file
   void SpectrumForRRI::draw_fft_now(); 
   void SpectrumForRRI::cal_band(); // cal the ratio of breathing band and apnea band
   void SpectrumForRRI::cal_spectrum(); 
   // for get value:
   void SpectrumForRRI::cal_apnea(string filename, char c = ";");
   // filename : path of csv
   int SpectrumForRRI::get_num_of_RR_peaks(); // return num of peaks included in one file
   int SpectrumForRRI::get_num_of_window();   // num of windows 
   int SpectrumForRRI::get_total_num_resample();
   double *get_RR_interval();
   double *get_RR_position(); // sum of RRI
   double *get_resample_t();
   double *get_new_RRI();
   TGraph *get_breathing_band();
   TGraph *get_apnea_band();
   TGraph *get_apnea();
   TGraph *get_ave_apnea();
   TH2D* get_spectrum();
   ```

3. rootlogon.C : root will load this file when it started.
