/*
 * =====================================================================================
 *
 *       Filename:  Read_CSV.h
 *
 *    Description:  from csv file read the apnea and hyponea data
 *
 *        Version:  1.0
 *        Created:  02/27/2018 15:28:38
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yaopeng Ma (), sdumyp@126.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#ifndef READ_CSV_H
#define READ_CSV_H

using std::vector;
using std::string;

struct mytime{
    int hour;
    int minuit;
    int second;
    mytime(int h, int m, int s){
        hour   = h;
        minuit = m;
        second = s;
    }
    mytime(const mytime &t){
        hour = t.hour;
        minuit = t.minuit;
        second = t.second;
    }
    int get_total_second(){
        return hour*3600 + minuit*60 + second;
    }
    void print(){
        printf("%2d:%2d:%2d",hour, minuit, second);
    }
    mytime operator + (mytime & t){
        int s, m, h;
        s = second + t.second;
        m = minuit + t.minuit + s/60;
        h = hour   + t.hour   + m/60;
        s = s%60;
        m = m%60;
        return mytime(h, m, s);
    }
    mytime operator - (mytime & t){
        int s, m, h;
        s = second - t.second;
        if( s < 0 ){
            m = minuit - t.minuit - 1;
            s += 60;
        }
        else
            m = minuit - t.minuit;
        if( m < 0 ){
            h = hour - t.hour - 1;
            m += 60;
        }
        else
            h = hour - t.hour;
        return mytime(h, m, s);
    }
};

/**********************************************************************
*                  basic class for reading csv file                  *
**********************************************************************/

class ReadCSV
{
    private:
        std::ifstream inFile;
        vector<vector<string> > strArray;
        char separate;

    public:
        ReadCSV(std::string filename, char c=';');
        void set_file(std::string filename, char c);
        vector<vector<string> > get_whole_file(){return strArray;}
        string get_value(int i, int j){
            if(i < 0 || i >= (int) strArray.size()){
                printf("i should between %d to %d\n", 0,(int)strArray.size()-1);
                return "ERROR!!";
            }
            if(j < 0 || j >= (int) strArray[0].size()){
                printf("j should between %d to %d\n",0,(int)strArray[0].size()-1);
                return "ERROR!!";
            }
            return strArray[i][j];
        }
        int get_colume(){
            return strArray.size();
        }
        int get_line(){
            return strArray[0].size();
        }
        virtual ~ReadCSV();
};
/**********************************************************************
*                    end of class def for ReadCSV                    *
**********************************************************************/

/**********************************************************************
*                           unity function                           *
**********************************************************************/

double str2num(string str);
mytime str2time(string str);
mytime second_to_time(int second);

/**********************************************************************
*             class for reading event file to get apnea              *
**********************************************************************/

class Event: private ReadCSV{
    private:
        vector<double> apnea_starttime;
        vector<double> apnea_duration;
        vector<double> apnea_central_starttime;
        vector<double> apnea_central_duration;
        vector<double> apnea_mixed_starttime;
        vector<double> apnea_mixed_duration;
        vector<double> apnea_obstructive_duration;
        vector<double> apnea_obstructive_starttime;
        vector<double> hypopnea_starttime;
        vector<double> hypopnea_duration;

        int num_of_apnea;
        int num_of_hypopnea;
    public:
        Event(string filename, char c = ';');
        double get_apnea_starttime(int i); // unit is second
        double get_apnea_duration(int i);  // unit is second
        double get_hypopnea_starttime(int i);
        double get_hypopnea_duration(int i);
        int get_number_of_apnea(){
            return num_of_apnea;
        }
        int get_number_of_hypopnea(){
            return num_of_hypopnea;
        }
        vector<double> get_apnea_starttime(){
            return apnea_starttime;
        }
        vector<double> get_apnea_duration(){
            return apnea_duration;
        }
        vector<double> get_hypopnea_starttime(){
            return hypopnea_starttime;
        }
        vector<double> get_hypopnea_duration(){
            return hypopnea_duration;
        }
};

#endif /* READ_CSV_H */
