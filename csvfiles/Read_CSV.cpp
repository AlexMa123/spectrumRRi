/*
 * =====================================================================================
 *
 *       Filename:  Read_CSV.cpp
 *
 *    Description:  from csv files read the apnea and hyponea data
 *
 *        Version:  1.0
 *        Created:  02/27/2018 15:27:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yaopeng Ma (), sdumyp@126.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

#include "Read_CSV.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::stringstream;
using std::ifstream;


ReadCSV::ReadCSV(string filename, char c)
{
    set_file(filename, c);
}
void ReadCSV::set_file(std::string filename, char c){
    separate = c;
    if(!strArray.empty()){
        strArray.clear();
    }
    inFile.open(filename,std::ios::in);
    string lineStr;
    while(getline(inFile,lineStr)){
        stringstream ss(lineStr);
        string str;
        vector<string> lineArray;
        while(getline(ss,str,separate)){
            str.erase(0,1);
            str.pop_back();
            lineArray.push_back(str);
        }
        strArray.push_back(lineArray);
    }

}

ReadCSV::~ReadCSV(){
    inFile.close();
    strArray.clear();
}

double str2num(string str){
    std::istringstream iss(str);
    double out;
    iss >> out;
    return out;
}

mytime str2time(string str){
    int h, m, s;
    char c;
    std::istringstream iss(str);
    iss >> h;
    iss >> c;
    iss >> m;
    iss >> c;
    iss >> s;
    return mytime(h, m, s);
}

Event::Event(string filename, char c)
    : ReadCSV(filename, c)
{
    int num_of_colume = get_colume();
    for(int i = 1; i < num_of_colume; i++){
        if(get_value(i, 0).find("Apnoe") < 30){
            if(get_value(i, 13) != "-"){
                mytime ti = str2time(get_value(i, 2));
                //cout << i << "th apnea start at time "<<ti.hour <<":" << ti.minuit << ":" << ti.second << std::endl;
                if(ti.hour < 12){
                    ti.hour += 24;
                }
                apnea_starttime.push_back(ti.get_total_second());
                apnea_duration.push_back(str2num(get_value(i,5)));
            }
        }
        else if(get_value(i, 0).find("Hypopnoe") < 30){
            if(get_value(i, 13) != "-"){
                mytime ti = str2time(get_value(i, 2));
                //cout << i << "th hypopnea start at time "<<ti.hour <<":" << ti.minuit << ":" << ti.second << std::endl;
                if(ti.hour < 12)
                    ti.hour += 24;
                hypopnea_starttime.push_back(ti.get_total_second());
                hypopnea_duration.push_back(str2num(get_value(i,5)));
            }
        }
    }
    num_of_apnea = apnea_duration.size();
    num_of_hypopnea = hypopnea_duration.size();
}

double Event::get_apnea_starttime(int i){
    if(i >= num_of_apnea){
        std::cout << "i should less than" << num_of_apnea << std::endl;
        return -1;
    }
    return apnea_starttime[i];
}
double Event::get_apnea_duration(int i){
    if(i >= num_of_apnea){
        std::cout << "i should less than" << num_of_apnea << std::endl;
        return -1;
    }
    return apnea_duration[i];
}
double Event::get_hypopnea_starttime(int i){
    if(i >= num_of_hypopnea){
        std::cout << "i should less than" << num_of_hypopnea << std::endl;
        return -1;
    }
    return hypopnea_starttime[i];
}
double Event::get_hypopnea_duration(int i){
    if(i >= num_of_hypopnea){
        std::cout << "i should less than" << num_of_hypopnea << std::endl;
        return -1;
    }
    return hypopnea_duration[i];
}
