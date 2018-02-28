/*
 * =====================================================================================
 *
 *       Filename:  rootlogon.C
 *
 *    Description:  root file
 *
 *        Version:  1.0
 *        Created:  02/27/2018 15:49:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yaopeng Ma (), sdumyp@126.com
 *   Organization:  
 *
 * =====================================================================================
 */

{
    gROOT->ProcessLine(".L csvfiles/Read_CSV.cpp+");
    gROOT->ProcessLine(".L windowsfft/spectrum_for_RRI.cpp+");
    gROOT->ProcessLine(".L utility/my_functions.cpp+");
}
