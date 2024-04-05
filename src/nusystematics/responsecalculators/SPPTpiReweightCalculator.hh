#pragma once

#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/simbUtility.hh"

#include <cmath>

namespace nusyst {

inline double GetSPPQ2TemplateReweight(double Q2_GeV2){

  double X = Q2_GeV2;
  if(Q2_GeV2>=3.0) X = 3.000000;

  double this_Q2RW = 1.;
  if( X < 0.025000) this_Q2RW = 1.253255;
  else if( X >= 0.025000 && X < 0.050000) this_Q2RW = 1.589738;
  else if( X >= 0.050000 && X < 0.100000) this_Q2RW = 1.733869;
  else if( X >= 0.100000 && X < 0.200000) this_Q2RW = 1.651728;
  else if( X >= 0.200000 && X < 0.300000) this_Q2RW = 1.659705;
  else if( X >= 0.300000 && X < 0.400000) this_Q2RW = 1.584229;
  else if( X >= 0.400000 && X < 0.500000) this_Q2RW = 1.703793;
  else if( X >= 0.500000 && X < 0.700000) this_Q2RW = 1.475510;
  else if( X >= 0.700000 && X < 1.000000) this_Q2RW = 1.456727;
  else if( X >= 1.000000 && X < 1.300000) this_Q2RW = 1.252215;
  else if( X >= 1.300000 && X < 2.000000) this_Q2RW = 1.048199;
  else if( X >= 2.000000 && X < 3.000000) this_Q2RW = 1.650489;
  else{
    this_Q2RW = 1.650489;
  }

  return this_Q2RW;

}

inline double GetSPPTpiReweight(double Tpi_GeV){

  static double landau_Cutoff = 0.225;

  double this_TpiRW = 1.;
  if(Tpi_GeV<=0.){
    // TODO Including zero or not?
    this_TpiRW = 1.;
  }
  else if(Tpi_GeV<landau_Cutoff){
    // Params for Function = norm * ROOT.TMath.Landau(value, mu, sigma)
    // norm, mpv, width
    static double LandauParams[3] = {6.70797696, 0.12235454, 0.05731087};
    this_TpiRW = LandauParams[0] * TMath::Landau(Tpi_GeV, LandauParams[1], LandauParams[2]);
  }
  else{
    if( landau_Cutoff <= Tpi_GeV && Tpi_GeV < 0.250000 ) this_TpiRW = 0.755932;
    else if( 0.250000 <= Tpi_GeV && Tpi_GeV < 0.275000 ) this_TpiRW = 0.638574;
    else if( 0.275000 <= Tpi_GeV && Tpi_GeV < 0.300000 ) this_TpiRW = 0.493987;
    else if( 0.300000 <= Tpi_GeV && Tpi_GeV < 0.325000 ) this_TpiRW = 0.391947;
    else if( 0.325000 <= Tpi_GeV && Tpi_GeV < 0.350000 ) this_TpiRW = 0.323265;
    else if( 0.350000 <= Tpi_GeV && Tpi_GeV < 0.400000 ) this_TpiRW = 0.452765;
    else if( 0.400000 <= Tpi_GeV && Tpi_GeV < 0.500000 ) this_TpiRW = 0.594541;
    else if( 0.500000 <= Tpi_GeV && Tpi_GeV < 0.700000 ) this_TpiRW = 0.768459;
    else if( 0.700000 <= Tpi_GeV && Tpi_GeV < 1.000000 ) this_TpiRW = 0.658024;
    else if( 1.000000 <= Tpi_GeV && Tpi_GeV < 2.000000 ) this_TpiRW = 0.873622;
    else this_TpiRW = 0.873622;
  }

  return this_TpiRW;

}

// CV correction
inline double GetSPPTpiCVCorrection(double Q2_GeV2, double Tpi_GeV){

  double this_Q2RW = GetSPPQ2TemplateReweight(Q2_GeV2);
  double this_TpiRW = GetSPPTpiReweight(Tpi_GeV);

  return this_Q2RW * this_TpiRW;

}

inline double GetSPPTpiCorrectionRW(double Q2_GeV2, double Tpi_GeV, double parameter_value){
  
  double CVCorr = GetSPPTpiCVCorrection(Q2_GeV2, Tpi_GeV);

  // 1/CVCorr is the correction back to nominal = 1sigma 
  double oneSigRW = 1./CVCorr;
  // Size of the one-sigma uncertainty obtained by subtracting 1
  double oneSigUnc = oneSigRW-1.;

  double this_rw = 1. + parameter_value * oneSigUnc;

  return this_rw;

}

};

