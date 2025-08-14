#ifndef weightAMUvariations_h
#define weightAMUvariations_h

#include <fstream>  //ifstream
#include <iostream> //cout
#include <cmath>

// Removed all dependencies on root
// for example TMath::stuff

//#include "math.h"
//#include "assert.h"

// The input arrays will be uniform and known at compile time
// If needed later, this can be refactored to be more flexible
#define XBINS 19
#define Q2BINS 6
#define NWEIGHTS 10


//For Compatibility with ROOT compiler
//uncomment the following:
//#define ROOT

// My Callan Gross relation is F1 = F2 / 2x and does not have TMC (1 + 4 MN x2 / Q2)/(1 + r) with r ~ FL
// There are two bits there, the TMC bit and the r bit.  See QPMDISStrucFuncBase.cxx

/*!
 *  Update this text with correct citations
 *  This is based on code by 
 *  Heidi Schellman (Oregon State) and Rik Gran (Minnesota Duluth)
 *  originally to read in a file and interpolate RPA uncertainties for QE
 *  for use in MINERvA experiment analysis
 *  must compile with the ROOT libraries
 *  g++ `root-config --glibs --cflags` -O3 weightAMUDIStest.cxx -o weightAMUDIStest
 *    where the text.cxx code instantiates and then calls this class for some given x, y, Enu, nu/antinu
 *
 *    This version instead gives variations of the DIS cross section
 *    Especially in the 1.0 < Q2 < 4.0 GeV2 SIS region but also Q2 < 1.0 GeV2 multiquark region.
 *    To be cute, I renamed it AMUvariations, homage to Tom Waits' 1999 album Mule Variations
 *
 *    As of this writing, the code returns the full effect of a single feature of the AMU model
 *    The calling code should decide what fraction to take as a systematic (0.25 for 1 sigma?)
 *    And whether to separate and/or uncorrelate the Q2 < 1.0 part which is really outside the model
 *    even though the model does calculate stuff there.
 *  
 *    Beyond the edge of the Q2 and x space (currently 0.05 < x < 0.95, 0.5 < Q2 < 20) the code
 *    uses the weight from a nearby point within the space rather than extrapolating or returning no weight.
 *    At the edge of this 2D space, there may be other shortcomings in the model to account for,
 *    which the user or future developers should consider with specific use-cases in mind.
 *
 *    Also, the inputs are isoscalar i.e. (p + n)/2 .  The ratios especially are certainly meaningful
 *    for all nuclei except hydrogen, especially if applied to a correct non-isoscalar model in the first place.    
 *
 *    The underlying model is from the Aligarh Muslim University group
 *    courtesy of Huma Haider and students
 *    see (public) Zaidi:2019mfd Zaidi, F., H. Haider, M. Sajjad Athar, S.K. Singh, and I. Ruiz Simo
 *    Nucleon and nuclear structure functions with nonperturbative and higher order perturbative QCD effects
 *    Phys. Rev. D 99 (2019) p. 093011
 *
 *    When this expands to nuclear effects I guess citation will include PRD 101 (2020) p.033001 
 *
 */


// NOTE UNITS ARE GEV in the calculation
// make sure you convert MeV to GeV before calling these functions

namespace nusyst {

class weightAMUvariations {
public:
  //Constructor: Read in params from a filename
  //weightAMUvariations(const TString filename1, const TString filename2, const TString filename3) { read(filename1,filename2,filename3); }    //Read in params from file
  weightAMUvariations(const char *filename1, const char *filename2, const char *filename3) { read(filename1,filename2,filename3); }    //Read in params from file


  double ar1LO[XBINS][Q2BINS];
  double ar2LO[XBINS][Q2BINS];
  double ar3LO[XBINS][Q2BINS];
  double ar1NLO[XBINS][Q2BINS];
  double ar2NLO[XBINS][Q2BINS];
  double ar3NLO[XBINS][Q2BINS];
  double ar1NLOTMC[XBINS][Q2BINS];
  double ar2NLOTMC[XBINS][Q2BINS];
  double ar3NLOTMC[XBINS][Q2BINS];
  double ar1NLOTMCHT[XBINS][Q2BINS];
  double ar2NLOTMCHT[XBINS][Q2BINS];
  double ar3NLOTMCHT[XBINS][Q2BINS];
  double ar1NLOTMCHTC12[XBINS][Q2BINS];
  double ar2NLOTMCHTC12[XBINS][Q2BINS];
  double ar3NLOTMCHTC12[XBINS][Q2BINS];
  double ar1NLOTMCHTCa40[XBINS][Q2BINS];
  double ar2NLOTMCHTCa40[XBINS][Q2BINS];
  double ar3NLOTMCHTCa40[XBINS][Q2BINS];
  double ar1GenieDefaultFree[XBINS][Q2BINS];
  double ar2GenieDefaultFree[XBINS][Q2BINS];
  double ar3GenieDefaultFree[XBINS][Q2BINS];
  double ar1GenieBugfixFree[XBINS][Q2BINS];
  double ar2GenieBugfixFree[XBINS][Q2BINS];
  double ar3GenieBugfixFree[XBINS][Q2BINS];
  double ar1GenieDefaultNuc[XBINS][Q2BINS];
  double ar2GenieDefaultNuc[XBINS][Q2BINS];
  double ar3GenieDefaultNuc[XBINS][Q2BINS];


  // MINERvA holds kinematics in MeV, but all these functions require GeV
  // So make sure you pass them in GeV.
  double getWeight(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double * weights);
  double getAbsolute(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double * weights);
  double getF1(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double * weights);
  double getF2(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double * weights);
  double getF3(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double * weights);
  double getF1withQ2(const double x, const double Q2, const double Enu, const double mlep, const int isAntinu, double * weights);
  double getF2withQ2(const double x, const double Q2, const double Enu, const double mlep, const int isAntinu, double * weights);
  double getF3withQ2(const double x, const double Q2, const double Enu, const double mlep, const int isAntinu, double * weights);

  double getIntegral(const double Enu, const double mlep, const int isAntinu, double * weights, double W = -1.0);


  //in GeV

  //Initializer
  //void read(const TString filename1, const TString filename2, const TString filename3);
  void read(const char *filename1, const char *filename2, const char *filename3);
  
  
  private:

  // These are hard coded and the input file needs to be known in advance
  // Thats fine for now and makes the code easy to read.
  double xArray[XBINS] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
  0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};

  // y bins are Q2 bins, and are variable sized, so bin center is tricksy.
  // I am using lower case q2 but its really the positive Q2 throughout
  double q2Array[Q2BINS] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0};


  // These are private members of the class
  // because I am computing them once for every x, y, Q2
  // and reusing them for all the ratios that get calculated.
  double mXforStructureFunction;
  double mQ2forStructureFunction;
  double mFront;
  double mF1term;
  double mF1termCG;
  double mF2term;
  double mF3term;
  double mF4term;
  double mF5term;
  double mSign;
  int mAbsolute;

  double MN = 0.5 * (0.93827208943 + 0.93956542052);

  // Also, the interpolator is going to figure out four xy points (really x, Q2 points)
  // and since the input arrays are uniform, these can be reused too.
  // might only need two array indexes and get the four points from those
  int mxarrayindex;  // this index and +1 gives mx1 and mx2
  int mqqarrayindex;  // this index and +1 gives mqq1 and mqq2
  double mx1;
  double mx2;
  double mqq1;  // is the "y" value in the interpolation, but really Q2 in physics
  double mqq2;  // is the "y" value in the interpolation, but really Q2 in physics

  
  double getWeightInternal(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double *weights, int absolute = 0);

  // these are going to use the pre-computed front and kinematic term factors stored as members of the class.
  double getOneCrossSection(double arrayf1[][Q2BINS], double arrayf2[][Q2BINS], double arrayf3[][Q2BINS], double x = -1.0);
  double getOneCrossSection(double F1, double F2, double F3, double x = -1.0);

  // the two bits of interpolate2D are separated
  // so that I can take the ratio first then interpolate, which might be more precise
  // or legacy code will simply interpolate first then take the ratio
  double getInterpolatedValueFromFour(double x, double y, double q11, double q12, double q21, double q22);
  double getInterpolatedValueFromFour(double x, double y, double array[4]);

  //double getFourValuesForInterpolation(double x, double y, const double array[][Q2BINS], double q11, double q12, double q21, double q22);
  double getFourValuesForInterpolation(double x, double y, const double array[][Q2BINS], double fourval[4]);
  double getFourXSForInterpolation(const double array1[][Q2BINS], const double array2[][Q2BINS], const double array3[][Q2BINS], double fourXS[4], double x = -1.0);

  double interpolate2D(double x, double y, const double array[][Q2BINS]);

};

#endif

  //void weightAMUvariations::read(const TString filename1, const TString filename2, const TString filename3){
  inline void weightAMUvariations::read(const char *filename1, const char *filename2, const char *filename3){
    //Read in the params doubles from a file
  //argument: valid filenames

  //fRPAratio = TFile::Open(f,"READONLY");

  std::ifstream in1;
  in1.open(filename1); //in1.open(filename1.Data());
  
  std::ifstream in2;
  in2.open(filename2); //in2.open(filename2.Data());

  std::ifstream in3;
  in3.open(filename3); // in3.open(filename3.Data());

  if(!in1 || !in2 || !in3){
    std::cout << "File " << filename1 << " or " << filename2 << " or " << filename3 << " not found.  exiting." << std::endl;
    std::cerr << "File " << filename1 << " or " << filename2 << " or " << filename3 << " not found.  exiting." << std::endl;
    return;
  }

  // The original implementation tried to stuff this into a TH2D
  // But in the end the binning of the model runs won't be a perfect match
  // To the need to put each of these points at bin center.
  // So switched to arrays after all.

  // These are hard coded and the input file needs to be known in advance
  // Thats fine for now and makes the code easy to read.
  // They are private members of the class, so copied here just for reference
  // Don't need std::vector functionality, just need a fixed array and size but didn't change it back.
  //
  // double xArray = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 
  //    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9};
    
  // y bins are Q2 bins, and are variable sized, so bin center is tricksy.
  // I am using lower case q2 but its really the positive Q2 throughout
  // double q2Array = {1, 2, 5, 10, 20};

  // use as an array not a histogram.
  // double f1array[XBINS][Q2BINS];
  // these are members of the class instantiated at compile time.
    
  while(1){
    double q2, x, a, b, c, d, f, g, h, k, m;
    in1 >> q2 >> x >> a >> b >> c >> d >> f >> g >> h >> k >> m;
    if(!in1.good()) break;
    //std::cout << q2 << " " << x << " " << a << " " << b << " " << c << " " << d << std::endl;

    // I'm expecting near exact match to the bins
    // add some checking for the possibility the right point is not found.
    bool notfound = true;
    int ix = 0;
    double tolerance = 0.00001;
    for(ix = 0; ix <= XBINS; ix++){
      if( std::abs(x - xArray[ix] ) <= tolerance){
        notfound = false;
        break;
      }
      // ix is set
    }
    if(notfound)std::cerr << " Problem f1 data file input x doesnt match the internal x bins" << std::endl;
    notfound = true;
    int iq2 = 0;
    for(iq2 = 0; iq2 <= Q2BINS; iq2++){
      if( std::abs(q2 - q2Array[iq2] ) <= tolerance){
        notfound = false;
        break;
      }
      // iq2 is set
    }
    if(notfound)std::cerr << " Problem f1 data file input Q2 doesnt match the internal Q2 bins" << std::endl;

    ar1LO[ix][iq2] = a;
    ar1NLO[ix][iq2] = b;
    ar1NLOTMC[ix][iq2] = c;
    ar1NLOTMCHT[ix][iq2] = d;
    ar1NLOTMCHTC12[ix][iq2] = f;
    ar1NLOTMCHTCa40[ix][iq2] = g;
    ar1GenieDefaultFree[ix][iq2] = h;
    ar1GenieBugfixFree[ix][iq2] = k;
    ar1GenieDefaultNuc[ix][iq2] = m;

    //std::cout << "ij " << ix << " " << iq2 << " " << x << " " << q2 << " " << a << std::endl;

  }

  // test does the array look right ?
  if(0){
    for(int ixx=0; ixx<19; ixx++){
      std::cout << " ix " << ixx << " " << xArray[ixx] << " array " 
      << ar1LO[ixx][0] << " " << ar1LO[ixx][1] << " " << ar1LO[ixx][2] 
      << " " << ar1LO[ixx][3] << " " << ar1LO[ixx][4] << std::endl;
      //if(ixx >= 20)break;
    }
  }

  while(1){
    double q2, x, a, b, c, d, f, g, h, k, m;
    in2 >> q2 >> x >> a >> b >> c >> d >> f >> g >> h >> k >> m;
    if(!in2.good()) break;
    //std::cout << q2 << " " << x << " " << a << " " << b << " " << c << " " << d << std::endl;

    // I'm expecting near exact match to the bins
    // add some checking for the possibility the right point is not found.
    bool notfound = true;
    int ix = 0;
    double tolerance = 0.00001;
    for(ix = 0; ix <= XBINS; ix++){
      if( std::abs(x - xArray[ix] ) <= tolerance){
        notfound = false;
        break;
      }
      // ix is set
    }
    if(notfound)std::cerr << " Problem f2 data file input x doesnt match the internal x bins" << std::endl;
    notfound = true;
    int iq2 = 0;
    for(iq2 = 0; iq2 <= Q2BINS; iq2++){
      if( std::abs(q2 - q2Array[iq2] ) <= tolerance){
        notfound = false;
        break;
      }
      // iq2 is set
    }
    if(notfound)std::cerr << " Problem f2 data file input Q2 doesnt match the internal Q2 bins" << std::endl;

    ar2LO[ix][iq2] = a;
    ar2NLO[ix][iq2] = b;
    ar2NLOTMC[ix][iq2] = c;
    ar2NLOTMCHT[ix][iq2] = d;
    ar2NLOTMCHTC12[ix][iq2] = f;
    ar2NLOTMCHTCa40[ix][iq2] = g;
    ar2GenieDefaultFree[ix][iq2] = h;
    ar2GenieBugfixFree[ix][iq2] = k;
    ar2GenieDefaultNuc[ix][iq2] = m;

    //std::cout << "ij " << ix << " " << iq2 << " " << x << " " << q2 << " " << a << std::endl;

  }

  while(1){
    double q2, x, a, b, c, d, f, g, h, k, m;
    in3 >> q2 >> x >> a >> b >> c >> d >> f >> g >> h >> k >> m;
    if(!in3.good()) break;
    //std::cout << q2 << " " << x << " " << a << " " << b << " " << c << " " << d << std::endl;
    bool notfound = true;
    // I'm expecting near exact match to the bins
    // add some checking for the possibility the right point is not found.
    int ix = 0;
    double tolerance = 0.00001;
    for(ix = 0; ix <= XBINS; ix++){
      if( std::abs(x - xArray[ix] ) <= tolerance){
        notfound = false;
        break;
      }
      // ix is set
    }
    if(notfound)std::cerr << " Problem f3 data file input x doesnt match the internal x bins" << std::endl;
    notfound = true;
    int iq2 = 0;
    for(iq2 = 0; iq2 <= Q2BINS; iq2++){
      if( std::abs(q2 - q2Array[iq2] ) <= tolerance){
        notfound = false;
        break;
      }
      // iq2 is set
    }
    if(notfound)std::cerr << " Problem f3 data file input Q2 doesnt match the internal Q2 bins" << std::endl;

    ar3LO[ix][iq2] = a;
    ar3NLO[ix][iq2] = b;
    ar3NLOTMC[ix][iq2] = c;
    ar3NLOTMCHT[ix][iq2] = d;
    ar3NLOTMCHTC12[ix][iq2] = f;
    ar3NLOTMCHTCa40[ix][iq2] = g;
    ar3GenieDefaultFree[ix][iq2] = h;
    ar3GenieBugfixFree[ix][iq2] = k;
    ar3GenieDefaultNuc[ix][iq2] = m;

    //std::cout << "ij " << ix << " " << iq2 << " " << x << " " << q2 << " " << f << std::endl;

  }

  in1.close();
  in2.close();
  in3.close();  // probably not necessary.

  mSign = 1.0;
  mFront = 0.0;
  mF1term = 0.0;
  mF1termCG = 0.0;
  mF2term = 0.0;
  mF3term = 0.0;
  mF4term = 0.0;
  mF5term = 0.0;
  mXforStructureFunction = -1.0;
  mQ2forStructureFunction = -1.0;
  //Needed to test whether recalculation is needed.
  //mInputX = -1.0;
  //mInputY = -1.0;
  //mInputEnu = -1.0;
  //mInputMlep = -1.0;
  //mInputNuAntinu = -1.0;
  
}


inline double weightAMUvariations::getWeightInternal(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double *weights, int absolute){

  //should this be x and Q2 or x and y ?

  // Put something here to store the previous mc_x and mc_y calculation
  // and simply return the same weight without recalculating.
  // the calling program is responsible for making sure the array is the right size
  if (weights != 0){
    for (int i = 0; i < NWEIGHTS; i++){
      weights[i] = 1.0;
    }
  }
  
  mAbsolute = absolute;

  // Actually I don't want the above, thats just for playing around.
  // I want to return a set of weights for each pair.
  //double getOneCrossSectionxy(double x, double y, double Enu, double mlep, int isAntinu, *TH2D hF2, *TH2D hF3){
  //q2,x, f2(free), f2(NLO), f2(NLO+TMC), f2(NLO+TMC+HT)
  //q2,x, f3(free), f3(NLO), f3(NLO+TMC), f3(NLO+TMC+HT)

  // Pre-calculate all the kinematic terms.  Hopefully the compiler optimizes.
  // Could save them as a private member of the class if I want,
  // then in a multi-universe environment they might be saved one to next.

  // Mnucleon = is a member of the class.
  //double mMNforAMU = 0.5 * (0.93827208943 + 0.93956542052);
  double mlep2 = mlep*mlep;

  // need to convert xy into Q2
  // I'm using the definition of x and y in terms of nu, Ehad, and Enu and Q2.
  // Could go through Q2 = 2 Enu (Emu + pmu cos(theta)) + mmu2 give or take a sign mistake.
  // Dont use 4 Enu Emu sin2(theta/2) which is only true for massless neutrinos and leptons.
  double Q2 = 2.0 * x * MN * (y * Enu);

  if(0){
    // not needed, but what is q3 ?   Turn this on for debugging or orienting.
    double Ehad = y * Enu;   // is basically energy transfer.
    double q3 = std::sqrt(Ehad*Ehad + Q2);
    std::cout << "x " << x << " " << y << " " << Enu << " " << mlep << " " << MN << " Q2 " << Q2 << " q3q0 " << q3 << " " << Ehad << std::endl;
  }

  //
  //  Bounds checking on the TH2D interpolation
  //  Two options, strict bounds and return something that will become no weight 1.0 outside bounds
  //  Catch the bounds and return the weight from some point back at the edge of the bounds
  //  The second option is acive, the first option is commented out.
  //  The second option assumes the TH2D Interpolate can extrapolate a little bit.
  //  and applies only to the histogram interpolate (x,Q2), not to other places where x,Q2 are used
  //
  //  Note, the calling routine needs to check that this was generated from a DIS model
  //  this function will return weights for QE and Delta kinematics, if asked.
  //  In fact, probably wise to separate the weights for Q2 < 1 and Q2 > 1
  //  uncorrelated or only partially correlated.

  //  if the inputs are ridiculous, sent back zeros.
  if(x < 0 || x > 1 || y < 0 || y > 1){
    for(int i=0; i<NWEIGHTS; i++)weights[i] = 0.0;
    return 0;
  }

  // Class data members, don't instantiate new ones.
  mQ2forStructureFunction = Q2;
  mXforStructureFunction = x;
  
  // And alternative to declaring it simply doesn't work is to force it to take the nearest reasonable value.

  if(Q2 <= 0.40)mQ2forStructureFunction = 0.40;
  if(Q2 >= 25.0)mQ2forStructureFunction = 25.0;
  if(x >= 0.960)mXforStructureFunction = 0.960;
  if(x <= 0.040)mXforStructureFunction = 0.040;
  
  // Front Matter GF2 M / pi
  mFront = Enu;  // Start without the constants and the 1/(1+Q2/W2)

  // if I am really only taking ratios, can stop here.
  // but the compiler should be able to optimize all the next lines
  // so not much harm in keeping them in.

  double GFermi = 1.166387e-5;  // actually over hbarc3 
  double hbarc2 = 0.1973269804*0.1973269804; // in GeV2 fm2
  double cm2overfm2 = 1.0E-26;
  double pi = 3.14159265359;
  double Wmass2 = 80.379*80.379;  // GeV
  double frontNOEnu = GFermi * GFermi * MN * hbarc2 * cm2overfm2 / pi;
  // I wish I remember why there is no Cabibbo factor here.
  double frontdenom = 1.0 / (1.0 + Q2/Wmass2);  // square it in next line.
  frontNOEnu *= frontdenom * frontdenom;
  // if I am only taking ratios, can skip the math above.
  // this next line is optional if I'm taking ratios and saves a couple optimizable calculations.
  // but turn it back on if I'm comparing absolute cross sections
  mFront *= frontNOEnu;
  // Units will be cm2, I hope.

  //Incidentally, this part of the formula doesn't look like Mark Thompsons version.
  //Its got an extra (MN/Enu)/(1+Q2/MW2)  Not sure why.  Hallsie's version has this term.

  //double F1 = hF1->Interpolate(xforStructureFunction,Q2forStructureFunction);
  mF1term = (y * y * x + (mlep2 * y * 0.5 / (Enu * MN)));
  mF1termCG = mF1term * (1.0 + ( 4.0 * MN * MN * x * x / Q2));

  //double F1F2term = F1term * (1.0 + 4.0 * MN * MN * x * x / Q2);
  //F1term *= ( (1.0 + 4.0 * MN * MN * x * x / Q2) * F2 - FL );
  // document more about the meaning of the F1, FL, and Whitlow term.

  //double F2 = hF2->Interpolate(xforStructureFunction,Q2forStructureFunction);
  mF2term = ( (1.0 - mlep2/(4.0*Enu*Enu)) - y * (1.0 + (0.5 * MN * x / Enu)));
  //F2term *= F2;

  //double F3 = hF3->Interpolate(xforStructureFunction,Q2forStructureFunction);
  mF3term = ( x * y * (1.0 - 0.5 * y) - mlep2 * y * 0.25 / (Enu * MN));
  //F3term *= F3;

  // there are additional terms that are almost always tiny as they are for QE and Delta
  // now I"m reading the formula out of Hallsie's paper PRD 108 113010
  // in the first case, it goes as mlep2/Enu2 in the second case as mlep2/Enu
  //double F4 = not available
  //mF4term = mlep2 * (mlep2 + Q2) / (4.0 * Enu * Enu * MN * MN * x);  // from Reno and Huma
  //mF4term = mlep2 * mlep2 / (4.0 * Enu * Enu * MN * MN * x) + mlep2 * y / (2.0 * Enu * MN) // expand out Reno Q2
  //mF4term = mlep2 * mlep2 / (4.0 * Enu * Enu * MN * MN) + mlep2 * x * y / (2.0 * Enu * MN)  // from Genie.  Same?  No
  // Reno and GENIE differ by a factor of x in this term, unless the convention is the x is in F4 itself, like xF4 .
  // For GENIE and BY it doesn't matter, F4 is zero.   I am looking at the code in QPMDISPXSec.cxx

  //double F5 = not available
  mF5term = -1.0 * mlep2 / (Enu * MN); 
  // This version assumes F5 = F2 / 2x or a similar convention as used in Jeong and Reno, also Albright and Jarlskog
  // GENIE sets F5 = F2 / x but puts the factor of 2 in the front term.   Just a convention.  Watch out.
  
  // In GENIE the front matter is the same as is the terms in front of F1, F2, F3
  // GENIE also applies a scaling based on number of P and N to get xs per nucleus, not per isoscalar nucleon
  // Then a factor of fScale "to reach well known asymmptotic value" a fudge factor, what value does it have?
  // It is simply "DIS-XSecScale"  in most CommonParam.xml its commented out with value 1.032
  // in the JTV tunes it is 1.062213 or possibly 1.019031
  // in QPMDISPXSec.xml it is active and 1.032 .  I guess this value is the default and the JTV tune overrides it.

  // One of these changes sign for anti-neutrino
  mSign = 1.0;
  if(isAntinu) mSign = -1.0;

  // end calculating things that are the same for all cross sections at this x, Q2

  // end calculating things that are the same for all cross sections at this x, Q2

  double XSLO = getOneCrossSection(ar1LO,ar2LO,ar3LO);
  double XSNLO = getOneCrossSection(ar1NLO,ar2NLO,ar3NLO);
  double XSNLOTMC = getOneCrossSection(ar1NLOTMC,ar2NLOTMC,ar3NLOTMC);
  double XSNLOTMCHT = getOneCrossSection(ar1NLOTMCHT,ar2NLOTMCHT,ar3NLOTMCHT);
  double XSNLOTMCHTCG = getOneCrossSection(ar1NLOTMCHT,ar2NLOTMCHT,ar3NLOTMCHT,x);
  // to use the Callan-Gross relation, put in a valid value for x
  //double XSNLOTMCHTCG = getOneCrossSection(h2NLOTMCHT,h2NLOTMCHT,h3NLOTMCHT,x);
  double XSNLOTMCHTC12 = getOneCrossSection(ar1NLOTMCHTC12,ar2NLOTMCHTC12,ar3NLOTMCHTC12);
  double XSNLOTMCHTCa40 = getOneCrossSection(ar1NLOTMCHTCa40,ar2NLOTMCHTCa40,ar3NLOTMCHTCa40);
  double XSGenieDefaultFree = getOneCrossSection(ar1GenieDefaultFree,ar2GenieDefaultFree,ar3GenieDefaultFree);
  double XSGenieBugfixFree = getOneCrossSection(ar1GenieBugfixFree,ar2GenieBugfixFree,ar3GenieBugfixFree);
  double XSGenieDefaultNuc = getOneCrossSection(ar1GenieDefaultNuc,ar2GenieDefaultNuc,ar3GenieDefaultNuc);
  
  // overload the interface so that instead of returning the cross section,
  // these functions return F1, F2, or xF3 .
  // that actually goes in getOneCrossSection, so it needs to see the value of absolute


  if(absolute){
    //std::cout << "is absolute " << XSLO << " " << XSNLO << " " << XSNLOTMC << std::endl;
    weights[0] = XSLO;
    weights[1] = XSNLO;
    weights[2] = XSNLOTMC;
    weights[3] = XSNLOTMCHT;
    weights[4] = XSNLOTMCHTCG;
    weights[5] = XSNLOTMCHTC12;
    weights[6] = XSNLOTMCHTCa40;
    weights[7] = XSGenieDefaultFree;
    weights[8] = XSGenieBugfixFree;
    weights[9] = XSGenieDefaultNuc;

  
  } else {

    // now form the ratios one at a time
    double ratioNLO = XSNLO/XSLO;
    double ratioTMC = XSNLOTMC/XSNLO;
    double ratioHT = XSNLOTMCHT/XSNLOTMC;
    double ratioCG = XSNLOTMCHTCG/XSNLOTMCHT;
    double ratioC12 = XSNLOTMCHTC12/XSNLOTMCHT;
    double ratioCa40 = XSNLOTMCHTCa40/XSNLOTMCHTC12;
    double ratioGenieAmu = XSGenieDefaultFree/XSNLOTMCHT;
    double ratioGenieBug = XSGenieBugfixFree/XSGenieDefaultFree;
    double ratioGenieNuc = XSGenieDefaultNuc/XSGenieDefaultFree;

    bool truncateHighWeights = false;
    double truncateThreshold = 10.0;
    if(truncateHighWeights){
      // it might be better to have a user choose/apply these, so turned off.
      // this is needed because sometimes the cross section goes very small or negative
      // and its set to 1E-43 cm2 .  Maybe 1E-43 is too small ?
      // If that is in the denominator but something else is in the numerator
      // it will lead to very high weights.
      // This is only an effect in anti-neutrino and always high y or high x or both
      // if they are few and random for some sample, truncate will be ok.
      // if instead they are central to the sample, then more careful treatment
      // to deal with weighting a model from something with missing kinematics
      // the user probably needs to know something like this is happening for their sample
      // because the resulting weighted events will always underpredict the target model.
      if(ratioNLO > truncateThreshold)ratioNLO = truncateThreshold;
      if(ratioTMC > truncateThreshold)ratioTMC = truncateThreshold;
      if(ratioHT > truncateThreshold)ratioHT = truncateThreshold;
      if(ratioCG > truncateThreshold)ratioCG = truncateThreshold;
      if(ratioC12 > truncateThreshold)ratioC12 = truncateThreshold;
      if(ratioCa40 > truncateThreshold)ratioCa40 = truncateThreshold;
      if(ratioGenieAmu > truncateThreshold)ratioGenieAmu = truncateThreshold;
      if(ratioGenieBug > truncateThreshold)ratioGenieBug = truncateThreshold;
      if(ratioGenieNuc > truncateThreshold)ratioGenieNuc = truncateThreshold;
    }

    //std::cout << "results " << ratioNLO << " " << ratioTMC << " " << ratioHT << " " << ratioCG << std::endl;

    weights[0] = 1.0;
    weights[1] = ratioNLO;
    weights[2] = ratioTMC;
    weights[3] = ratioHT;
    weights[4] = ratioCG;
    weights[5] = ratioC12;
    weights[6] = ratioCa40;
    weights[7] = ratioGenieAmu;
    weights[8] = ratioGenieBug;
    weights[9] = ratioGenieNuc;

  }



  /*
  if(0){
    // This code has a problem, so its commented out.
    // It was intended to take ratios first then interpolate.
    // there is another problem, the absolute cross section isn't different
    // but I think it should be and I don't know what the problem is.
    // is it because if I use the CV x, y, Enu for the four corners before interpolation
    // it cancels out in the interpolation ?  Even despite the presence of a sum?  Maybe.
    // To check, I need to print out the four corners, see they are different
    // but then confirm the resulting interpolated answer is the same.
    // assuming I am about to do that, explicitly test the extrapolated behavior too.
    // 
    // before doing all of these, make sure to set mx1 mx2 mqq1 mqq2 to -1.0
    // these are members that will take on the same values for all calculations.
    // have not done so yet, but can move this set and test elsewhere
    // so it would be preserved between calls to this function
    // (like applying the same CV to many universes)
    mx1 = -1.0; mx2 = -1.0; mqq1 = -1.0; mqq2 = -1.0;

    // to pipe these things around without crazy amounts of cut and paste, store these in an array.
    double dummy = 0.0;
    // don't need to pass x,y,Q2, its going to use class member variables.

    double fourXSLO[4] = {0.0,0.0,0.0,0.0};  dummy = getFourXSForInterpolation(ar1LO,ar2LO,ar3LO,fourXSLO);
    double fourXSNLO[4] = {0.0,0.0,0.0,0.0};  dummy = getFourXSForInterpolation(ar1NLO,ar2NLO,ar3NLO,fourXSNLO);
    double fourXSNLOTMC[4] = {0.0,0.0,0.0,0.0};  dummy = getFourXSForInterpolation(ar1NLOTMC,ar2NLOTMC,ar3NLOTMC,fourXSNLOTMC);
    double fourXSNLOTMCHT[4] = {0.0,0.0,0.0,0.0};  dummy = getFourXSForInterpolation(ar1NLOTMCHT,ar2NLOTMCHT,ar3NLOTMCHT,fourXSNLOTMCHT);
    double fourXSNLOTMCHTCG[4] = {0.0,0.0,0.0,0.0};  dummy = getFourXSForInterpolation(ar1NLOTMCHT,ar2NLOTMCHT,ar3NLOTMCHT,fourXSNLOTMCHTCG,x);
    double fourXSNLOTMCHTC12[4] = {0.0,0.0,0.0,0.0};  dummy = getFourXSForInterpolation(ar1NLOTMCHTC12,ar2NLOTMCHTC12,ar3NLOTMCHTC12,fourXSNLOTMCHTC12);
    double fourXSNLOTMCHTCa40[4] = {0.0,0.0,0.0,0.0};  dummy = getFourXSForInterpolation(ar1NLOTMCHTCa40,ar2NLOTMCHTCa40,ar3NLOTMCHTCa40,fourXSNLOTMCHTCa40);

    double fourXSGenieDefaultFree[4] = {0.0,0.0,0.0,0.0}; dummy = getFourXSForInterpolation(ar1GenieDefaultFree,ar2GenieDefaultFree,ar3GenieDefaultFree,fourXSGenieDefaultFree);
    double fourXSGenieBugfixFree[4] = {0.0,0.0,0.0,0.0}; dummy = getFourXSForInterpolation(ar1GenieBugfixFree,ar2GenieBugfixFree,ar3GenieBugfixFree,fourXSGenieBugfixFree);
    double fourXSGenieDefaultNuc[4] = {0.0,0.0,0.0,0.0}; dummy = getFourXSForInterpolation(ar1GenieDefaultNuc,ar2GenieDefaultNuc,ar3GenieDefaultNuc,fourXSGenieDefaultNuc);

    double myx = x;
    double myQ2 = Q2;
    if(!absolute){
      // make ratios of the four points we will interpolate
      // order matters, I'm overwriting things, not keeping safe copies.
      for(int i=0; i<4; i++)fourXSNLOTMCHTCG[i] /= fourXSNLOTMCHT[i];
      for(int i=0; i<4; i++)fourXSNLOTMCHTCa40[i] /= fourXSNLOTMCHTC12[i];
      for(int i=0; i<4; i++)fourXSNLOTMCHTC12[i] /= fourXSNLOTMCHT[i];
    
      for(int i=0; i<4; i++)fourXSGenieBugfixFree[i] /= fourXSGenieDefaultFree[i];
      for(int i=0; i<4; i++)fourXSGenieDefaultNuc[i] /= fourXSGenieDefaultFree[i];
      for(int i=0; i<4; i++)fourXSGenieDefaultFree[i] /= fourXSNLOTMCHT[i];    

      for(int i=0; i<4; i++)fourXSNLOTMCHT[i] /= fourXSNLOTMC[i];
      for(int i=0; i<4; i++)fourXSNLOTMC[i] /= fourXSNLO[i];
      for(int i=0; i<4; i++)fourXSNLO[i] /= fourXSLO[i];
      for(int i=0; i<4; i++)fourXSLO[i] /= fourXSLO[i];

      myx = mXforStructureFunction;
      myQ2 = mQ2forStructureFunction;
  
    }

    if(0)std::cout << "fourXSNLOTMCHT " << fourXSNLOTMCHT[0] << " " << fourXSNLOTMCHT[1] << " " << fourXSNLOTMCHT[2] << " " << fourXSNLOTMCHT[3] << std::endl;
      //std::cout << "is absolute " << XSLO << " " << XSNLO << " " << XSNLOTMC << std::endl;

    // not sure.  If I am outside the bounds, do I still want to use the structure function bounds?
    // I think probably yes.   Something we should test.   Not sure how much it will matter.
    myx = mXforStructureFunction;
    myQ2 = mQ2forStructureFunction;

    weights[0] = getInterpolatedValueFromFour(myx,myQ2,fourXSLO);
    weights[1] = getInterpolatedValueFromFour(myx,myQ2,fourXSNLO);
    weights[2] = getInterpolatedValueFromFour(myx,myQ2,fourXSNLOTMC);
    weights[3] = getInterpolatedValueFromFour(myx,myQ2,fourXSNLOTMCHT);
    weights[4] = getInterpolatedValueFromFour(myx,myQ2,fourXSNLOTMCHTCG);
    weights[5] = getInterpolatedValueFromFour(myx,myQ2,fourXSNLOTMCHTC12);
    weights[6] = getInterpolatedValueFromFour(myx,myQ2,fourXSNLOTMCHTCa40);
    weights[7] = getInterpolatedValueFromFour(myx,myQ2,fourXSGenieDefaultFree);
    weights[8] = getInterpolatedValueFromFour(myx,myQ2,fourXSGenieBugfixFree);
    weights[9] = getInterpolatedValueFromFour(myx,myQ2,fourXSGenieDefaultNuc);

    // There used to be code that interpolated the absolute cross section first,
    // then took the ratio of those values.   The function calls themselves still exist.
  }
  */
  
  return weights[1];


}
 

inline double weightAMUvariations::getWeight(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double *weights){
  // this is just a wrapper.   At the moment, I don't need any such thing, consider eliminating it.
  
  return getWeightInternal(x, y, Enu, mlep, isAntinu, weights);

}

inline double weightAMUvariations::getAbsolute(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double *weights){
  // this is just a wrapper.   At the moment, I don't need any such thing, consider eliminating it.
  
  int absolute = 1;
  return getWeightInternal(x, y, Enu, mlep, isAntinu, weights, absolute);

}

inline double weightAMUvariations::getF1(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double *weights){
  int absolute = 11;
  return getWeightInternal(x, y, Enu, mlep, isAntinu, weights, absolute);
}

inline double weightAMUvariations::getF2(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double *weights){
  int absolute = 12;
  return getWeightInternal(x, y, Enu, mlep, isAntinu, weights, absolute);
}

inline double weightAMUvariations::getF3(const double x, const double y, const double Enu, const double mlep, const int isAntinu, double *weights){
  int absolute = 13;
  return getWeightInternal(x, y, Enu, mlep, isAntinu, weights, absolute);
}

inline double weightAMUvariations::getF1withQ2(const double x, const double Q2, const double Enu, const double mlep, const int isAntinu, double *weights){
  int absolute = 11;
  //double MN = 0.5 * (0.93827208943 + 0.93956542052);
  double y = Q2 / (2.0 * x * MN * Enu);
  //double Q2 = 2.0 * x * MN * (y * Enu);
  return getWeightInternal(x, y, Enu, mlep, isAntinu, weights, absolute);
}

inline double weightAMUvariations::getF2withQ2(const double x, const double Q2, const double Enu, const double mlep, const int isAntinu, double *weights){
  int absolute = 12;
  //double MN = 0.5 * (0.93827208943 + 0.93956542052);
  double y = Q2 / (2.0 * x * MN * Enu);
  return getWeightInternal(x, y, Enu, mlep, isAntinu, weights, absolute);
}

inline double weightAMUvariations::getF3withQ2(const double x, const double Q2, const double Enu, const double mlep, const int isAntinu, double *weights){
  int absolute = 13;
  //double MN = 0.5 * (0.93827208943 + 0.93956542052);
  double y = Q2 / (2.0 * x * MN * Enu);
  return getWeightInternal(x, y, Enu, mlep, isAntinu, weights, absolute);
}


inline double weightAMUvariations::getIntegral(const double Enu, const double mlep, const int isAntinu, double *weights, double W){

  double xstep = 0.001;  // guessing this is enough precision
  double ystep = 0.001;
  
  double tempweights[NWEIGHTS];
  double integral[NWEIGHTS];
  for(int i=0; i<NWEIGHTS; i++){tempweights[i] = 0.0; integral[i] = 0.0;};
  
  for(double x = 0.0 + 0.5*xstep; x < 1.0; x += xstep){
    for (double y = 0.0 + 0.5*ystep; y < 1.0; y += ystep){
      double onestep = getAbsolute(x,y,Enu,mlep,isAntinu,tempweights);
      //std::cout << "onestep Exy " << " " << Enu << " " << x << " " << y << " " 
      //          << tempweights[0] << " " << tempweights[1] << std::endl;

      // To do, add an optional W cut here.
      bool cutW = false;
      if(W > 0){
        double thisQ2 = 2.0 * MN * x * y * Enu;
        double W2 = MN * MN + 2.0 * MN * y * Enu - thisQ2;
        if(sqrt(W2) < W)cutW = true;
      }

      if(!cutW)for(int i=0; i<NWEIGHTS; i++){integral[i] += tempweights[i];}
    }
  }
  // multiply by the step size
  for(int i=0; i<NWEIGHTS; i++){integral[i] *= xstep * ystep;}

  for(int i=0; i<NWEIGHTS; i++){weights[i] = integral[i];}

  return weights[0];

}


inline double weightAMUvariations::getOneCrossSection(double arrayf1[][Q2BINS], double arrayf2[][Q2BINS], double arrayf3[][Q2BINS], double x){

  double F1 = interpolate2D(this->mXforStructureFunction, this->mQ2forStructureFunction, arrayf1);
  double F2 = interpolate2D(this->mXforStructureFunction, this->mQ2forStructureFunction, arrayf2);
  double F3 = interpolate2D(this->mXforStructureFunction, this->mQ2forStructureFunction, arrayf3);

  return getOneCrossSection(F1,F2,F3,x);

}

inline double weightAMUvariations::getOneCrossSection(double F1, double F2, double F3, double x){

  //double weightAMUvariations::getOneCrossSection(double arrayf1[][Q2BINS], double arrayf2[][Q2BINS], double arrayf3[][Q2BINS], double x){

  //  double F1 = interpolate2D(this->mXforStructureFunction, this->mQ2forStructureFunction, arrayf1);
  //  double F2 = interpolate2D(this->mXforStructureFunction, this->mQ2forStructureFunction, arrayf2);
  //  double F3 = interpolate2D(this->mXforStructureFunction, this->mQ2forStructureFunction, arrayf3);

  // The interpolation can return a negative number beyond the edges
  // Catch it here just in case the interpolate2D doesn't do it.
  // an example is x = 0.95 (restricted to 0.92) and y = 0.99 and Enu = 50.0
  // Best to do a better job of restricting the mQ2forStructureFunction range.
  if(F1 < 0.0)F1 = 0.0;
  if(F2 < 0.0)F2 = 0.0;
  if(F3 < 0.0)F3 = 0.0;

  double answer = 0.0;

  
  if(x >= 0.0 && x <= 1.0){
    // use Callan Gross F1 = F2 / 2x if user puts in a valid value for 0 < x < 1 
    // the input x defaults to -1.0 in the definition
    
    // this first one does not include the TMC-looking term Huma calls gamma
    // equivalently Huma says gamma = 1.0 
    //answer = mFront * (mF1term*F2*0.5/x + mF2term*F2 + mSign * mF3term*F3);
    // This one with the TMC-looking term but still no FL
    answer = mFront * (mF1termCG*F2*0.5/x + mF2term*F2 + mSign * mF3term*F3 + mF5term*F2*0.5/x); 

  } else {
    //answer = mFront * (mF1term*F1 + mF2term*F2 + mSign * mF3term*F3);
    answer = mFront * (mF1term*F1 + mF2term*F2 + mSign * mF3term*F3 + mF5term*F2*0.5/x); 
    // set F4 to zero for now, mF5term includes 1/2x here.
  }

  // if absolute == 1, no change
  // these bypass the calculation and just return the structure function
  // so the user can validate the interpolation 
  // or make direct comparisons like the theory papers do.
  //if(absolute == 2)answer = F1;
  //else if(absolute == 3)answer = F2;
  //else if(absolute == 4)answer = F3;  // traditionally xF3 .

  if(answer < 0)std::cout << "negative " << answer << " F1 " << mF1term << " " << F1 << " F2 " << mF2term << " " << F2 
    << " F3 " << mSign << " " << mF3term << " " << F3 << std::endl;

  // To be fixed or at least to be understood
  // either a bug or could be natural, the cross section goes zero and negative for a pair of these inputs
  // the base model (no NLO, no nothing) and the GENIE bugfix.  And only for 

  // use this to trap values that are less than zero.
  // don't set it to zero exactly, because that will give divide by zeros.
  // later on code probably also wants to trap obscenely large weights ?
  double threshold = 1.0e-43;
  if(answer < threshold)answer = 1.0e-43;

  return answer;

}

//
//  Code up my own interpolation so I don't have to use the TH2D for these
//  This version just needs a rectangle grid of points in an array
//  And doesn't have any of the rest of the TH2D overhead
//

inline double weightAMUvariations::getFourXSForInterpolation(const double array1[][Q2BINS], const double array2[][Q2BINS], const double array3[][Q2BINS], double fourXS[4], double x){

  // turn F1 F2 F3 arrays into the four points to interpolate, then get the cross section for each interpolation, then interpolate.
  // these are going to use a common set of mx1, mx2, mqq1, mqq2 set the first time getFourValues is called.

  double F1[4] = {0.0,0.0,0.0,0.0};
  double F2[4] = {0.0,0.0,0.0,0.0};
  double F3[4] = {0.0,0.0,0.0,0.0};

  // the return values are bogus, the filled array is the thing.
  double retFourF1 = getFourValuesForInterpolation(mXforStructureFunction,mQ2forStructureFunction,array1,F1);//F1[0],F1[1],F1[2],F1[3]);
  double retFourF2 = getFourValuesForInterpolation(mXforStructureFunction,mQ2forStructureFunction,array2,F2);//F1[0],F2[1],F2[2],F2[3]);
  double retFourF3 = getFourValuesForInterpolation(mXforStructureFunction,mQ2forStructureFunction,array3,F3);//F1[0],F3[1],F3[2],F3[3]);

  if(0)std::cout << "getFourXSForInterpolationFourF1 " << F1[0] << " " << F1[1] << " " << F1[2] << " " << F1[3] << " xQ2 " << mXforStructureFunction << " " << mQ2forStructureFunction << std::endl;

  if(mAbsolute <= 1){
    //double fourXS[4] = {0.0,0.0,0.0,0.0};
    for(int i=0; i<4; i++){fourXS[i] = getOneCrossSection(F1[i],F2[i],F3[i],x);}
    //fourXS[0] = getOneCrossSection(F1[0],F2[0],F3[0],x);
    //fourXS[1] = getOneCrossSection(F1[1],F2[1],F3[1],x);
    //fourXS[2] = getOneCrossSection(F1[2],F2[2],F3[2],x);
    //fourXS[3] = getOneCrossSection(F1[3],F2[3],F3[3],x);
  } else {
    // user requested just one of the structure functions.
    // will only return absolutes, user will need to take ratios themselves.
    if(mAbsolute == 11){
      for(int i=0; i<4; i++){fourXS[i] = F1[i];}
    } else if(mAbsolute == 12){
      for(int i=0; i<4; i++){fourXS[i] = F2[i];}
    } else if(mAbsolute == 13){
      for(int i=0; i<4; i++){fourXS[i] = F3[i];}
    }
  }

  if(0)std::cout << "getFourXSForInterpolationFourXS " << fourXS[0] << " " << fourXS[1] << " " << fourXS[2] << " " << fourXS[3] << std::endl;

  return fourXS[0];   // return the first thing, but it doesn't matter really.

}


//double weightAMUvariations::getFourValuesForInterpolation(double x, double y, const double array[][Q2BINS], double q11, double q12, double q21, double q22){
inline double weightAMUvariations::getFourValuesForInterpolation(double x, double y, const double array[][Q2BINS], double fourval[4]){

  // This bilinear interpolation function is lightweight enough
  // its not important to me if I'm reimplementing something
  // it saves a library dependency too, though gsl has bilinear and bicubic
  //
  // this one explicitly puts in the array binning
  // pass in x = xbj and y = Q2
  // also pass in the array and the number of elements in the array
  // or simply use a c++ vector which knows both at once.


  // first step, find the point xy in the array
  // that is equal to or higher than the input xy
  // then step back one point in both x and y
  // and that defines the four points to use in interpolation.
  // special treatment of bound here.


  // test if this is already done
  if(mx1 < 0.0 || mx2 < 0.0 || mqq1 < 0.0 || mqq2 < 0.0){

    int xmaxindex = XBINS - 1;  // array index max = size - 1
    // initialize these to the highest or one before
    // if the input x is higher than the highest, this will be the result
    int indexx = xmaxindex - 1;
    double x2 = xArray[xmaxindex]; double x1 = xArray[xmaxindex-1];
    //x2 = xArray[xmaxindex]; x1 = xArray[xmaxindex-1];

    for(int i=0; i<XBINS; i++){
      //std::cout << "i " << i << " " << x << " " << xArray[i] << " " << std::endl;
       if(x <= xArray[i]){
        indexx = i - 1;
        if(indexx < 0)indexx = 0;
        x1 = xArray[indexx];
        x2 = xArray[indexx +1];
       break;
      }
   }
    int ymaxindex = Q2BINS - 1;
    int indexy = ymaxindex - 1;
    double y2 = q2Array[ymaxindex]; double y1 = q2Array[ymaxindex-1];
    //y2 = q2Array[ymaxindex]; y1 = q2Array[ymaxindex-1];

    // if the value presented is smaller than the first Q2 point
    // the interpolation routine will extrapolate below the array.
    // I guess that happens above too.
    for(int i=0; i < Q2BINS; i++){
      if(y <= q2Array[i]){
        indexy = i - 1;
        if(indexy < 0)indexy = 0;
        y1 = q2Array[indexy];
        y2 = q2Array[indexy + 1];
        break;
      }
    }
    mxarrayindex = indexx;
    mqqarrayindex = indexy;
    mx1 = x1;
    mx2 = x2;
    mqq1 = y1;
    mqq2 = y2;

    if(0)std::cout << "setting xQ2 " << mx1 << " " << mx2 << " " << mqq1 << " " << mqq2 << " index " << mxarrayindex << " " << mqqarrayindex << " fromxQ2 " << x << " " << y << std::endl;

  }

  // now with the lower bound in x and y, I've got four points
  // q11 is that first point
  // q12, q21, and q22 are the others.
  // this will extrapolate too, though probably want to set bounds to prevent it.

  //q11 = array[indexx][indexy];
  //q12 = array[indexx][indexy+1];
  //q21 = array[indexx+1][indexy];
  //q22 = array[indexx+1][indexy+1];

  fourval[0] = array[mxarrayindex][mqqarrayindex]; //q11
  fourval[1] = array[mxarrayindex][mqqarrayindex+1]; //q12
  fourval[2] = array[mxarrayindex+1][mqqarrayindex];  //q21
  fourval[3] = array[mxarrayindex+1][mqqarrayindex+1];  //q22

  if(0){
    double q11 = fourval[0];  double q12 = fourval[1]; double q21 = fourval[2]; double q22 = fourval[3];
    std::cout << "fourpoints " << q11 << " " << q12 << " " << q21 << " " << q22 << " index " << mxarrayindex << " " << mqqarrayindex << " fromxQ2 " << x << " " << y << std::endl;
  }

/*

  // this is the formula in wikipedia and Andy Mastbaum's code
  double denom = (x2-x1)*(y2-y1);  // the common denominator
  double interpvalue =  q11*(x2-x)*(y2-y) +  q21*(x-x1)*(y2-y)
                     +  q12*(x2-x)*(y-y1) +  q22*(x-x1)*(y-y1);
  interpvalue = interpvalue / denom;

  if(0){
    std::cout << "indexes are x " << indexx << " " << x1 << " " << x2 << " y " << indexy << " " << y1 << " " << y2 
              << " for xQ2 " << x << " " <<  y << std::endl; 
    std::cout << "qs are " << q11 << " " << q12 << " " << q21 << " " << q22 << " interp " << interpvalue << std::endl;
  }
  return interpvalue;

  */
  return fourval[0];

}

inline double weightAMUvariations::getInterpolatedValueFromFour(double x, double y, double array[4]){
  // hmm.  x and y (really x and Q2) mean different things depending what order I do this
  // because if I am only interpolating structure functions, there is bounds checking on whether I'm inside the input arrays.
  return getInterpolatedValueFromFour(x,y,array[0],array[1],array[2],array[3]);
}

inline double weightAMUvariations::getInterpolatedValueFromFour(double x, double y, double q11, double q12, double q21, double q22){
  // Trusting the user has found the best four values for the interpolation
  // and has submitted them in order.   No sanity checking done here.

  // Still using mx1 mx2 my1 my2 saved as class members
  // because they are presumably the same, if the F1 F2 F3 arrays are the same.
  if(mx1 < 0.0 || mx2 < 0.0 || mqq1 < 0.0 || mqq2 < 0.0){

    return -1.0;   // actually should throw an error.
  }

  // this is the formula in wikipedia and Andy Mastbaum's code
  double denom = (mx2-mx1)*(mqq2-mqq1);  // the common denominator
  double interpvalue =  q11*(mx2-x)*(mqq2-y) +  q21*(x-mx1)*(mqq2-y)
                     +  q12*(mx2-x)*(y-mqq1) +  q22*(x-mx1)*(y-mqq1);
  interpvalue = interpvalue / denom;

  if(0){
    std::cout << "indexes are x " << " " << mx1 << " " << mx2 << " y " << " " << mqq1 << " " << mqq2 
              << " for xQ2 " << x << " " <<  y << std::endl; 
    std::cout << "qs are " << q11 << " " << q12 << " " << q21 << " " << q22 << " interp " << interpvalue << std::endl;  
  }
  return interpvalue;

}

//double weightAMUvariations::interpolate2D(double x, double y, const double array[XBINS][Q2BINS], const std::vector<double> xbins, const std::vector<double> ybins){
  inline double weightAMUvariations::interpolate2D(double x, double y, const double array[][Q2BINS]){

    //  This function is still here and an active member of the class, but
    //  I've separated the parts because numerical precision is probably better
    //  if we get the q11 q12 q21 q22 , then take the ratio of them, then interpolate
    //  Not interpolate first then take the ratio.

    // This bilinear interpolation function is lightweight enough
    // its not important to me if I'm reimplementing something
    // it saves a library dependency too, though gsl has bilinear and bicubic
    //
    // this one explicitly puts in the array binning
    // pass in x = xbj and y = Q2
    // also pass in the array and the number of elements in the array
    // or simply use a c++ vector which knows both at once.
  
  
    // first step, find the point xy in the array
    // that is equal to or higher than the input xy
    // then step back one point in both x and y
    // and that defines the four points to use in interpolation.
    // special treatment of bound here.
    int xmaxindex = XBINS - 1;  // array index max = size - 1
    // initialize these to the highest or one before
    // if the input x is higher than the highest, this will be the result
    int indexx = xmaxindex - 1;
    double x2 = xArray[xmaxindex]; double x1 = xArray[xmaxindex-1];
  
    for(int i=0; i<XBINS; i++){
      //std::cout << "i " << i << " " << x << " " << xArray[i] << " " << std::endl;
       if(x <= xArray[i]){
        indexx = i - 1;
        if(indexx < 0)indexx = 0;
        x1 = xArray[indexx];
        x2 = xArray[indexx +1];
        break;
       }
    }
    int ymaxindex = Q2BINS - 1;
    int indexy = ymaxindex - 1;
    double y2 = q2Array[ymaxindex]; double y1 = q2Array[ymaxindex-1];
    // if the value presented is smaller than the first Q2 point
    // the interpolation routine will extrapolate below the array.
    // I guess that happens above too.
    for(int i=0; i < Q2BINS; i++){
      if(y <= q2Array[i]){
        indexy = i - 1;
        if(indexy < 0)indexy = 0;
        y1 = q2Array[indexy];
        y2 = q2Array[indexy + 1];
        break;
      }
    }
  
    // now with the lower bound in x and y, I've got four points
    // q11 is that first point
    // q12, q21, and q22 are the others.
    // this will extrapolate too, though probably want to set bounds to prevent it.
  
    double q11 = array[indexx][indexy];
    double q12 = array[indexx][indexy+1];
    double q21 = array[indexx+1][indexy];
    double q22 = array[indexx+1][indexy+1];
  
  
  
    // this is the formula in wikipedia and Andy Mastbaum's code
    double denom = (x2-x1)*(y2-y1);  // the common denominator
    double interpvalue =  q11*(x2-x)*(y2-y) +  q21*(x-x1)*(y2-y)
                       +  q12*(x2-x)*(y-y1) +  q22*(x-x1)*(y-y1);
    interpvalue = interpvalue / denom;
  
    if(0){
      std::cout << "indexes are x " << indexx << " " << x1 << " " << x2 << " y " << indexy << " " << y1 << " " << y2 
                << " for xQ2 " << x << " " <<  y << std::endl; 
      std::cout << "qs are " << q11 << " " << q12 << " " << q21 << " " << q22 << " interp " << interpvalue << std::endl;
    }
    return interpvalue;
  
  }

}
