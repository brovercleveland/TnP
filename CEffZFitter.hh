#ifndef CEFFZFITTER_HH
#define CEFFZFITTER_HH

//================================================================================================
//
// Signal Extraction
//-------------------
//  0: probe counting
//  1: Breit-Wigner convolved with Crystal Ball function
//  2: MC template convolved with Gaussian
//  3: Phil's Crystal Ball based "Voigtian" shape
//  4: Unbinned MC data convolved with Gaussian
//
// Background Model
//------------------
//  0: no background
//  1: exponential model
//  2: erfc*exp model
//  3: double exponential model
//  4: linear*exp model
//  5: quadratic*exp model
//
//________________________________________________________________________________________________

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

class TTree;
class TCanvas;
class TGraphAsymmErrors;
class TH1D;
class TH2D;
class RooFitResult;

class CEffZFitter
{
public:
  CEffZFitter();
  ~CEffZFitter();
  
  void initialize(const std::string conf, const int sigpass, const int bkgpass, const int sigfail, const int bkgfail,
                  const std::string infname, const std::string outdir, const std::string temfname,
                  const double massLo, const double massHi, const double fitMassLo, const double fitMassHi, 
		  const int uncMethod, const std::string pufname, const int charge,
		  const unsigned int runNumLo, const unsigned int runNumHi); 
  
  void computeEff();
  
    
protected:

  // parse binning configuration file
  void parseConf(const std::string conf);
  
  // generate MC-based signal templates
  void makeBinnedTemplates(const std::string temfname, const int charge, TH1D *puWeights);
  void makeUnbinnedTemplates(const std::string temfname, const int charge);
  
  // make efficiency graph
  TGraphAsymmErrors* makeEffGraph(const std::vector<double> &edgesv, 
                                  const std::vector<TTree*> &passv,
				  const std::vector<TTree*> &failv,
                                  const std::string name);
  
  // make 2D efficiency map
  void makeEffHist2D(TH2D *hEff, TH2D *hErrl, TH2D *hErrh, 
                     const std::vector<TTree*> &passv,
		     const std::vector<TTree*> &failv,
                     const std::string name);
  
  // efficiency computation
  void performCount(double &resEff, double &resErrl, double &resErrh,
                    const int ibin,
		    const double xbinLo, const double xbinHi,
		    const double ybinLo, const double ybinHi,
                    TTree *passTree, TTree *failTree, 
                    const std::string name, TCanvas *cpass, TCanvas *cfail);
  void performFit(double &resEff, double &resErrl, double &resErrh,
                  const int ibin,
		  const double xbinLo, const double xbinHi,
		  const double ybinLo, const double ybinHi,
                  TTree *passTree, TTree *failTree,
                  const std::string name, TCanvas *cpass, TCanvas *cfail);
  
  // print correlations of fit parameters
  void printCorrelations(std::ostream& os, RooFitResult* res);
  
  // parse fit results file
  void parseFitResults(std::ifstream &ifs, double &eff, double &errl, double &errh);

  // create HTML page
  void makeHTML(const std::string name, const unsigned int nbins);
  
  ///// data members /////
  
  bool fIsInitialized;
  
  // signal and background models
  int fSigPass, fBkgPass, fSigFail, fBkgFail;
  
  double fMassLo, fMassHi;        // signal extraction mass window  
  double fFitMassLo, fFitMassHi;  // fit mass window
  
  // efficiency uncertainty calculation method
  // method: 0 -> Clopper-Pearson
  //         1 -> Feldman-Cousins
  int fUncMethod;
  
  // output directory for results
  std::string fOutputDir;
      
  // bin edges for kinematic
  std::vector<double> fPtBinEdgesv;
  std::vector<double> fEtaBinEdgesv;
  std::vector<double> fPhiBinEdgesv;
  std::vector<double> fNPVBinEdgesv;
  
  // flags for |eta| and |phi| binning
  bool fDoAbsEta, fDoAbsPhi;
  
  // flags for binnings to compute efficiencies for
  bool fDoPt, fDoEta, fDoPhi, fDoEtaPt, fDoEtaPhi, fDoNPV;
  
  // trees for pass/fail samples
  std::vector<TTree*> fPassTreePtv,     fFailTreePtv;
  std::vector<TTree*> fPassTreeEtav,    fFailTreeEtav;
  std::vector<TTree*> fPassTreePhiv,    fFailTreePhiv;
  std::vector<TTree*> fPassTreeEtaPtv,  fFailTreeEtaPtv;
  std::vector<TTree*> fPassTreeEtaPhiv, fFailTreeEtaPhiv;
  std::vector<TTree*> fPassTreeNPVv,    fFailTreeNPVv;
};

#endif
