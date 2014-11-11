#include "CEffZFitter.hh"
//#include "CEffPlotter.hh"
#include "KStyle.hh"

#include <iostream>
#include <string>
#include <cstdlib>

int main(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings and constants
  //==============================================================================================================

  // handle input arguments
  const std::string conf     = argv[1];         // input configuration file
  const int         sigpass  = atoi(argv[2]);	// signal model for PASS sample
  const int         bkgpass  = atoi(argv[3]);	// background model for PASS sample
  const int         sigfail  = atoi(argv[4]);	// signal model for FAIL sample
  const int         bkgfail  = atoi(argv[5]);	// background model for FAIL sample
  const std::string infname  = argv[6];         // input ROOT file of probes
  const std::string outdir   = argv[7];         // output directory
  const int         doPU     = atoi(argv[8]);	// PU re-weighting mode
  const int         charge   = atoi(argv[9]);	// probe charge requirement (0, -1, +1)
  const std::string temfname = argv[10];        // ROOT file for generating MC-based templates
  
  // other settings
  const double       massLo    = 60;
  const double       massHi    = 120;
  const double       fitMassLo = massLo;
  const double       fitMassHi = massHi;
  const int          uncMethod = 0;
  const std::string  pufname   = doPU ? "PUWeights_2012.root" : "none";
  const unsigned int runNumLo  = 0;
  const unsigned int runNumHi  = 999999;
  
  std::cout << std::endl;
  std::cout << " <> Using binning configuration: " << conf << std::endl;
  std::cout << " <> Processing probes file: " << infname << std::endl;
  std::cout << std::endl;

  KStyle();  
  
  //--------------------------------------------------------------------------------------------------------------
  // Measure efficiency using Z events
  //==============================================================================================================
/*
  CEffZFitter
  initialize()
    ==> parseConf()
    ==> set up output dirs, PU weights, pass/fail trees, templates
    ==> makeBinnedTemplates() / makeUnbinnedTemplates()
  computeEff()
    ==> makeEffGraph()
      ==> loop over bins
      ==> performCount()
      ==> parseFitResults()
      ==> performFit()
        ==> printCorrelations
    ==> makeEffHist2D()
      ==> loop over bins
      ==> performCount()
      ==> parseFitResults()
      ==> performFit()
        ==> printCorrelations
    ==> produce output
*/
  CEffZFitter fitter;
  fitter.initialize(conf, sigpass, bkgpass, sigfail, bkgfail,
                    infname, outdir, temfname,
                    massLo, massHi, fitMassLo, fitMassHi, 
		    uncMethod, pufname, charge, runNumLo, runNumHi);
  fitter.computeEff();
  

  //--------------------------------------------------------------------------------------------------------------
  // Produce summary plots
  //==============================================================================================================
  
//  CEffPlotter plotter;
//  plotter.drawPlots(outdir);
  
  std::cout << " <> Output saved in " << outdir << std::endl;

  return 0;
}
