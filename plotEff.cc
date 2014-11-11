#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH2D.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>      // graph class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O

#include "CPlot.hh"	            // helper class for plots
#include "KStyle.hh"                // style settings for drawing

using namespace std;

void main()
{ 
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // handle input arguments
  const string outdir   = argv[1];
  const string mcfname  = argv[2];
  const string datfname = argv[3];

  vector<string> etalabelv;
  vector<string> ptlabelv;
  
  // eta bins
  etalabelv.push_back("0 < |#eta| < 1.5");
  etalabelv.push_back("1.5 < |#eta| < 2.5");
  
  CPlot::sOutDir = outdir.c_str();
  const string format = "png";

  // y-axis ranges
  double efflow   = 0.10, effhigh   = 1.20;
  double scalelow = 0.95, scalehigh = 1.05;
  
  KStyle();  


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TFile mcfile(mcfname);      
  TFile datafile(datfname);
  
  TH2D *hMCEff=0,   *hMCErrl=0,   *hMCErrh=0;
  TH2D *hDataEff=0, *hDataErrl=0, *hDataErrh=0;
  
  TGraphAsymmErrors *grMCEffEta=0, *grDataEffEta=0, *grScaleEta=0;
  TGraphAsymmErrors *grMCEffPt=0,  *grDataEffPt=0,  *grScalePt=0;
  
  vector<TGraphAsymmErrors*> mceff_vs_pt_per_etav;
  vector<TGraphAsymmErrors*> eff_vs_pt_per_etav;
  vector<TGraphAsymmErrors*> mceff_vs_eta_per_ptv;
  vector<TGraphAsymmErrors*> eff_vs_eta_per_ptv;
  vector<TGraphAsymmErrors*> scale_vs_pt_per_etav;
  vector<TGraphAsymmErrors*> scale_vs_eta_per_ptv;
  
  grMCEffEta = (TGraphAsymmErrors*)mcfile.Get("grEffEta");
  grMCEffPt  = (TGraphAsymmErrors*)mcfile.Get("grEffPt");
  if(grMCEffPt) {
    grMCEffPt->SetPoint(grMCEffPt->GetN()-1, 125, grMCEffPt->GetY()[grMCEffPt->GetN()-1]);
    grMCEffPt->SetPointError(grMCEffPt->GetN()-1,25,25, 
                             grMCEffPt->GetErrorYlow(grMCEffPt->GetN()-1),
			     grMCEffPt->GetErrorYhigh(grMCEffPt->GetN()-1));
  }  			               
  hMCEff  = (TH2D*)mcfile.Get("hEffEtaPt");
  hMCErrl = (TH2D*)mcfile.Get("hErrlEtaPt");
  hMCErrh = (TH2D*)mcfile.Get("hErrhEtaPt");  

  grDataEffEta = (TGraphAsymmErrors*)datafile.Get("grEffEta");
  grDataEffPt  = (TGraphAsymmErrors*)datafile.Get("grEffPt");
  if(grDataEffPt) {
    grDataEffPt->SetPoint(grDataEffPt->GetN()-1, 125, grDataEffPt->GetY()[grDataEffPt->GetN()-1]);
    grDataEffPt->SetPointError(grDataEffPt->GetN()-1,25,25, 
                                grDataEffPt->GetErrorYlow(grDataEffPt->GetN()-1),
			        grDataEffPt->GetErrorYhigh(grDataEffPt->GetN()-1));
  }
  
  hDataEff  = (TH2D*)datafile.Get("hEffEtaPt");
  hDataErrl = (TH2D*)datafile.Get("hErrlEtaPt");
  hDataErrh = (TH2D*)datafile.Get("hErrhEtaPt");
    
  if(grMCEffEta && grDataEffEta) {
    grScaleEta = new TGraphAsymmErrors(grMCEffEta->GetN());
    for(int i=0; i<grMCEffEta->GetN(); i++) {
      double mcval   = grMCEffEta->GetY()[i];
      double dataval = grDataEffEta->GetY()[i];
      double scale   = dataval/mcval;
      grScaleEta->SetPoint(i,grMCEffEta->GetX()[i],scale);
      
      double mcerrl   = grMCEffEta->GetErrorYlow(i);
      double mcerrh   = grMCEffEta->GetErrorYhigh(i);
      double dataerrl = grDataEffEta->GetErrorYlow(i);
      double dataerrh = grDataEffEta->GetErrorYhigh(i);
      grScaleEta->SetPointError(i, 0, 0,
	  		        scale*sqrt(mcerrl*mcerrl/mcval/mcval + dataerrl*dataerrl/dataval/dataval),
			        scale*sqrt(mcerrh*mcerrh/mcval/mcval + dataerrh*dataerrh/dataval/dataval));
    }
  }
  
  if(grMCEffPt && grDataEffPt) {
    grScalePt = new TGraphAsymmErrors(grMCEffPt->GetN());
    for(int i=0; i<grMCEffPt->GetN(); i++) {
      Double_t mcval   = grMCEffPt->GetY()[i];
      Double_t dataval = grDataEffPt->GetY()[i];
      Double_t scale   = dataval/mcval;
      grScalePt->SetPoint(i, grMCEffPt->GetX()[i],scale);
      if(i==grMCEffPt->GetN()-1)
        grScalePt->SetPoint(i,165,scale);
    
      Double_t mcerrl   = grMCEffPt->GetErrorYlow(i);
      Double_t mcerrh   = grMCEffPt->GetErrorYhigh(i);
      Double_t dataerrl = grDataEffPt->GetErrorYlow(i);
      Double_t dataerrh = grDataEffPt->GetErrorYhigh(i);
      grScalePt->SetPointError(i, 0, 0,
			        scale*sqrt(mcerrl*mcerrl/mcval/mcval + dataerrl*dataerrl/dataval/dataval),
			        scale*sqrt(mcerrh*mcerrh/mcval/mcval + dataerrh*dataerrh/dataval/dataval));
    }
  }
  
  if(hMCEff->GetEntries()>0 && hDataEff->GetEntries()>0) {
    const int nx = hMCEff->GetNbinsX();
    const int ny = hMCEff->GetNbinsY(); 
/*    
    for(int iy=1; iy<=ny; iy++) {
      double xval[nx], xerr[nx];
      double effval[nx],   efferrl[nx],   efferrh[nx];
      double scaleval[nx], scaleerrl[nx], scaleerrh[nx];
      
      for(int ix=1; ix<=nx; ix++) {
        xval[ix-1] = 0.5*(hMCEff->GetXaxis()->GetBinLowEdge(ix) + hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
        xerr[ix-1] = 0;
        
        double mceff  = hMCEff->GetCellContent(ix,iy);
        double mcerrl = hMCErrl->GetCellContent(ix,iy);
        double mcerrh = hMCErrh->GetCellContent(ix,iy);
        
        double dataeff  = hDataEff->GetCellContent(ix,iy);
        double dataerrl = hDataErrl->GetCellContent(ix,iy);
        double dataerrh = hDataErrh->GetCellContent(ix,iy);
        effval[ix-1]  = dataeff;
        efferrl[ix-1] = dataerrl;
        efferrh[ix-1] = dataerrh;
        
        double scale = dataeff/mceff;
        scaleval[ix-1]  = scale;
        scaleerrl[ix-1] = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
        scaleerrh[ix-1] = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      }
      
      eff_vs_eta_per_ptv.push_back(new TGraphAsymmErrors(nx,xval,effval,xerr,xerr,efferrl,efferrh));
      scale_vs_eta_per_ptv.push_back(new TGraphAsymmErrors(nx,xval,scaleval,xerr,xerr,scaleerrl,scaleerrh));
    }
*/    
    for(int ix=1; ix<=nx; ix++) {
      double xval[ny], xerr[ny];
      double mceffval[ny], mcefferrl[ny], mcefferrh[ny];
      double effval[ny],   efferrl[ny],   efferrh[ny];
      double scaleval[ny], scaleerrl[ny], scaleerrh[ny];
      
      for(int iy=1; iy<=ny; iy++) {
        xval[iy-1] = 0.5*(hMCEff->GetYaxis()->GetBinLowEdge(iy+1) + hMCEff->GetYaxis()->GetBinLowEdge(iy));
        xerr[iy-1] = 0.5*(hMCEff->GetYaxis()->GetBinLowEdge(iy+1) - hMCEff->GetYaxis()->GetBinLowEdge(iy));
	if(iy==ny) {
	  xval[iy-1] = 125;
          xerr[iy-1] = 25;
        }
	
        double mceff  = hMCEff->GetCellContent(ix,iy);
        double mcerrl = hMCErrl->GetCellContent(ix,iy);
        double mcerrh = hMCErrh->GetCellContent(ix,iy);
	mceffval[iy-1]  = mceff;
        mcefferrl[iy-1] = mcerrl;
        mcefferrh[iy-1] = mcerrh;
        
        double dataeff  = hDataEff->GetCellContent(ix,iy);
        double dataerrl = hDataErrl->GetCellContent(ix,iy);
        double dataerrh = hDataErrh->GetCellContent(ix,iy);
        effval[iy-1]  = dataeff;
        efferrl[iy-1] = dataerrl;
        efferrh[iy-1] = dataerrh;
        
        double scale = dataeff/mceff;
        scaleval[iy-1]  = scale;
        scaleerrl[iy-1] = (scale>0) ? scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff) : 0;
        scaleerrh[iy-1] = (scale>0) ? scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff) : 0;
      }
      
      mceff_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,mceffval,xerr,xerr,mcefferrl,mcefferrh));
      eff_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,effval,xerr,xerr,efferrl,efferrh));
      scale_vs_pt_per_etav.push_back(new TGraphAsymmErrors(ny,xval,scaleval,0,0,scaleerrl,scaleerrh));
    }  
  }

    
  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *c = MakeCanvas("c","c",800,600);
  
  if(grMCEffEta && grDataEffEta) {
    CPlot plotEffEta("effeta","","|#eta|","#varepsilon");
    plotEffEta.AddGraph(grMCEffEta,   "MC","",  kRed, kOpenSquare);
    plotEffEta.AddGraph(grDataEffEta,"data","",kBlue,kFullDotLarge);
    plotEffEta.SetYRange(efflow,effhigh);
    plotEffEta.TransLegend(0,-0.5);
    plotEffEta.Draw(c,kTRUE,format);
    
    CPlot plotScaleEta("scaleeta","","|#eta|","scale factor");
    plotScaleEta.AddGraph(grScaleEta,"",kBlue,kFullDotLarge);
    plotScaleEta.AddLine(0,1.0,2.7,1.0,kBlack,7);
    plotScaleEta.SetYRange(scalelow,scalehigh);
    plotScaleEta.SetXRange(0,2.7);
    plotScaleEta.Draw(c,kTRUE,format);
  }
  
  if(grMCEffPt && grDataEffPt) {
    CPlot plotEffPt("effpt","","p_{T} [GeV/c]","#varepsilon");
    plotEffPt.AddGraph(grMCEffPt,   "MC","",  kRed, kOpenSquare);
    plotEffPt.AddGraph(grDataEffPt,"data","",kBlue,kFullDotLarge);
    plotEffPt.SetYRange(efflow,effhigh);
    plotEffPt.SetXRange(0,150);
    plotEffPt.TransLegend(0,-0.5);
    plotEffPt.Draw(c,kTRUE,format);
  
    CPlot plotScalePt("scalept","","p_{T} [GeV/c]","scale factor");
    plotScalePt.AddGraph(grScalePt,"",kBlue,kFullDotLarge);
    plotScalePt.AddLine(0,1.0,150,1.0,kBlack,7);
    plotScalePt.SetXRange(0,150);
    plotScalePt.SetYRange(scalelow,scalehigh);
    plotScalePt.Draw(c,kTRUE,format);
  }
      
  if(mceff_vs_pt_per_etav.size()>0 && eff_vs_pt_per_etav.size()>0) {
    for(unsigned int ig=0; ig<mceff_vs_pt_per_etav.size(); ig++) {
      char pname[100];
      sprintf(pname,"effpt_eta%i",ig);
      CPlot plotEffPt_perEta(pname,"","p_{T} [GeV/c]","#varepsilon");
      plotEffPt_perEta.AddGraph(mceff_vs_pt_per_etav[ig],"MC","",   kRed,  kOpenSquare);
      plotEffPt_perEta.AddGraph(eff_vs_pt_per_etav[ig], "data","", kBlue, kFullDotLarge);
      plotEffPt_perEta.AddTextBox(etalabelv[ig],0.2,0.82,0.4,0.88,0);
      plotEffPt_perEta.TransLegend(0,-0.5); 
      plotEffPt_perEta.SetXRange(0,150);
      plotEffPt_perEta.SetYRange(efflow,effhigh);
      plotEffPt_perEta.Draw(c,kTRUE,format);
    }    
    
    for(unsigned int ig=0; ig<scale_vs_pt_per_etav.size(); ig++) {
      char pname[100];
      sprintf(pname,"scalept_eta%i",ig);
      CPlot plotScalePt_perEta(pname,"","p_{T} [GeV/c]","scale factor");
      plotScalePt_perEta.AddGraph(scale_vs_pt_per_etav[ig],"",kBlue,kFullDotLarge);
      plotScalePt_perEta.AddTextBox(etalabelv[ig],0.7,0.82,0.9,0.88,0);
      plotScalePt_perEta.AddLine(0,1,150,1,kBlack,7);
      plotScalePt_perEta.SetXRange(0,150);
      plotScalePt_perEta.SetYRange(scalelow,scalehigh);
      plotScalePt_perEta.Draw(c,kTRUE,format);
    }
  }

  if(mceff_vs_eta_per_ptv.size()>0 && eff_vs_eta_per_ptv.size()>0) {
    for(unsigned int ig=0; ig<mceff_vs_eta_per_ptv.size(); ig++) {
      char pname[100];
      sprintf(pname,"effeta_pt%i",ig);
      CPlot plotEffEta_perPt(pname,"","|#eta|","#varepsilon");
      plotEffEta_perPt.AddGraph(mceff_vs_eta_per_ptv[ig],"MC","",  kRed,  kOpenSquare);
      plotEffEta_perPt.AddGraph(eff_vs_eta_per_ptv[ig], "data","", kBlue, kFullDotLarge);
      plotEffEta_perPt.AddTextBox(ptlabelv[ig],0.2,0.82,0.4,0.88,0);
      plotEffEta_perPt.SetYRange(efflow,effhigh);
      plotEffEta_perPt.Draw(c,kTRUE,format);
    }
    
    for(unsigned int ig=0; ig<scale_vs_eta_per_ptv.size(); ig++) {
      char pname[100];
      sprintf(pname,"scaleeta_pt%i",ig);
      CPlot plotScaleEta_perPt(pname,"","|#eta|","scale factor");
      plotScaleEta_perPt.AddGraph(scale_vs_eta_per_ptv[ig],"", kBlue, kFullDotLarge);
      plotScaleEta_perPt.AddLine(-2.7,1.0,2.7,1.0,kBlack,7);
      plotScaleEta_perPt.AddTextBox(ptlabelv[ig],0.7,0.82,0.9,0.88,0);
      plotScaleEta_perPt.SetYRange(scalelow,scalehigh);
      plotScaleEta_perPt.SetXRange(-2.7,2.7);
      plotScaleEta_perPt.Draw(c,kTRUE,format);
    }
  }
}
