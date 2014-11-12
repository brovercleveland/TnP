#!/bin/bash

#================================================================================================
#
# Signal Extraction
#-------------------
#  0: probe counting
#  1: Breit-Wigner convolved with Crystal Ball function
#  2: MC template convolved with Gaussian
#  3: Phil's Crystal Ball based "Voigtian" shape
#  4: Unbinned MC data convolved with Gaussian
#
# Background Model
#------------------
#  0: no background
#  1: exponential model
#  2: erfc*exp model
#  3: double exponential model
#  4: linear*exp model
#  5: quadratic*exp model
#________________________________________________________________________________________________
#
#================================================================================================
#
# effZFit usage
#---------------
### Need to do: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
# ./effZFit <conf> <sigpass> <bkgpass> <sigfail> <bkgfail> <infname> <outdir> <doPU> <charge> <temfile> 
#________________________________________________________________________________________________

#BACONBITSDIR=/tthome/ksung/cms/Analysis/01d/CMSSW_5_3_14_patch2/src/DMSAna/TagAndProbe/macros
NOOBDIR=/tthome/bpollack/CMSSW_6_1_1/src/HZG_Analyzer/HiggsZGAnalyzer/tnpFiles

#./effZFit smu.bins 0 0 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_smubits.root smu_mc   1 0 none
#./effZFit smu.bins 0 0 0 0 ${BACONBITSDIR}/SingleMu_2012-22Jan2013_smubits.root             smu_data 0 0 none

#./effZFit musel.bins 0 0 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_muselbits.root musel_mc   1 0 none
#./effZFit musel.bins 2 1 2 2 ${BACONBITSDIR}/SingleMu_2012-22Jan2013_muselbits.root             musel_data 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_muselbits.root

#./effZFit muiso.bins 0 0 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_muisobits.root muiso_mc   1 0 none
#./effZFit muiso.bins 2 1 2 2 ${BACONBITSDIR}/SingleMu_2012-22Jan2013_muisobits.root             muiso_data 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_muisobits.root

#./effZFit muid.bins 0 0 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_muidbits.root muid_mc   1 0 none
#./effZFit muid.bins 2 1 2 2 ${BACONBITSDIR}/SingleMu_2012-22Jan2013_muidbits.root             muid_data 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_muidbits.root


#./effZFit sel.bins 0 0 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_selbits.root sel_mc   1 0 none
#./effZFit sel.bins 0 0 0 0 ${BACONBITSDIR}/SingleElectron_2012-22Jan2013_selbits.root       sel_data 0 0 none

#./effZFit elsel.bins 0 0 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_elselbits.root elsel_mc   1 0 none
#./effZFit elsel.bins 2 1 2 2 ${BACONBITSDIR}/SingleElectron_2012-22Jan2013_elselbits.root       elsel_data 0 0 ${BACONBITSDIR}/Summer12_DYJetsToLL_M-50_TuneZ2Star_elselbits.root

#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root elsel_MC 1 0 none 
./effZFit sel.bins 2 1 2 2 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root elsel_data 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root
