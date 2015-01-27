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
### Need to do: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`  # note: path setting added to bashrc
# ./effZFit <conf> <sigpass> <bkgpass> <sigfail> <bkgfail> <infname> <outdir> <doMC> <charge> <temfile> 
#________________________________________________________________________________________________

NOOBDIR=/tthome/bpollack/CMSSW_6_1_1/src/HZG_Analyzer/HiggsZGAnalyzer/tnpFiles


# Brian's settings
#./effZFit sel.bins 3 4 3 2 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root el_ID_data 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root

#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root el_ID_MC 1 0 none 
#./effZFit sel.bins 2 4 2 2 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root el_ID_data 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root

#./effZFit sel.binsLow 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root el_ID_MC_low 1 0 none 
#./effZFit sel.binsLow 2 4 2 2 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root el_ID_data_low 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-11-14_ID.root

##
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-13-14_Iso.root el_Iso_MC 1 0 none 
./effZFit sel.bins 2 4 2 2 ${NOOBDIR}/tnpFile_EE2012ABCD_11-13-14_Iso.root el_Iso_data 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-13-14_Iso.root

#./effZFit sel.binsLow 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-13-14_Iso.root el_Iso_MC_low 1 0 none 
#./effZFit sel.binsLow 2 4 2 2 ${NOOBDIR}/tnpFile_EE2012ABCD_11-13-14_Iso.root el_Iso_data_low 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_11-13-14_Iso.root
##
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_12-04-14_LeadingTrigger.root el_LeadingTrigger_MC 1 0 none 
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_12-04-14_LeadingTrigger.root el_LeadingTrigger_data 0 0 none 
##
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_12-04-14_TrailingTrigger.root el_TrailingTrigger_MC 1 0 none 
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_12-04-14_TrailingTrigger.root el_TrailingTrigger_data 0 0 none 

#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_ID_syst.root el_ID_MC_syst 1 0 none 
#./effZFit sel.bins 1 5 1 1 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_ID_syst.root el_ID_data_syst 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_ID_syst.root
#
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_Iso_syst.root el_Iso_MC_syst 1 0 none 
#./effZFit sel.bins 1 5 1 1 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_Iso_syst.root el_Iso_data_syst 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_Iso_syst.root
#
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_LeadingTrigger_syst.root el_LeadingTrigger_MC_syst 1 0 none 
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_LeadingTrigger_syst.root el_LeadingTrigger_data_syst 0 0 none 
#
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_TrailingTrigger_syst.root el_TrailingTrigger_MC_syst 1 0 none 
#./effZFit sel.bins 0 0 0 0 ${NOOBDIR}/tnpFile_EE2012ABCD_01-15-14_TrailingTrigger_syst.root el_TrailingTrigger_data_syst 0 0 none 
