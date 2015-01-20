#! /usr/bin/env python

import csv


def getVals(path):
  output = []
  with open(path+'/summary.txt') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in spamreader:
      if 'Efficiency' in row: continue
      if '-----------------' in row: continue
      if len(row)==0: break
      if '' == row[0]: continue
      for i,number in enumerate(row):
        #print number
        if(']' not in number): continue
        #print ' '.join(row[i+1:-1])
        #print row[i+1:-1]
        output.append([float(x) for x in row[i+1:-1] if x !=''])
        break

      #print ' '.join(row)
      #print row[2]
  return output

def outputDiffs(realVals, systVals):
  ptBins = ['[  10,  20]','[  20,  40]','[  40,  60]','[  60, 100]','[ 100,8000]']
  if len(realVals) != len(systVals): raise Exception('realVals and systVals not same length what the fuck?!')
  print '[   0, 0.8][ 0.8,1.4442][1.4442,1.556][1.556,   2][   2, 2.5]'
  for i,row in enumerate(realVals):
    print(ptBins[i]),
    for j,val in enumerate(row):
      diff = round(100*abs(realVals[i][j]-systVals[i][j])/realVals[i][j],3)
      print(str(diff)+'     '),
    print

if __name__ == "__main__":
  realPaths = ['el_ID_MC','el_Iso_MC','el_LeadingTrigger_MC','el_TrailingTrigger_MC']
  systPaths = ['el_ID_MC_syst','el_Iso_MC_syst','el_LeadingTrigger_MC_syst','el_TrailingTrigger_MC_syst']

  for i, path in enumerate(realPaths):
    realVals = getVals(realPaths[i])
    systVals = getVals(systPaths[i])

    print realPaths[i].split('_')[1]
    outputDiffs(realVals, systVals)
    print



