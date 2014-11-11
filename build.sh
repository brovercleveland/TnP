#!/bin/bash

echo "Building CPlot..."
g++ -O3 -g -c -Wall -fPIC -o CPlot.o CPlot.cc  `root-config --cflags` `root-config --libs` -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64
g++ -O3 -shared -Wl,-soname,libCPlot.so -o libCPlot.so CPlot.o

echo "Building KStyle..."
g++ -O3 -g -c -Wall -fPIC -o KStyle.o KStyle.cc `root-config --cflags` `root-config --libs` -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64
g++ -O3 -shared -Wl,-soname,libKStyle.so -o libKStyle.so KStyle.o

echo "Building RooCMSShape..."
g++ -O3 -g -c -Wall -fPIC -o RooCMSShape.o RooCMSShape.cc `root-config --cflags` `root-config --libs` -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64
g++ -O3 -shared -Wl,-soname,libRooCMSShape.so -o libRooCMSShape.so RooCMSShape.o

echo "Building RooVoigtianShape..."
g++ -O3 -g -c -Wall -fPIC -o RooVoigtianShape.o RooVoigtianShape.cc `root-config --cflags` `root-config --libs` -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64
g++ -O3 -shared -Wl,-soname,libRooVoigtianShape.so -o libRooVoigtianShape.so RooVoigtianShape.o

echo "Building CEffUser1D..."
g++ -O3 -g -c -Wall -fPIC -o CEffUser1D.o CEffUser1D.cc `root-config --cflags` `root-config --libs` -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64
g++ -O3 -shared -Wl,-soname,libCEffUser1D.so -o libCEffUser1D.so CEffUser1D.o

echo "Building CEffUser2D..."
g++ -O3 -g -c -Wall -fPIC -o CEffUser2D.o CEffUser2D.cc `root-config --cflags` `root-config --libs` -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64
g++ -O3 -shared -Wl,-soname,libCEffUser2D.so -o libCEffUser2D.so CEffUser2D.o

echo "Building CEffZFitter..."
g++ -O3 -g -c -Wall -fPIC -o CEffZFitter.o CEffZFitter.cc `root-config --cflags` `root-config --libs` -l Minuit -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore -L. -l CPlot -l KStyle -l RooCMSShape -l RooVoigtianShape -l CEffUser1D -l CEffUser2D -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64
g++ -O3 -shared -Wl,-soname,libCEffZFitter.so -o libCEffZFitter.so CEffZFitter.o

echo "Building effZFit..."
g++ -O3 -g -o effZFit effZFit.cc `root-config --cflags` `root-config --libs` -l Minuit -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore -L. -l CEffZFitter -l CPlot -l KStyle -l RooCMSShape -l RooVoigtianShape -l CEffUser1D -l CEffUser2D -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64 -l c_nonshared

#echo "Building printEff..."
#g++ -O3 -g -o printEff printEff.cc `root-config --cflags` `root-config --libs` -L. -l CEffUser1D -l CEffUser2D -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64 -l c_nonshared

#echo "Building plotEff..."
#g++ -O3 -g -o plotEff plotEff.cc `root-config --cflags` `root-config --libs` -I $ROOFITSYS/include -L $ROOFITSYS/lib -l RooFit -l RooFitCore -L. -l CPlot -l KStyle -l CEffUser1D -l CEffUser2D -B /usr/lib/x86_64-redhat-linux5E/lib64 -I /usr/lib/x86_64-redhat-linux5E/include -L /usr/lib/x86_64-redhat-linux5E/lib64 -l c_nonshared
