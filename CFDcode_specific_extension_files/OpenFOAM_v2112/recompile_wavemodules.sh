#!/bin/sh

rm -R build/linux64GccDPInt32Opt/applications/utilities/preProcessing/setWavesCFDwavemaker
rm -R build/linux64GccDPInt32Opt/src/waveModels
rm -R src/waveModels/Make/linux64GccDPInt32Opt
rm -R applications/utilities/preProcessing/setWavesCFDwavemaker/Make/linux64GccDPInt32Opt

cd src 
./Allwmake waveModels
cd ..
cd applications
./Allwmake utilities/preProcessing/setWavesCFDwavemaker
cd ..