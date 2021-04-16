MASSEDGES=($1 $2 $3 $4)
PTEDGES=(1.5 3.0 5.0 10.0)
PTLIMITforLIKELIHOOD=2
INTERPREGION=2
CENTMIN=30.
CENTMAX=50.
ANHISTPATH="/Users/jonare/LikelihoodFit/Cent10to30/InputRootFiles/AnalysisHistograms_jpsi2eeMerged.root"

CANDTYPE="FF;FS"
cp -r emptyDir/inputFiles_1.5_10.0 .

binptlimits="${PTEDGES[0]}_${PTEDGES[${#PTEDGES[@]}-1]}"
echo  "$1 $2 $3 $4" >> inputFiles_$binptlimits/massEdges.txt
#### fit mass signal & background 
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2SigMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"

for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFF,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FF\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFS,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FS\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolSS,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"SS\",$CENTMIN,$CENTMAX")"
done

## RUN FOR MIXED EVENT BOTH FF AND FF;FS CANDIDATES
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2")"
root -b -q composeMEBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2")"
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2")"
root -b -q composeMEBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2")"

cp -r mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}/inputFiles_1.5_10.0/bkgPtIntegrated inputFiles_1.5_10.0/
#### fit x resolution and background (vs pt / mass)
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
  for ((j =  0; j < ${#MASSEDGES[@]}-1; ++j)); do
  if [[ $j == $INTERPREGION ]]; then
  continue
  fi
  #### fit all candidate's type simultaneously to determine bkg parameters
  ### to fit the background parameters separately for each candidate type (depending on the available stats.) the last argument
  ### should be changed accordingly, e.g.
  if [[ $i == $PTLIMITforLIKELIHOOD ]]; then
  root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FF\;FS\;SS\",$CENTMIN,$CENTMAX")"
  else
  root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FF\",$CENTMIN,$CENTMAX")"
  root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FS\",$CENTMIN,$CENTMAX")"
 #echo "do nothing"
  fi
  #### fit likelihood unbinned
  #root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FF\"")"
  #root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FS\"")"
  #root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"SS\"")"
  done
done
###

root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE")"

#WITH MIXED EVENTS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"

CANDTYPE="FF"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE")"
###

#WITH MIXED EVENTS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"


mv mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}/inputFiles_1.5_10.0 mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}/inputFiles_1.5_10.0_ptIntegratedShapes
mv inputFiles_1.5_10.0 mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}
