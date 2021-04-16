# ######
## mass / pt edges and extrapolation region. 
## Set same values in the FitCDFLikelihoodPbPb.C macro
#######
MASSEDGES=($1 $2 $3 $4 $5 $6)  
PTEDGES=(1.5 3.0 5.0 10.0)
PTLIMITforLIKELIHOOD=2
#PTEDGES=(5.0 10.0)
INTERPREGION=2
CENTMIN=10.
CENTMAX=30.
ANHISTPATH="/Users/jonare/LikelihoodFit/Cent10to30/InputRootFiles/AnalysisHistograms_jpsi2ee_tigthPID.root"
#######
CANDTYPE="FF;FS"
######
# mass bkg options
1DfitOPT=0 
MEOPT=1
LSOPT=2
POLYNORD=3

cp -r emptyDir/inputFiles* .
binptlimits="${PTEDGES[0]}_${PTEDGES[${#PTEDGES[@]}-1]}"
echo  "$1 $2 $3 $4 $5 $6" >> inputFiles_$binptlimits/massEdges.txt
                       
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2SigMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
#exit
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFF,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FF\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFS,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FS\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolSS,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"SS\",$CENTMIN,$CENTMAX")"

## RUN FOR MIXED EVENT BOTH FF AND FF;FS CANDIDATES
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2,$POLYNORD")"
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2,$POLYNORD")"

## RUN FOR LS EVENT BOTH FF AND FF;FS CANDIDATES
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2,kTRUE")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2,$POLYNORD,\"LS\"")"
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2,kTRUE")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2,$POLYNORD,\"LS\"")"

for ((j =  0; j < ${#MASSEDGES[@]}-1; ++j)); do
  if [[ $j == $INTERPREGION ]]; then
  continue
  fi
  ### to fit the background parameters separately for each candidate type (depending on the available stats.) the last argument 
  ### should be changed accordingly, e.g. 
  root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FF\",$CENTMIN,$CENTMAX")"
  root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FS\",$CENTMIN,$CENTMAX")"
  #root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"SS\"")" 
  #### fit likelihood unbinned
  #root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FF\"")"
  #root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FS\"")"
  #root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"SS\"")" 
  done


  mkdir -p inputFiles_1.5_10.0/bkgPtIntegrated
  cp  inputFiles_1.5_10.0/XBkg*root inputFiles_1.5_10.0/bkgPtIntegrated

root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE")"

#WITH MIXED EVENTS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$MEOPT")"
  
#WITH LS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
  
CANDTYPE="FF"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE")"
#exit
#WITH MIXED EVENTS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$MEOPT")"

#WITH LS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"


CANDTYPE="FF;FS"
## ADD MIXED EVENT FOR PT DEPENDENT CASES

### resolutions vs pt
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
binptlimits="${PTEDGES[$i]}_${PTEDGES[$i+1]}"
echo  "$1 $2 $3 $4 $5 $6" >> inputFiles_$binptlimits/massEdges.txt
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2SigMass,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFF,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FF\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFS,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FS\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolSS,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"SS\",$CENTMIN,$CENTMAX")"
### ME
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[$i]},${PTEDGES[$i+1]},$CENTMIN,$CENTMAX,1,2")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2,$POLYNORD")"
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[$i]},${PTEDGES[$i+1]},$CENTMIN,$CENTMAX,2,2")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2,$POLYNORD")"
### LS
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[$i]},${PTEDGES[$i+1]},$CENTMIN,$CENTMAX,1,2,kTRUE")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2,$POLYNORD,\"LS\"")"
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[$i]},${PTEDGES[$i+1]},$CENTMIN,$CENTMAX,2,2,kTRUE")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2,$POLYNORD,\"LS\"")"
#### cp resolutions for the pt integrated case study
cp inputFiles_$binptlimits/XResol*root inputFiles_1.5_10.0/
done

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

for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
#	if [[ $i == $PTLIMITforLIKELIHOOD ]]; then
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE")"
#WITH MIXED EVENTS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$MEOPT")"
#        fi
#WITH LS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
done

CANDTYPE="FF"
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
#	if [[ $i == $PTLIMITforLIKELIHOOD ]]; then
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE")"
#WITH MIXED EVENTS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$MEOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$MEOPT")"
#WITH LS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
#        fi
done

mkdir mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
mv inputFiles* mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
mv *txt mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
mv MixedEventFiles mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
mv LSFiles mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
##### lkelihood fit to extract fb / fsig
### root LoadLib.C 
### root[1].L FitCDFLikelihoodPbPb.C++
### root[2] FitCDFLikelihoodMVA("FF;FS",1.3,10.,2.2,4.0,kFALSE)  
### if last argument is kTRUE the interpolation is done using linear weights
