# ######
## mass / pt edges and extrapolation region. 
## Set same values in the FitCDFLikelihoodPbPb.C macro
#######
MASSEDGES=($1 $2 $3 $4)  
PTEDGES=(1.5 3.0 5.0 10.0)
PTLIMITforLIKELIHOOD=2
#PTEDGES=(5.0 10.0)
INTERPREGION=1
CENTMIN=30.
CENTMAX=50.
ANHISTPATH="/Users/jonare/LikelihoodFit/Cent30to50/AnalysisHistograms_jpsi2eeMerged.root"
#######
CANDTYPE="FF;FS"

cp -r emptyDir/inputFiles* .
binptlimits="${PTEDGES[0]}_${PTEDGES[${#PTEDGES[@]}-1]}"
echo  "$1 $2 $3 $4" >> inputFiles_$binptlimits/massEdges.txt


:''
                       
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2SigMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
#exit
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFF,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FF\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFS,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FS\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolSS,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"SS\",$CENTMIN,$CENTMAX")"

## RUN FOR MIXED EVENT BOTH FF AND FF;FS CANDIDATES
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2")"
root -b -q composeMEBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2")"
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2")"
root -b -q composeMEBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2")"

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
#exit
#WITH MIXED EVENTS
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"


CANDTYPE="FF;FS"
## ADD MIXED EVENT FOR PT DEPENDENT CASES

### resolutions vs pt
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
binptlimits="${PTEDGES[$i]}_${PTEDGES[$i+1]}"
echo  "$1 $2 $3 $4" >> inputFiles_$binptlimits/massEdges.txt
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2SigMass,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFF,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FF\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFS,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"FS\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolSS,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"SS\",$CENTMIN,$CENTMAX")"
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[$i]},${PTEDGES[$i+1]},$CENTMIN,$CENTMAX,1,2")"
root -b -q composeMEBackgroundforLikelihoodFit.C"("${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2")"
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[$i]},${PTEDGES[$i+1]},$CENTMIN,$CENTMAX,2,2")"
root -b -q composeMEBackgroundforLikelihoodFit.C"("${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2")"

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
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"
#        fi
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
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,kTRUE")"
#        fi
done

mkdir mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}
mv inputFiles* mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}
mv *txt mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}

##### lkelihood fit to extract fb / fsig
### root LoadLib.C 
### root[1].L FitCDFLikelihoodPbPb.C++
### root[2] FitCDFLikelihoodMVA("FF;FS",1.3,10.,2.2,4.0,kFALSE)  
### if last argument is kTRUE the interpolation is done using linear weights
