# ######
## mass / pt edges and extrapolation region. 
## Set same values in the FitCDFLikelihoodPbPb.C macro
#######
MASSEDGES=(2.5 2.8 3.2 3.6)
PTEDGES=(1.5 10.0)
PTLIMITforLIKELIHOOD=0
INTERPREGION=1
CENTMIN=30.
CENTMAX=50.

#######
CANDTYPE="FF;FS"
#CANDTYPE="FF"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[1]},${MASSEDGES[3]},${MASSEDGES[4]},\"FF\",$CENTMIN,$CENTMAX")"
##### fit mass signal & background 
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2SigMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
#exit; 
#### resolutions vs pt
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFF,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]}")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolFS,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]}")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2XResolSS,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]}")"
done
exit
#### fit x resolution and background (vs pt / mass)
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
  for ((j = 0; j < ${#MASSEDGES[@]}-1; ++j)); do  
  if [[ $j == $INTERPREGION ]]; then
  continue
  fi

  if [[ $i == $PTLIMITforLIKELIHOOD ]]; then
  #### fit all candidate's type simultaneously to determine bkg parameters
  root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FF\;FS\;SS\",$CENTMIN,$CENTMAX")"
  else
  ### to fit the background parameters separately for each candidate type (depending on the available stats.) the last argument 
  ### should be changed accordingly, e.g. 
#  echo "do nothing" 
  root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FF\",$CENTMIN,$CENTMAX")"
  root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FS\",$CENTMIN,$CENTMAX")"
  #root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"SS\",$CENTMIN,$CENTMAX")" 
  fi
  done
done
exit

for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
#       if [[ $i == $PTLIMITforLIKELIHOOD ]]; then
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE")"
#        fi
done
exit
CANDTYPE="FF"
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
#       if [[ $i == $PTLIMITforLIKELIHOOD ]]; then
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE")"
#        fi
done


##### lkelihood fit to extract fb / fsig
### root LoadLib.C 
### root[1].L FitCDFLikelihoodPbPb.C++
### root[2] FitCDFLikelihoodMVA("FF;FS",1.3,10.,2.2,4.0,kFALSE)  
### if last argument is kTRUE the interpolation is done using linear weights
