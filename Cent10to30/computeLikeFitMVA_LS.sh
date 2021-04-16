###
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
POLYNORD=3
LSOPT=2

RUNPTINT=1

 if [[ $RUNPTINT == 1 ]]; then

mv mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}/inputFiles_* .

#WITH LS EVENTS
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2,kTRUE")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2,$POLYNORD,\"LS\"")"
####
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"

CANDTYPE="FF"
#WITH LS EVENTS
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2,kTRUE")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2,$POLYNORD,\"LS\"")"
#
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"

mv LSFiles mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
mv inputFiles_* mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
mv fb*txt mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
 exit

 fi

### pt dep
mv mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}/inputFiles_* .
mv mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}/LSFiles .

CANDTYPE="FF;FS"
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[$i]},${PTEDGES[$i+1]},$CENTMIN,$CENTMAX,1,2,kTRUE")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,1,2,$POLYNORD,\"LS\"")"
###
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
#### cp resolutions for the pt integrated case study
done

CANDTYPE="FF"
for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
root -b -q ExtractJpsiSignal.C"("\"$ANHISTPATH\",${PTEDGES[$i]},${PTEDGES[$i+1]},$CENTMIN,$CENTMAX,2,2,kTRUE")"
root -b -q composeMEandLSBackgroundforLikelihoodFit.C"("${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},$CENTMIN,$CENTMAX,2,2,$POLYNORD,\"LS\"")"
###
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kFALSE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kTRUE,$CENTMIN,$CENTMAX,kTRUE,$LSOPT")"
#### cp resolutions for the pt integrated case study
done

mv LSFiles mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
mv inputFiles_* mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
mv fb*txt mass_${MASSEDGES[0]}_${MASSEDGES[1]}_${MASSEDGES[2]}_${MASSEDGES[3]}_${MASSEDGES[4]}_${MASSEDGES[5]}
