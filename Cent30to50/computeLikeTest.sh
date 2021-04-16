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

for ((i = 0; i < ${#PTEDGES[@]}-1; ++i)); do
root -b -q LoadLib.C doLikelihoodFitMVA.C"("\"$CANDTYPE\",${PTEDGES[$i]},${PTEDGES[$i+1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},kFALSE,$CENTMIN,$CENTMAX,kFALSE,kTRUE")"
done
