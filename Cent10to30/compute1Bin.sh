# ######
## mass / pt edges and extrapolation region. 
## Set same values in the FitCDFLikelihoodPbPb.C macro
#######
MASSEDGES=(2.6 2.7 2.8 3.2 3.4 3.6)
PTEDGES=(3.0 5.0)
PTLIMITforLIKELIHOOD=2
INTERPREGION=1
CENTMIN=10.
CENTMAX=30.


#######
CANDTYPE="FS"
#CANDTYPE="FF"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[1]},${MASSEDGES[3]},${MASSEDGES[4]},\"FF\",$CENTMIN,$CENTMAX")"
##### fit mass signal & background 
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2SigMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
#exit 

#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[1]},${MASSEDGES[2]},${MASSEDGES[3]},\"FS\",$CENTMIN,$CENTMAX")"
root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},3.22,3.4,\"FF\",$CENTMIN,$CENTMAX")"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[$j]},${MASSEDGES[$j+1]},\"FS\",$CENTMIN,$CENTMAX")"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,5.0,10.0,2.6,2.7,\"FF\;FS\;SS\",$CENTMIN,$CENTMAX")"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,5.0,10.0,2.7,2.8,\"FF\;FS\;SS\",$CENTMIN,$CENTMAX")"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,5.0,10.0,3.2,3.4,\"FF\;FS\;SS\",$CENTMIN,$CENTMAX")"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVA,5.0,10.0,3.4,3.6,\"FF\;FS\;SS\",$CENTMIN,$CENTMAX")"
