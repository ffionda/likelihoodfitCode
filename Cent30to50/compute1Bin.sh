# ######
## mass / pt edges and extrapolation region. 
## Set same values in the FitCDFLikelihoodPbPb.C macro
#######
MASSEDGES=(2.5 2.8 3.2 3.6)
PTEDGES=(3.0 5.0)
PTLIMITforLIKELIHOOD=0
INTERPREGION=1
CENTMIN=30.
CENTMAX=50.


#######
CANDTYPE="FS"
#CANDTYPE="FF"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[1]},${MASSEDGES[3]},${MASSEDGES[4]},\"FF\",$CENTMIN,$CENTMAX")"
##### fit mass signal & background 
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2SigMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
#root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitChi2BkgMass,${PTEDGES[0]},${PTEDGES[${#PTEDGES[@]}-1]},${MASSEDGES[0]},${MASSEDGES[${#MASSEDGES[@]}-1]},\"$CANDTYPE\",$CENTMIN,$CENTMAX")"
#exit 

root -b -q LoadLib.C FitCDFLikelihoodPbPb.C++"("kFitXBkgMVAChi2,${PTEDGES[0]},${PTEDGES[1]},${MASSEDGES[2]},${MASSEDGES[3]},\"FS\",$CENTMIN,$CENTMAX")"

