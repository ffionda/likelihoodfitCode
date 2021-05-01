MASSEDGES=(2.6 2.7 2.8 3.2 3.4 3.6)
MREGIONS=("L1" "L2" "S" "H1" "H2")
PT="1.5_10.0"
INTERPREGION=2

mkdir startingParametersBkg
for ((j = 0; j < 5; ++j)); do
  if [[ $j == $INTERPREGION ]]; then
  continue
  fi
  cp mass_${MASSEDGES[0]}0_${MASSEDGES[1]}0_${MASSEDGES[2]}0_${MASSEDGES[3]}0_${MASSEDGES[4]}0_${MASSEDGES[5]}0/inputFiles_$PT/bkgPtIntegrated/XBkgParameters_mass${MASSEDGES[j]}_${MASSEDGES[j+1]}_resTypeFF_pt$PT.root startingParametersBkg/XBkgParameters_mass${MREGIONS[j]}_FF_pt$PT.root
  cp mass_${MASSEDGES[0]}0_${MASSEDGES[1]}0_${MASSEDGES[2]}0_${MASSEDGES[3]}0_${MASSEDGES[4]}0_${MASSEDGES[5]}0/inputFiles_$PT/bkgPtIntegrated/XBkgParameters_mass${MASSEDGES[j]}_${MASSEDGES[j+1]}_resTypeFS_pt$PT.root startingParametersBkg/XBkgParameters_mass${MREGIONS[j]}_FS_pt$PT.root
done

