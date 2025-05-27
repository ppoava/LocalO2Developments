
#! /bin/bash

rm -f stop
rm -f ./AO2D.root
#rm -rf AnalysisResults
mkdir -p AnalysisResults

#AO2DLIST=AO2D_list-LHC24l7.txt
#AO2DLIST=AO2D_list-LHC22p.txt
#AO2DLIST=AO2D_list-LHC23h-apass4_skimmed.txt
#AO2DLIST=AO2D_list-LHC23zk.txt
AO2DLIST=AO2D_list-LHC24am.txt
#AO2DLIST=AO2D_list-LHC24aq-apass1_muon_matching.txt

NFILES=$(cat "${AO2DLIST}" | wc -l)
I=1

while [ $I -le $NFILES ]; do

    if [ -e stop ]; then
        break;
    fi

    if [ -e AnalysisResults/AnalysisResults-${I}.root ]; then
        I=$((I+1))
        continue
    fi

    AO2DFILE=$(cat "${AO2DLIST}" | head -n $I | tail -n 1)
    echo "I=${I} => alien_cp \"alien://${AO2DFILE}\" file://./AO2D.root"
    rm -f ./AO2D.root && alien_cp "alien://${AO2DFILE}" file://./AO2D.root

    echo "I=${I} => processing \"${AO2DFILE}\" ..."
    rm -f AnalysisResults.root && ./run-qa.sh >& AnalysisResults/log-${I}.txt
    gzip -f AnalysisResults/log-${I}.txt
    echo "... done."
    
    cp AnalysisResults.root AnalysisResults/AnalysisResults-${I}.root

    I=$((I+1))

    #break;

    #rm -f ./AO2D.root

done
