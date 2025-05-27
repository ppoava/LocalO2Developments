
#! /bin/bash

#o2-analysis-dq-task-muon-dca --configuration json://config.json -b | o2-analysis-event-selection --configuration json://config.json -b | o2-analysis-dq-table-maker-with-assoc --configuration json://config.json -b | o2-analysis-tracks-extra-v002-converter --configuration json://config.json -b | o2-analysis-track-propagation --configuration json://config.json -b | o2-analysis-timestamp --configuration json://config.json -b

OPTION="-b --configuration json://config.json --pipeline track-propagation:8 -b"

o2-analysis-dq-qa-muon ${OPTION} | o2-analysis-event-selection ${OPTION} | o2-analysis-fwdtrackextension ${OPTION} | o2-analysis-tracks-extra-v002-converter ${OPTION} | o2-analysis-track-propagation ${OPTION} | o2-analysis-timestamp ${OPTION}

#o2-analysis-dq-qa-muon ${OPTION} | o2-analysis-event-selection ${OPTION} | o2-analysis-fwdtrackextension ${OPTION} | o2-analysis-track-propagation ${OPTION} | o2-analysis-timestamp ${OPTION}

