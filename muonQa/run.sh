#!/bin/bash
echo "Running QA analysis pipeline..."

OPTION="-b --configuration json://configuration.json --aod-file AO2D.root"

# o2-analysis-muon-qa ${OPTION} | o2-analysis-event-selection ${OPTION} | o2-analysis-tracks-extra-v002-converter ${OPTION} | o2-analysis-timestamp ${OPTION} > stdout.log 2>&1

o2-analysis-muon-qa ${OPTION} | o2-analysis-event-selection ${OPTION} | o2-analysis-timestamp ${OPTION} > stdout.log 2>&1
