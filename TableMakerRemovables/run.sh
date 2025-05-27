#!/bin/bash
o2-analysis-dq-table-maker -b --configuration json://configuration.json | o2-analysis-timestamp -b --configuration json://configuration.json | o2-analysis-event-selection -b --configuration json://configuration.json | o2-analysis-fwdtrackextension -b --configuration json://configuration.json | o2-analysis-tracks-extra-v002-converter -b --configuration json://configuration.json | o2-analysis-fwdtrack-to-collision-associator -b --configuration json://configuration.json --aod-file AO2D.root

