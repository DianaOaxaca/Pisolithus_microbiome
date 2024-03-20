#!/usr/bin/bash
#DianaOaxaca 230823

#import data
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'\
 --input-path data/manifest.csv \
 --output-path results/01.demux_import_data.qza\
 --input-format 'PairedEndFastqManifestPhred33'

#convert qza to visual file
qiime demux summarize --i-data results/01.demux_import_data.qza --o-visualization results/01.demux_import_data.qzv

#view summarize
#qiime tools view results/demux_import_data.qzv
