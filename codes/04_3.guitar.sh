#!/bin/bash
mkdir -p ./data/04.Tx_quantification/03.Guitar
# cp all files in data/04.Tx_quantification/02.Stringtie_with_impute/bam_impute/GS689_2Aligned.sortedByCoord.out/ end with transcript.expression.txt
cp ./data/04.Tx_quantification/02.Stringtie_with_impute/bam_impute/*/*transcript.expression.txt ./data/04.Tx_quantification/03.Guitar/
