# RBProximity-CLIP-Analysis
# T-to-C conversion filtering script

This script filters BAM alignments for diagnostic T>C (forward strand)
or A>G (reverse strand) substitutions characteristic of PAR-CLIP/fPAR-CLIP
crosslinking events.

## Requirements
- Python 3
- pysam

## Usage
python filter_tc_conversions.py input.bam output.bam

## Notes
Input BAM files must contain MD tags.
