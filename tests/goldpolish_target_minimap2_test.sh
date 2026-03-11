#!/bin/bash
set -eux -o pipefail

# Download the test data
curl -L --output test_reads.fq https://www.bcgsc.ca/downloads/btl/goldrush/test/test_reads.fq

# Run this demo to test your GoldRush installation
echo "Launching GoldPolish-Target with minimap2"

goldpolish --target --target_dev --minimap2 goldpolish_target_test_golden_path.fa test_reads.fq goldpolish_target_test_golden_path.targeted.fa

if cmp --silent -- goldpolish_target_test_golden_path.targeted.fa $(pwd)/expected_files/minimap2/goldpolish_target_test_golden_path.targeted.polished.expected.fa && \
   cmp --silent -- goldpolish_target_test_golden_path.targeted.gaps.i.fa $(pwd)/expected_files/minimap2/goldpolish_target_test_golden_path.targeted.gaps.expected.fa && \
   cmp --silent -- goldpolish_target_test_golden_path.targeted.gaps.goldpolished.i.fa $(pwd)/expected_files/minimap2/goldpolish_target_test_golden_path.targeted.gaps.goldpolished.expected.fa; then 
    echo "Test successful"
else
  echo "Final polishing file doesn't match expected result - please check your installation"
  exit 1
fi
exit 0