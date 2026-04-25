#!/usr/bin/env bash
set -euo pipefail

plink \
  --file toy \
  --make-bed \
  --out toy \
  --allow-no-sex

echo "Created toy.bed, toy.bim, and toy.fam"
