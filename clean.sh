VCF=YRI.cov.vcf

# Find the 1-based column index where the first duplicate sample name appears
dup_col=$(awk -F'\t' '
  /^#CHROM/{
    delete seen;
    for (i=10; i<=NF; i++) {
      if ($i in seen) { print i; exit }  # first duplicate position
      seen[$i]=1
    }
    print 0; exit   # no duplicates at all
  }' "$VCF")

# Build the cut field list and run cut (or bgzip)
if [ "$dup_col" -eq 0 ]; then
  echo "No duplicate start found; keeping all columns."
  cut -f 1-"$(awk -F'\t' "/^#CHROM/{print NF; exit}" "$VCF")" "$VCF" > YRI.cov.untildup.vcf
else
  if [ "$dup_col" -le 10 ]; then
    # duplicate starts at first sample -> keep only fixed VCF cols 1..9
    cut -f 1-9 "$VCF" > YRI.cov.untildup.vcf
  else
    # keep 1..9 and 10..(dup_col-1)
    cut -f "1-9,10-$(($dup_col-1))" "$VCF" > YRI.cov.uniq.vcf
  fi
fi