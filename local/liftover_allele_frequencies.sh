HG37_FASTA=$HOME/GRCh37/human_g1k_v37.fasta
HG38_FASTA=$HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

for CHROM in {1..22};
do
  python -c "import pandas as pd; pd.read_parquet('~/Downloads/chr${CHROM}_afr_unadmixed_with_rsid.parquet').assign(BP = lambda df: df.SNP.str.split(':', expand = True).iloc[:,1]).assign(SNP = lambda df: df.dbsnp_rsid).to_csv('~/Downloads/chr${CHROM}_afr_unadmixed_with_rsid.tsv', index = False, sep  = '\t')"
  /Users/jor6523/bin/bcftools +munge --no-version -Ou -C colheaders.tsv -f ${HG38_FASTA} ~/Downloads/chr${CHROM}_afr_unadmixed_with_rsid.tsv| \
  /Users/jor6523/bin/bcftools +liftover --no-version -Ou -- -s ${HG38_FASTA} \
    -s ${HG38_FASTA} \
    -f ${HG37_FASTA} \
    -c $HOME/GRCh38/hg38ToHg19.over.chain.gz | \
  /Users/jor6523/bin/bcftools norm -m- | \
  /Users/jor6523/bin/bcftools sort -o ~/Downloads/chr${CHROM}_afr_unadmixed_with_rsid.bcf -Ob --write-index
  /Users/jor6523/bin/bcftools view -H ~/Downloads/chr${CHROM}_afr_unadmixed_with_rsid.bcf > ~/Downloads/chr${CHROM}_afr_unadmixed_with_rsid_lifted.tsv
done
scp ~/Downloads/chr*_afr_unadmixed_with_rsid_lifted.tsv jor6523@o2.hms.harvard.edu:/n/groups/price/jordan/h2xancestry/data/af/afr_unadmixed/
