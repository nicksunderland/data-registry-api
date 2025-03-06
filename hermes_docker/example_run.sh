## from command line
#python ./hermes_docker/main.py \
#  --gwaspath "/Users/$(whoami)/Desktop/test_data/gwas_1/test_data.tsv.gz" \
#  --refdir "/Users/$(whoami)/Desktop/test_data/refdata" \
#  --ancestry "EUR" \
#  --phenotype "DCM" \
#  --colmap "varId=rsid,chromosome=chr,position=bp,alt=ea,reference=oa,eaf=eaf,beta=beta,stdErr=se,pValue=p,n=n,ncases=ncase" \
#  --dichotomous \
#  --outputdir "/Users/$(whoami)/Desktop/test_data/gwas_1/output"


# from command line using the docker image
docker run -it \
  -v '/Users/xx20081/Desktop/test_data/refdata:/mnt/refdata' \
  -v '/Users/xx20081/Desktop/test_data/gwas_1:/mnt/testdata' \
  nicksunderland/gwas_qc:latest \
  --gwaspath "/mnt/testdata/test_data.tsv.gz" \
  --refdir "/mnt/refdata" \
  --ancestry "EUR" \
  --phenotype "DCM" \
  --colmap "varId=rsid,chromosome=chr,position=bp,alt=ea,reference=oa,eaf=eaf,beta=beta,stdErr=se,pValue=p,n=n,ncases=ncase" \
  --dichotomous \
  --outputdir "/mnt/testdata/output"