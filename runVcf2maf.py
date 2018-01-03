
import os
import pandas as pd
import synapseclient
syn = synapseclient.login()

vcf2mafPath = "/home/ubuntu/vcf2maf-1.6.14"
veppath = "/home/ubuntu/vep"
vepdata = "/home/ubuntu/.vep"
reference = "/home/ubuntu/.vep/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
tesla_path = "/home/ubuntu/Tesla/TESLA_VCF"
modifiedVCFPath = "/home/ubuntu/Tesla"
ncbibuild = "GRCh38"

mafFiles = []
files = os.walk(tesla_path)
vcffiles = []
for dirpath, dirname, filenames in files:
	vcffiles.extend([os.path.join(dirpath,filename) for filename in filenames if filename.endswith(".vcf")])

for path in vcffiles:
	vcfName = os.path.basename(path)
	metadata = vcfName.split("_")
	team = metadata[1]
	sample = metadata[0]
	newVCFPath = os.path.join(modifiedVCFPath, vcfName)
	os.system("sed 's/^chr//' %s > %s" % (path, newVCFPath))
	os.system("sed -i 's/\t\t/\t.\t/g' %s" % newVCFPath)
	os.system("sed -i 's/ p\./,p./' %s" % newVCFPath)

	if os.path.isfile(newVCFPath+".maf"):
		mafFiles.append(newVCFPath+".maf")
	else:
		os.system(('perl %s/vcf2maf.pl '
				   '--input-vcf %s ' 
				   '--output-maf %s.maf '
				   '--vep-path %s ' 
				   '--vep-data %s '
				   '--vep-forks 8 '
				   '--tumor-id %s '
				   '--maf-center %s '
				   '--custom-enst %s/data/isoform_overrides_uniprot '
				   '--ncbi-build %s '
				   '--ref-fasta %s') % (vcf2mafPath, newVCFPath, newVCFPath, veppath, vepdata, sample, team, vcf2mafPath, ncbibuild, reference))
		if (os.path.isfile(newVCFPath+".maf")):
			mafFiles.append(newVCFPath + ".maf")

maf = pd.DataFrame()
for i in mafFiles:
	temp = pd.read_csv(i,sep="\t",comment="#")
	if len(temp)>0:
		maf = maf.append(temp)
maf.to_csv("TESLA_round2_vcf.maf",sep="\t",index=False)
syn.store(synapseclient.File("TESLA_round2_vcf.maf",parentId = "syn8123644"))
#syn.store(synapseclient.File("uncleanedVCF_vcf2maf_logs.txt",parentId = "syn8123644"))

