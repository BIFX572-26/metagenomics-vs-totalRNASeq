print("Processing files in python...")

import sys,re,pandas as pd, numpy as np, csv

# Set variable names
SILVA_VERSION="138.2"
NCBI_scientific_input="NCBI_staxids_scientific.txt"
NCBI_non_scientific_input="NCBI_staxids_non_scientific.txt"
taxmap_input=f"taxmap_slv_ssu_lsu_ref_nr_{SILVA_VERSION}.txt"

# Read in NCBI scientific names file as dictionary
NCBI_scientific = open(NCBI_scientific_input,'r')
# Set all file content to lowercase, so that matching of the files later is not
# depending on upper- or lowercase:
NCBI_scientific_lower = (line.lower() for line in NCBI_scientific)
reader1=csv.reader(NCBI_scientific_lower, delimiter='\t')
NCBI_scientific_dict={}
for row in reader1:
	NCBI_scientific_dict[row[0]]=row[1]

# Read in NCBI non-scientific names file as dictionary
NCBI_non_scientific = open(NCBI_non_scientific_input,'r')
# Set all file content to lowercase, so that matching of the files later is not
# depending on upper- or lowercase:
NCBI_non_scientific_lower = (line.lower() for line in NCBI_non_scientific)
reader2=csv.reader(NCBI_non_scientific_lower, delimiter='\t')
NCBI_non_scientific_dict={}
for row in reader2:
	NCBI_non_scientific_dict[row[0]]=row[1]

# ======  1. makeblast

# Read in SILVA file as df
taxmap_df=pd.read_csv(taxmap_input, delimiter = "\t")
# Extract second last rank of taxonomy paths and turn them into df col "second_last_rank"
second_last_rank_lst=taxmap_df["path"].str.split(";")
taxmap_df["second_last_rank"]=second_last_rank_lst.str[-2]
# Extract organism name of taxmap df, replace unwanted characters,
# keep only first two words (idealy genus+species) and turn them into df col "species"
name=taxmap_df["organism_name"].str.replace(r'[-()\[\]\\/]', '', regex=True).str.split()
taxmap_df['species'] = name.str[0].str.cat(name.str[1], sep=" ", na_rep='' ).str.strip()
# Replace all species names by NA if their first word (ideally genus) doesn't
# match the second last rank (which should be genus ideally)
taxmap_df['species'] = np.where(taxmap_df['second_last_rank'] == taxmap_df['species'].str.split().str[0], taxmap_df['species'], "NA")
# Merge path and species
taxmap_df["path"] = taxmap_df["path"] + taxmap_df["species"]
# Add ; at end of path to keep format consistent
taxmap_df["path"]=taxmap_df["path"] + ";"
# Keep all rows that contain species (used later)
taxmap_df_sp=taxmap_df.loc[taxmap_df['species'] != "NA"]
# Drop unwanted cols
taxmap_df=taxmap_df.drop(["organism_name", "taxid", "second_last_rank", "species"], axis=1)
# Turn all path leters into lowercase
taxmap_df["path"] = taxmap_df["path"].str.lower()
# Turn df into dict
taxmap_dict={}
for i in range(len(taxmap_df)):
	taxmap_dict[taxmap_df["primaryAccession.start.stop"][i]]=taxmap_df["path"][i]
# Split up SILVA taxonomy and invert it
for key,value in taxmap_dict.items():
	taxmap_dict[key] = value.split(";")[::-1][1:]

# Loop taxmap taxonomy dictionary over NCBI taxonomy dictionary for each line and
# match with NCBI staxid when a hit is found
output_dict={} # Make empty dictionary for matching lines
for key,value in taxmap_dict.items():
	l=0
	while l < len(value):
		if value[l] in NCBI_scientific_dict:
			output_dict[key.upper()]=NCBI_scientific_dict[value[l]]
			hit=True
			break
		elif value[l] in NCBI_non_scientific_dict:
			output_dict[key.upper()]=NCBI_non_scientific_dict[value[l]]
			hit=True
			break
		else:
			l += 1
	if hit==False:
		output_dict[key.upper()]='0'

# Saves the merged dictionary in a defined output file, delimited by tab
with open("SILVA_accession_numbers_and_NCBI_taxids.txt", 'w') as f:
    for key in output_dict.keys():
        f.write("%s\t%s\n"%(key,output_dict[key]))


# ======= 2. kraken2

# Set variable names:
tax_slv_input=f"tax_slv_ssu_{SILVA_VERSION}.txt"

# Make nodes.dmp:

## Read in tax file (note: we don't need the LSU tax file as well since all taxa
## from the LSU tax file are also part of the SSU tax file):
tax_slv_df_child=pd.read_csv(tax_slv_input, sep='\t', names=['path','taxID','rank','remark','release'], usecols=['path','taxID','rank'])
## Add paths from species and set their rank to species
tax_slv_df_child_sp=pd.concat([tax_slv_df_child, taxmap_df_sp[["path"]]])
tax_slv_df_child_sp['rank']=tax_slv_df_child_sp['rank'].fillna("species")
## Get NCBI ID for every child and parent path, this now includes all ranks
child_ncbi_id=[]
parent_ncbi_id=[]
for child_path in tax_slv_df_child_sp["path"]:
	parent_path=re.sub('[^;]*;$', '', child_path)
	child_lst=child_path.lower().split(";")[::-1][1:]
	parent_lst=parent_path.lower().split(";")[::-1][1:]
	l=0
	hit_child=False
	while l < len(child_lst):
		if child_lst[l] in NCBI_scientific_dict:
			child_id=NCBI_scientific_dict[child_lst[l]]
			hit_child=True
			break
		elif child_lst[l] in NCBI_non_scientific_dict:
			child_id=NCBI_non_scientific_dict[child_lst[l]]
			hit_child=True
			break
		else:
			l += 1
	if hit_child==False:
		child_id='1'
	m=0
	hit_parent=False
	while m < len(parent_lst):
		if parent_lst[m] in NCBI_scientific_dict:
			if NCBI_scientific_dict[parent_lst[m]]==child_id:
				m+=1
				continue
			parent_id=NCBI_scientific_dict[parent_lst[m]]
			hit_parent=True
			break
		elif parent_lst[m] in NCBI_non_scientific_dict:
			if NCBI_non_scientific_dict[parent_lst[m]]==child_id:
				m+=1
				continue
			parent_id=NCBI_non_scientific_dict[parent_lst[m]]
			hit_parent=True
			break
		else:
			m += 1
	if hit_parent==False:
		parent_id='1'
	child_ncbi_id.append(child_id)
	parent_ncbi_id.append(parent_id)

## Add NCBI IDs to df
tax_slv_df_child_sp["child_ncbi_id"]=child_ncbi_id
tax_slv_df_child_sp["parent_ncbi_id"]=parent_ncbi_id
## Generate column indicating # of ranks
tax_slv_df_child_sp["numranks"]=tax_slv_df_child_sp["path"].str.count(";")
## Add a row with root as first line, otherwise we can't make the nodes.dmp
## file as there is no parent ID for domains:
tax_slv_df_child_sp.loc[-1] = ['root;', 1, 'no rank', '1', '1', 1]
## Only keep one path per child_ncbi_id, therefore sort df by child_ncbi_id and
## then numranks and keep only first occurence of child_ncbi_id
## (= path with the lowest # of ranks per child_ncbi_id)
tax_slv_df_child_sp=tax_slv_df_child_sp.sort_values(['child_ncbi_id', 'numranks']).drop_duplicates(subset=['child_ncbi_id']).reset_index().drop(["index"], axis=1)
## Columns have to be separated by "\t|\t" or "\t-\t", which doesn't work as delimiter, so
## we generate a spacer vector made of "|" and "-" and a spacer for the last
## column containing "scientific name" to fit the format:
spacer1=np.repeat("|", len(tax_slv_df_child_sp), axis=0)
spacer2=np.repeat("-", len(tax_slv_df_child_sp), axis=0)
spacer3=np.repeat("scientific name", len(tax_slv_df_child_sp), axis=0)

## Generate the final dataframe
nodes_SILVA_SSU_LSU_df=pd.DataFrame({'child_tax_id':tax_slv_df_child_sp["child_ncbi_id"], 'spacer1':spacer1, \
'parent_tax_id':tax_slv_df_child_sp["parent_ncbi_id"], 'spacer2':spacer1, 'rank':tax_slv_df_child_sp["rank"], \
'spacer3':spacer1, 'spacer4':spacer2, 'spacer5':spacer1})

## Write df
nodes_SILVA_SSU_LSU_df.to_csv("nodes.dmp", sep='\t', na_rep='NA', index=False, \
header=False)


# Make names.dmp
## Get last rank of every path and remove the last ";"
last_rank=tax_slv_df_child_sp["path"].str.replace(r';$', '').str.replace(r'^.*;', '')

## Generate the final dataframe
names_SILVA_SSU_LSU_df=pd.DataFrame({'tax_id':tax_slv_df_child_sp["child_ncbi_id"], 'spacer1':spacer1, \
'taxon':last_rank, 'spacer2':spacer1, 'spacer3':spacer2, 'spacer4':spacer1, \
'spacer5':spacer3, 'spacer6':spacer1})

## Write df
names_SILVA_SSU_LSU_df.to_csv("names.dmp", sep='\t', na_rep='NA', index=False, \
header=False)

print("Python done!")
