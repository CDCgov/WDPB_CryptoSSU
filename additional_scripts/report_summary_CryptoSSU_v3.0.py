## Edits made by Shatavia Morrison, xxh5 - 20240916
import os
import re
import pandas as pd
import argparse
from functools import reduce


def read_transposed_tsv(file_path):
    """ Reads and transposes a TSV file where columns are represented as rows. """
    df = pd.read_csv(file_path, sep='\t', header=None, index_col=0).transpose()
    return df


def get_sample_names(samplesheetInfo):
    """ Extracts sample names from samplesheet valid from pipeline info directory. """
    sample_names = pd.DataFrame()
    for filename in os.listdir(samplesheetInfo):
        if filename.endswith(".valid.csv"):
            sample_names = pd.read_csv(os.path.join(samplesheetInfo,filename), sep=',')
            sample_names['sample'] = sample_names['sample'].str.replace("_L001_T1","")

    return sample_names

def get_bbmap_skesa_metrics(bbmap_skesaDir):
    """ Get BBMap metrics for skesa. """
    results_bbmap_skesa = {}
    for filename in os.listdir(bbmap_skesaDir):
        if filename.endswith(".bbmap.log"):
            base = filename.split('.')[0]
            with open(os.path.join(bbmap_skesaDir,filename), 'r') as file:
                content = file.read()
                mapped_reads = re.search(r"Mapped reads:\s+(\d+)",content)
                average_coverage = re.search(r"Average coverage:\s+([\d\.]+)",content)
                number_of_bases = re.search(r"Reads Used:\s+\d+\s+\((\d+) bases\)", content)
                mapped_reads = mapped_reads.group(1) if mapped_reads else "NA"
                average_coverage = average_coverage.group(1) if average_coverage else "NA"
                number_of_bases = number_of_bases.group(1) if number_of_bases else "NA"
                results_bbmap_skesa[base] = {"Number of Reads - Skesa": mapped_reads, "Average Coverage - Skesa": average_coverage, "Number of Bases - Skesa": number_of_bases}
    df = pd.DataFrame.from_dict(results_bbmap_skesa,orient='index')
    df = df.reset_index()
    df = df.rename(columns={'index': 'Sample Name'})
    #print(results_bbmap_skesa)
    return df


def get_bbmap_unicycler_metrics(bbmap_unicyclerDir):
    """ Get BBMap metrics for unicycler. """
    results_bbmap_unicycler = {}
    for filename in os.listdir(bbmap_unicyclerDir):
        if filename.endswith(".bbmap.log"):
            base = filename.split('.')[0]
            with open(os.path.join(bbmap_unicyclerDir,filename), 'r') as file:
                content = file.read()
                mapped_reads = re.search(r"Mapped reads:\s+(\d+)",content)
                average_coverage = re.search(r"Average coverage:\s+([\d\.]+)",content)
                number_of_bases = re.search(r"Reads Used:\s+\d+\s+\((\d+) bases\)", content)
                mapped_reads = mapped_reads.group(1) if mapped_reads else "NA"
                average_coverage = average_coverage.group(1) if average_coverage else "NA"
                number_of_bases = number_of_bases.group(1) if number_of_bases else "NA"
                results_bbmap_unicycler[base] = {"Number of Reads - Unicycler": mapped_reads, "Average Coverage - Unicycler": average_coverage, "Number of Bases - Unicycler": number_of_bases}
    #print(results_bbmap_unicycler)
    df = pd.DataFrame.from_dict(results_bbmap_unicycler, orient='index')
    df = df.reset_index()
    df = df.rename(columns={'index': 'Sample Name'})
    #print(df)
    return df

def get_bbmap_ref_metrics(bbmap_refDir):
    """ Get BBMap metrics for ref. """
    results_bbmap_refDir = {}
    for filename in os.listdir(bbmap_refDir):
        if filename.endswith(".bbmap.log"):
            #print(filename)
            base = filename.split('.')[0]
            with open(os.path.join(bbmap_refDir,filename), 'r') as file:
                content = file.read()
                mapped_reads = re.search(r"Mapped reads:\s+(\d+)",content)
                average_coverage = re.search(r"Average coverage:\s+([\d\.]+)",content)
                number_of_bases = re.search(r"Reads Used:\s+\d+\s+\((\d+) bases\)", content)
                mapped_reads = mapped_reads.group(1) if mapped_reads else "NA"
                average_coverage = average_coverage.group(1) if average_coverage else "NA"
                number_of_bases = number_of_bases.group(1) if number_of_bases else "NA"
                results_bbmap_refDir[base] = {"Number of Reads - withRef": mapped_reads, "Average Coverage - withRef": average_coverage, "Number of Bases - withRef": number_of_bases}
    #print(results_bbmap_refDir)
    df = pd.DataFrame.from_dict(results_bbmap_refDir, orient='index')
    df = df.reset_index()
    df = df.rename(columns={'index': 'Sample Name'})
    #print(df)
    return df




def process_uni_quast(quast_uni_dir):
    """ Processes the Quast Unicycler output, extracting and renaming necessary columns, and cleaning sample names. """
        
    if quast_uni_dir.endswith('quast'):
        quast_uni_dir = quast_uni_dir.rsplit('/', 1)[0]
    # Directly construct the full path to the report.tsv within the assembly/unicycler/quast directory
    report_path = os.path.join(quast_uni_dir,"quast" ,"report.tsv")
    #print(report_path)
    if os.path.exists(report_path):
        df = read_transposed_tsv(report_path)
        required_columns = ["Assembly", "# contigs", "Total length", "N50", "L50", "GC (%)","# N's per 100 kbp"]
        if all(col in df for col in required_columns):
            # Copy the relevant columns for further processing
            df_filtered = df[required_columns].copy()
            # Remove the '.scaffolds' suffix from the entries in the 'Assembly' column
            df_filtered['Assembly'] = df_filtered['Assembly'].apply(lambda x: x.replace('.scaffolds', ''))
            # Rename columns as specified
            df_filtered.rename(columns={
                "Assembly": "Sample Name",
                "# contigs": "Number of Contigs_Unicycler",
                "Total length": "Total Length_Unicycler",
                "L50": "L50_Unicycler",
                "N50": "N50_Unicycler",
                "GC (%)":"GC% Content_Unicycler",
                "# N's per 100 kbp":"# N's per 100 kbp_Unicycler"
            }, inplace=True)
            #df_filtered["Sample Name"]=df_filtered["Sample Name"].str.split('_skesa').str[0]
            #print(df_filtered["Sample Name"])
            return df_filtered
        else:
            return pd.DataFrame(columns=required_columns)
    else:
        print(f"report.tsv not found in {report_path}")
        return pd.DataFrame()

def process_skesa_quast(quast_skesa_dir):
    """ Processes the Quast Unicycler output, extracting and renaming necessary columns, and cleaning sample names. """
        
    if quast_skesa_dir.endswith('quast'):
        quast_skesa_dir = quast_skesa_dir.rsplit('/', 1)[0]
    # Directly construct the full path to the report.tsv within the assembly/unicycler/quast directory
    report_path = os.path.join(quast_skesa_dir,"quast" ,"report.tsv")
    #(report_path)
    if os.path.exists(report_path):
        df = read_transposed_tsv(report_path)
        required_columns = ["Assembly", "# contigs", "Total length", "N50", "L50", "GC (%)","# N's per 100 kbp"]
        if all(col in df for col in required_columns):
            # Copy the relevant columns for further processing
            df_filtered = df[required_columns].copy()
            # Remove the '.scaffolds' suffix from the entries in the 'Assembly' column
            df_filtered['Assembly'] = df_filtered['Assembly'].apply(lambda x: x.replace('.scaffolds', ''))
            # Rename columns as specified
            df_filtered.rename(columns={
                "Assembly": "Sample Name",
                "# contigs": "Number of Contigs_Skesa",
                "Total length": "Total Length_Skesa",
                "L50": "L50_Skesa",
                "N50": "N50_Skesa",
                "GC (%)":"GC% Content_Skesa",
                "# N's per 100 kbp":"# N's per 100 kbp_Skesa"
            }, inplace=True)
            #df_filtered["Sample Name"]=df_filtered["Sample Name"].str.split('_skesa').str[0]
#            print(df_filtered["Sample Name"])
            return df_filtered
        else:
            return pd.DataFrame(columns=required_columns)
    else:
        print(f"report.tsv not found in {report_path}")
        return pd.DataFrame()

def process_kraken(kraken_dir):
#""" Processes Kraken2 results to extract percentage match for 'Cryptosporidium'. """
    kraken2_results={}
    #kraken2_DF=pd.DataFrame()
    for filename in os.listdir(kraken_dir):
        if filename.endswith(".kraken2.report.txt"):
            sample_name = filename.replace(".kraken2.report.txt", "")
            with open(os.path.join(kraken_dir, filename), 'r') as file:
                for line in file:
                    if "Cryptosporidium" in line:
                        parts = line.strip().split('\t')
                        if len(parts) >= 6:
                            percentage_match = parts[0]
                            kraken2_results[sample_name] = percentage_match
                            break
    kraken2_DF = pd.DataFrame(list(kraken2_results.items()),columns=['Sample Name','Kraken Percentage'])
    return kraken2_DF
    
def process_18S(species_18S):
    results_18s = pd.DataFrame()
    for filename in os.listdir(species_18S):
        if filename.endswith("_18S_mqc.csv"):
            results_18s = pd.read_csv(os.path.join(species_18S,filename), sep=',',header=0)
        
          #  'Sample Name', 'query_genome', 'db_bestmatch', 'species', 'pident',
       #'alignment_length', 'coverage', 'query_start', 'query_end',
       #'subject_start', 'subject_end', 'bitscore', 'NCE']
            results_18s.rename(columns ={'species':'species_unicycler','NCE':'NCE_unicycler_18s'},inplace=True)
            
    return results_18s

def process_gp60(species_gp60):
    results_gp60 = pd.DataFrame()
    for filename in os.listdir(species_gp60):
        if filename.endswith("_gp60_mqc.tsv"):
            results_gp60 = pd.read_csv(os.path.join(species_gp60,filename), sep='\t',header=0)
            results_gp60.rename(columns ={'gp60_Subtype':'gp60_Subtype_unicycler','NCE':'NCE_unicycler_gp60'},inplace=True)
    return results_gp60

def process_18S_skesa(species_18S_skesa):
    results_18s = pd.DataFrame()
    for filename in os.listdir(species_18S_skesa):
        if filename.endswith("_18S_mqc.csv"):
            #print(filename)
            results_18s = pd.read_csv(os.path.join(species_18S_skesa,filename), sep=',',header=0)
            #print(results_18s.columns)
            results_18s.rename(columns ={'species':'species_skesa','NCE':'NCE_skesa_18s'},inplace=True)
    return results_18s

def process_gp60_skesa(species_gp60_skesa):
    results_gp60 = pd.DataFrame()
    for filename in os.listdir(species_gp60_skesa):
        if filename.endswith("_gp60_mqc.tsv"):
            #print(filename)
            results_gp60 = pd.read_csv(os.path.join(species_gp60_skesa,filename), sep='\t',header=0)
            #print(results_gp60.columns)
            results_gp60.rename(columns ={'gp60_Subtype':'gp60_Subtype_skesa','NCE':'NCE_skesa_gp60'},inplace=True)

    return results_gp60

def process_kma(kma_results):
    results_kma = pd.DataFrame()
    tempDic = {}
    for filename in os.listdir(kma_results):
        if filename.endswith("_cgmlst.csv"):
           with open(os.path.join(kma_results,filename),'r') as file:
               for line in file:
                   lines = line.rstrip('\n')
                   if lines.startswith('sample_id'):
                      temp = lines.split(',')
                      tempDic[temp[0]]=temp[1:]
                      if temp[0] in tempDic.keys():
                         continue
                   else:
                      temp =lines.split(',')
                      tempDic[temp[0]] =temp[1:]
    results_kma = pd.DataFrame.from_dict(tempDic, orient='index')
    return results_kma

def process_etoki_skesa(etoki_skesa,allAlleles):
    results_etoki_skesa = pd.DataFrame()
    tempDic ={}
    for filename in os.listdir(etoki_skesa):
        if filename.endswith(".etoki.fasta"):
            tempAllAlleles=[]
            yesAlleles=[]
            with open(os.path.join(etoki_skesa,filename),'r') as file:
                sample = filename.split(".etoki.fasta")[0]
                for line in file:
                    temp=[]
                    if line.startswith('>'):
                       alleleInfo = line.split(' ')
                    else: 
                       alleleSeq = line.rstrip('\n')
                       temp.append(sample)
                       temp.append(alleleSeq)
                       nameCheck = alleleInfo[0].replace(">","")
                       if nameCheck in allAlleles:
                          yesAlleles.append(nameCheck)
                          temp.append(alleleInfo[0].replace(">",""))
                          temp.append(alleleInfo[1].replace("value_md5=",""))
                          temp.append(alleleInfo[6].replace("identity=",""))
                          tempAllAlleles.append(temp)
                    tempDic[sample]=tempAllAlleles
            differences = list(set(allAlleles) - set(yesAlleles)) + list(set(yesAlleles) - set(allAlleles))
            noNameSample = filename.split(".etoki.fasta")[0]
            for j in differences:
                tempNo =[]
                tempNo.append(noNameSample)
                tempNo.append('No Sequence')
                tempNo.append(j)
                tempNo.append('No Hash Value')
                tempNo.append('No Percent Identity Allele Match')
                if sample in tempDic:
                    tempDic[sample].append(tempNo)
    columnstartNames = ['Sample Name','Sequence', 'Loci Name','MD5 Hash','Percent Identity']
    allDF = []
    for samples in tempDic:
        df = pd.DataFrame(tempDic[samples], columns = columnstartNames)
        allDF.append(df)    
    results = pd.concat(allDF)
    #print(results)
    return results

def process_etoki_unicycler(etoki_unicycler,allAlleles):
    results_etoki_unicycler = pd.DataFrame()
    tempDic = {}
    for filename in os.listdir(etoki_unicycler):
        if filename.endswith(".etoki.fasta"):
            tempAllAlleles = []
            yesAlleles =[]
            with open(os.path.join(etoki_unicycler,filename),'r') as file:
                sample = filename.split(".etoki.fasta")[0]
                for line in file:
                    temp=[]
                    if line.startswith('>'):
                       alleleInfo = line.split(' ')
                    else: 
                       alleleSeq = line.rstrip('\n')
                       temp.append(sample)
                       temp.append(alleleSeq)
                       nameCheck = alleleInfo[0].replace(">","")
                       if nameCheck in allAlleles:
                          yesAlleles.append(nameCheck)
                          temp.append(alleleInfo[0].replace(">",""))
                          temp.append(alleleInfo[1].replace("value_md5=",""))
                          temp.append(alleleInfo[6].replace("identity=",""))
                          tempAllAlleles.append(temp)
                    tempDic[sample]=tempAllAlleles
            differences = list(set(allAlleles) - set(yesAlleles)) + list(set(yesAlleles) - set(allAlleles))
            noNameSample = filename.split(".etoki.fasta")[0]
            for j in differences:
                tempNo =[]
                tempNo.append(noNameSample)
                tempNo.append('No Sequence')
                tempNo.append(j)
                tempNo.append('No Hash Value')
                tempNo.append('No Percent Identity Allele Match')
                if sample in tempDic:
                    tempDic[sample].append(tempNo)
    columnstartNames = ['Sample Name','Sequence', 'Loci Name','MD5 Hash','Percent Identity']
    allDF = []
    for samples in tempDic:
        df = pd.DataFrame(tempDic[samples], columns = columnstartNames)
        allDF.append(df)    
    results = pd.concat(allDF)
    return results

def process_assem_locate_uni(assembly_location):
    assemDF_uni = pd.DataFrame(columns=['Sample Name','Assembly Location Unicycler'])
    #assemDF_skesa = pd.DataFrame(columns=['Sample Name','Assembly Location Skesa'])
    temp={}
    for filename in os.listdir(assembly_location):
        if filename.endswith(".scaffolds.fa.gz") or filename.endswith(".fasta.gz"):
           pathwayInfo = assembly_location+"/"+filename
           if filename.endswith(".scaffolds.fa.gz"):
              sample = filename.split(".scaffolds.fa.gz")[0]
              temp[sample]=pathwayInfo
    assemDF_uni = pd.DataFrame(temp.items(), columns=['Sample Name','Unicycler Assembly'])
    return assemDF_uni

def process_assem_locate_skesa(assembly_location):
    assemDF_skesa = pd.DataFrame(columns=['Sample Name','Assembly Location Skesa'])
    temp={}
    for filename in os.listdir(assembly_location):
        if filename.endswith(".scaffolds.fa.gz") or filename.endswith(".fasta.gz"):
           pathwayInfo = assembly_location+"/"+filename
           if filename.endswith(".fasta.gz"):
              sample = filename.split(".fasta.gz")[0]
              temp[sample]=pathwayInfo
    assemDF_skesa = pd.DataFrame(temp.items(), columns=['Sample Name','Skesa Assembly'])
    return assemDF_skesa

def get_ref_loci_name(refAlleles):
    refLoci = []
    with open(refAlleles,'r') as file:
        for line in file:
            if line.startswith(">"):
               lines = line.replace(">","")
               loci = lines.split("_")[0]
               refLoci.append(loci)
    return refLoci

def create_initial_dataframe(sample_names):
    """ Creates an initial DataFrame from a list of sample names. """
    df = pd.DataFrame(sample_names, columns=['Sample Name'])
    return df

def generate_summary(sample_batch, output_file,reference_alleles):
    """ Generates a summary report by processing and integrating all data sources. """
    # Base directory for all data
    BASE_DIR = os.path.abspath(sample_batch)
    #print(BASE_DIR)
    base_dir = BASE_DIR
    #print(base_dir)

    # Need reference alleles for Etoki DF set up in a list
    etoki_reference_alleles = get_ref_loci_name(reference_alleles)

    # Create dataframe with samples from samplesheet_valid
    samplesheetInfo = os.path.join(base_dir, "pipeline_info")
    sample_names = get_sample_names(samplesheetInfo)
    #print(sample_names)

    databaseFiles = "/scicomp/groups-pure/OID/NCEZID/DFWED/WDPB/CryptoSSU_WaterGLASS_version_1.0_results/"
    runName =os.path.basename(base_dir)
    #print(databaseFiles,runName)
    print(runName)
    output_file_path = os.path.join(os.getcwd(), output_file)
    #sample_names.to_csv(output_file_path, index=False)
    sample_names['run date'] = runName
    sample_names.to_csv(databaseFiles+runName+"/samplename_table1.csv",index=False) 
    print(sample_names)

    # Process BBMap data ** Need to figure out the metrics to record ***
    bbmap_dir_skesa = os.path.join(base_dir,"bbmap","skesa")
    bbmap_metrics_skesa = get_bbmap_skesa_metrics(bbmap_dir_skesa)
    #print(bbmap_metrics_skesa)
    
    # Process BBMap data ** Need to figure out the metrics to record ***
    bbmap_dir_unicycler = os.path.join(base_dir,"bbmap","unicycler")
    bbmap_metrics_unicycler = get_bbmap_unicycler_metrics(bbmap_dir_unicycler)
    #print(bbmap_metrics_unicycler)

    # Process BBMap data ** Need to figure out the metrics to record ***
    bbmap_dir_ref = os.path.join(base_dir,"bbmap","withReference")
    #print(bbmap_dir_ref)
    bbmap_metrics_ref = get_bbmap_ref_metrics(bbmap_dir_ref)
    #print(bbmap_metrics_ref)

    # Process Assembly data
    quast_uni_dir = os.path.join(base_dir,"assembly","unicycler","quast")
    quast_uni = process_uni_quast(quast_uni_dir)
    #print(quast_uni)

    quast_skesa_dir = os.path.join(base_dir,"assembly","skesa","quast")
    quast_skesa = process_skesa_quast(quast_skesa_dir)
    #print(quast_skesa)
 
     

    kraken_dir = os.path.join(base_dir,"kraken2")
    kraken_results = process_kraken(kraken_dir)
    #print(kraken_results)

    kraken_bbmap_quast = [bbmap_metrics_ref, quast_skesa, quast_uni, kraken_results]
    merged_df = reduce(lambda left, right: pd.merge(left, right, on='Sample Name', how='outer'), kraken_bbmap_quast)
    #print(merged_df['Kraken Percentage'])



    species_18S = os.path.join(base_dir,"18S_species_classification")
    species_18Sresult = process_18S(species_18S)
    #print(species_18Sresult)

    species_gp60 = os.path.join(base_dir,"gp60_subtype_classification")
    species_gp60results = process_gp60(species_gp60)
    #print(species_gp60results)

    species_18S_skesa = os.path.join(base_dir,"18S_species_classification_skesa")
    species_18Sresult_skesa = process_18S_skesa(species_18S_skesa)
    #print(species_18Sresult_skesa)

    species_gp60_skesa = os.path.join(base_dir,"gp60_subtype_classification_skesa")
    species_gp60_skesa_result = process_gp60_skesa(species_gp60_skesa)
    #print(species_gp60) 

    #print(merged_df['Kraken Percentage'])
    #print(species_18Sresult)
    t = pd.merge(merged_df, species_18Sresult[["species_unicycler", "NCE_unicycler_18s","Sample Name"]], on="Sample Name", how="left")
    skesa_18s = species_18Sresult_skesa[["species_skesa", "NCE_skesa_18s","Sample Name"]]
    skesa_gp60 = species_gp60_skesa_result[["gp60_Subtype_skesa","NCE_skesa_gp60","Sample Name"]]
    unicycler_gp60 = species_gp60results[["gp60_Subtype_unicycler","NCE_unicycler_gp60","Sample Name"]]

    regex_pattern = r'([I]+([a-z])+)'
    finalDF = [t, skesa_18s, skesa_gp60, unicycler_gp60]
    csv_df = reduce(lambda left, right: pd.merge(left, right, on='Sample Name', how='outer'), finalDF)
    csv_df['Sample Name'] = csv_df['Sample Name'].str.replace("_L001","")
    csv_df.rename(columns={'Sample Name': 'sample'},inplace=True)
    test = csv_df.gp60_Subtype_unicycler.str.extract(regex_pattern)
    csv_df['Subtype Family'] = csv_df.gp60_Subtype_unicycler.str.extract(regex_pattern)[0]
    csv_df['run date'] = runName
    #print(test[0])
    csv_df.to_csv(databaseFiles+runName+"/summaryresults_table4.csv",index=False)

    print(csv_df)


    kma_MLST = os.path.join(base_dir,"kma")
    kma_results = process_kma(kma_MLST)
    #print(kma_results)

    etoki_skesa_dir= os.path.join(base_dir,"etoki_skesa")
    etoki_skesa_results = process_etoki_skesa(etoki_skesa_dir,etoki_reference_alleles)
    etoki_skesa_results['Sample Name'] = etoki_skesa_results['Sample Name'].str.replace("_L001","")
    etoki_skesa_results.rename(columns={'Sample Name': 'sample'},inplace=True)
    etoki_skesa_results['run date'] =runName
    etoki_skesa_results.to_csv(databaseFiles+runName+"/etokiSkesa_table3a.csv",index=False)
    #print(etoki_skesa_results)
     
    etoki_unicycler_dir = os.path.join(base_dir,"etoki_unicycler")
    etoki_unicycler_results = process_etoki_unicycler(etoki_unicycler_dir,etoki_reference_alleles)
    etoki_unicycler_results['Sample Name'] = etoki_unicycler_results['Sample Name'].str.replace("_L001","")
    etoki_unicycler_results.rename(columns={'Sample Name': 'sample'},inplace=True)
    etoki_unicycler_results['run date'] = runName
    etoki_unicycler_results.to_csv(databaseFiles+runName+"/etokiUnicycler_table3b.csv",index=False)
    #print(etoki_unicycler_results)

    assembly_uni = os.path.join(base_dir,"assembly","unicycler")
    assembly_uni_path = process_assem_locate_uni(assembly_uni)
    #print(assembly_uni_path)

    assembly_skesa = os.path.join(base_dir,"assembly","skesa")
    assembly_skesa_path = process_assem_locate_skesa(assembly_skesa)
    #print(assembly_skesa_path)

    assembly_locations =  assembly_uni_path.merge(assembly_skesa_path, on='Sample Name', how='inner')
    assembly_locations['Sample Name'] = assembly_locations['Sample Name'].str.replace("_L001","")
    assembly_locations.rename(columns={'Sample Name': 'sample'},inplace=True)
    assembly_locations['run date'] = runName
    assembly_locations.to_csv(databaseFiles+runName+"/assemblylocations_table2.csv",index=False) 
    #print(assembly_locations)

    #skesa_dir = os.path.join(base_dir,"assembly","skesa")
    #sample_names_skesa = get_sample_names_from_skesa(skesa_dir)

    # Process BBMap data
    #bbmap_dir = os.path.join(base_dir, "assembly","bbmap","logs")
    #bbmap_results = process_bbmap(bbmap_dir)


    # Save the final DataFrame to CSV
    output_file_path = os.path.join(os.getcwd(), output_file)
    sample_names.to_csv(output_file_path, index=False)
    print(f"Summary report written to {output_file_path}")


if __name__ == "__main__":
#    import argparse
    parser = argparse.ArgumentParser(description="Generate a summary report from various genomic data sources.")
    parser.add_argument("-i", "--sample", required=True, help="Sample batch directory name")
    parser.add_argument("-o", "--output", required=True, help="Output file path for the summary report")
    parser.add_argument("-r", "--refAlleles",required=True,help="Reference Allele files")
    args = parser.parse_args()
    generate_summary(args.sample, args.output, args.refAlleles)
