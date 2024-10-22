import re
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statistics as st

def intersection(*lists):
# Finds the intersection between multiple sets
    intersection = set(lists[0]).intersection(*lists[1:])
    
    return intersection

def union(*lists):
# Finds the union between multiple sets
    union = set(lists[0]).union(*lists[1:])
    
    return union  

def missing_gene_in_gpls(difference, gpl_dict):
# Returns a dictionary saying how many genes are missing in each GPL
    missing_gene_count = {}
    for gpl_name, gpl_contents in gpl_dict.items():
        missing_gene_count[gpl_name] = 0
        for gene in difference:
            if gene not in gpl_contents:
                missing_gene_count[gpl_name] += 1
                
    return missing_gene_count

def remove_duplicates(list):
# Takes a list in as an argument and returns a new list with its duplicates removed
    new_list = []
    for item in list:
        if item not in new_list:
            new_list.append(item)
            
    return new_list

def rv_genome():
# Returns the list of Rv numbers found in the genome annotation
    rv_genome_list = []
    
    with open("RVintersection/ncbi_dataset1.tsv") as f:
        firstLine = True
        for line in f:
            data = line.split('\t')
            locus_tag = data[14]
        
            if firstLine:
                firstLine = False
                continue
            else:
                rv_genome_list.append(locus_tag.strip())
    
        return rv_genome_list
    
def total_sample_gene_count():
# Gets the count of Rv genes from all samples
    sample_counts = []
    directory = "RVintersection"
    for series in os.listdir(directory): # iterates through each series folder
        if series.startswith('GSE'):
            series_path = f"RVintersection/{series}"
            gpl_gene_count = 0
            
            # First we try look for the Rv genes in the series, since each sample shares the same # of Rv genes as the GPL file
            for sample in os.listdir(series_path):
                if sample.startswith('GPL'):
                    gpl_path = f"RVintersection/{series}/{sample}"
                    gpl_gene_count = len(rv_extraction_complete(gpl_path))
                else:
                    continue
            
            # Then we count how many samples there are
            for sample in os.listdir(series_path):
                if sample.startswith(('log2', 'GSM')):
                    sample_counts.append(gpl_gene_count)
                else:
                    continue
        else:
            continue
        
    return sample_counts

"""
def rv_intersection_matrix(intersection_set):
# Parses through all samples and collects the metadata and associates it with the Rv intersection set to create a matrix
    index = []
    intersection_dict = {}
    for rv_id in intersection_set:
        index.append(rv_id)
        
    directory = "RVintersection"
    for series in os.listdir(directory):
        series_path = f"RVintersection/{series}"
        
        if "GPL7477-tbl-1.txt" in os.listdir(series_path):
            mtub_hash_map = rv_num_hash_mtub_num_GPL7477()
            for sample in os.listdir(series_path):
                sample_name = sample[sample.find("GSM"):sample.find("-")]
                intersection_dict[sample_name] = []
                
                sample_path = f"RVintersection/{series}/{sample}"
                with open(sample_path) as file:
                    for line in file:
                        data = line.split(' ')
                        mtub_tag = data[0]
                        rv_tag = mtub_hash_map[mtub_tag]
                        meta_data = data[1]
                        
                        intersection_dict[sample_name].append()
                
            return
            
    
    rv_matrix = pd.DataFrame(data=intersection_dict, index=index)
    return rv_matrix
"""
            
    
def symbols_hash_rv_num():
# Hashes symbols to the Rv number and stores it in a dictionary
    hash_map = {}
    
    with open("RVintersection/ncbi_dataset1.tsv") as f:
        firstLine = True
        for line in f:
            data = line.split('\t')
            symbol = data[6]
            locus_tag = data[14]
        
            if firstLine:
                firstLine = False
                continue
            else:
                if symbol != '':
                    hash_map[symbol] = locus_tag.strip()
                else:
                    pass
        
        return hash_map
    
def gene_id_hash_rv_num():
# Hashes gene ID to the Rv number and stores it in a dictionary
    hash_map = {}
    
    with open("RVintersection/ncbi_dataset1.tsv") as f:
        firstLine = True
        for line in f:
            data = line.split('\t')
            gene_id = data[7]
            locus_tag = data[14]
        
            if firstLine:
                firstLine = False
                continue
            else:
                hash_map[gene_id] = locus_tag.strip()
        
        return hash_map
    
def genomic_location_hash_rv_num():
# Hashes genomic location to the Rv number and stores it in a dictionary
    hash_map = {}
    
    with open("RVintersection/ncbi_dataset1.tsv") as f:
        firstLine = True
        for line in f:
            data = line.split('\t')
            accession = data[0]
            begin = data[1]
            end = data[2]
            
            genomic_location = f"{accession}:{begin}-{end}"
            locus_tag = data[14]
        
            if firstLine:
                firstLine = False
                continue
            else:
                hash_map[genomic_location] = locus_tag.strip()
        
        return hash_map

def rv_extraction_complete(file):
# Combines all the Rv extraction methods below into one function while ignoring duplicates
    rv_nums = [] 
    rv_genome_list = rv_genome()
    genomic_location_hash_map = genomic_location_hash_rv_num()
    symbol_hash_map = symbols_hash_rv_num() 
    gene_id_hash_map = gene_id_hash_rv_num()
    
    # Tries to search for an Rv tag first
    with open(file) as f:
        for line in f:
            # Tries to search for an Rv number first
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            try:
                if rv[0] not in rv_nums and rv[0] in rv_genome_list:
                    rv_nums.append(rv[0])
                    continue
            except:
                pass
            
            # If that doesn't work we look for a genomic location
            genomic_location = re.findall("[c0-9]+-[0-9]+", line)
            try:
                # If 'c' is in the genomic location, the beginning and end numbers are flipped for some reason
                if 'c' in genomic_location[0]: 
                    c_index = genomic_location[0].find('c')
                    dash_index = genomic_location[0].find('-')
                    
                    begin_num = genomic_location[0][dash_index+1:]
                    end_num = genomic_location[0][c_index+1:dash_index]
                    
                    accession = f"NC_000962.3:{begin_num}-{end_num}"
                    if genomic_location_hash_map[accession] not in rv_nums:
                        rv_nums.append(genomic_location_hash_map[accession])
                        continue
                    
                # You don't have to flip the numbers if there's no 'c'
                else: 
                    accession = f"NC_000962.3:{genomic_location[0]}"
                    if genomic_location_hash_map[accession] not in rv_nums:
                        rv_nums.append(genomic_location_hash_map[accession])
                        continue
            except:
                pass
            
            # If that doesn't work either, then we try looking for a symbol
            for symbol, rv in symbol_hash_map.items():
                if symbol in line and rv not in rv_nums:
                    rv_nums.append(rv)
                    
            # If nothing above works, we look for a gene ID
            for gene_id, rv in gene_id_hash_map.items():
                if gene_id in line and rv not in rv_nums:
                    rv_nums.append(rv)          
                    
    return rv_nums

def duplicates(file):
# Counts the number of duplicates and unique duplicates for a series
    duplicates = 0
    unique_duplicates = []
    rv_nums = [] 
    rv_genome_list = rv_genome()
    genomic_location_hash_map = genomic_location_hash_rv_num()
    symbol_hash_map = symbols_hash_rv_num() 
    gene_id_hash_map = gene_id_hash_rv_num()
    
    # Tries to search for an Rv tag first
    with open(file) as f:
        for line in f:
            # Tries to search for an Rv number first
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            try:
                if rv[0] not in rv_nums and rv[0] in rv_genome_list:
                    rv_nums.append(rv[0])
                    continue
                elif rv[0] in rv_nums and rv[0] in rv_genome_list:
                    duplicates += 1
                    if rv[0] not in unique_duplicates:
                        unique_duplicates.append(rv[0])
                    continue
            except:
                pass
            
            # If that doesn't work we look for a genomic location
            genomic_location = re.findall("[c0-9]+-[0-9]+", line)
            try:
                # If 'c' is in the genomic location, the beginning and end numbers are flipped for some reason
                if 'c' in genomic_location[0]: 
                    c_index = genomic_location[0].find('c')
                    dash_index = genomic_location[0].find('-')
                    
                    begin_num = genomic_location[0][dash_index+1:]
                    end_num = genomic_location[0][c_index+1:dash_index]
                    
                    accession = f"NC_000962.3:{begin_num}-{end_num}"
                    if genomic_location_hash_map[accession] not in rv_nums:
                        rv_nums.append(genomic_location_hash_map[accession])
                        continue
                    else:
                        dupilcates += 1
                        if genomic_location_hash_map[accession] not in unique_duplicates:
                            unique_duplicates.append(genomic_location_hash_map[accession])
                        continue
                    
                # You don't have to flip the numbers if there's no 'c'
                else: 
                    accession = f"NC_000962.3:{genomic_location[0]}"
                    if genomic_location_hash_map[accession] not in rv_nums:
                        rv_nums.append(genomic_location_hash_map[accession])
                        continue
                    else:
                        duplicates += 1
                        if genomic_location_hash_map[accession] not in unique_duplicates:
                            unique_duplicates.append(genomic_location_hash_map[accession])
                        continue
            except:
                pass
            
            # If that doesn't work either, then we try looking for a symbol
            symbolFound = False
            for symbol, rv in symbol_hash_map.items():
                if symbol in line and rv not in rv_nums:
                    rv_nums.append(rv)
                    symbolFound = True
                    break
                elif symbol in line and rv in rv_nums:
                    duplicates += 1
                    if rv not in unique_duplicates:
                        unique_duplicates.append(rv)
                        symbolFound = True
                    break
            
            if symbolFound:
                continue
                    
            # If nothing above works, we look for a gene ID
            for gene_id, rv in gene_id_hash_map.items():
                if gene_id in line and rv not in rv_nums:
                    rv_nums.append(rv) 
                    break
                elif gene_id in line and rv in rv_nums:
                    duplicates += 1    
                    if rv not in unique_duplicates:
                        unique_duplicates.append(rv)   
                        break  
                    
    return (duplicates, len(unique_duplicates))
      
def rv_extraction(file):
# Extracts the Rv number from a sample
    rv_nums = []
    rv_genome_list = rv_genome()
    
    with open(file) as f:
        for line in f:
            # Tries to search for an Rv number first
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            try:
                if rv[0] not in rv_nums and rv[0] in rv_genome_list:
                    rv_nums.append(rv[0])
                    continue
            except:
                pass
                
    return rv_nums

def gene_id_to_rv_extraction(file):
# Extracts the Rv number in a sample from a gene ID
    rv_nums = []
    gene_id_hash_map = gene_id_hash_rv_num()
    
    with open(file) as f:
        for line in f:
            for gene_id, rv in gene_id_hash_map.items():
                if gene_id in line and rv not in rv_nums:
                    rv_nums.append(rv)          
            
    return rv_nums

def symbol_to_rv_extraction(file):
# Extracts the Rv number in a sample from a symbol
    rv_nums = []
    symbol_hash_map = symbols_hash_rv_num()
    
    with open(file) as f:
        for line in f:
            for symbol, rv in symbol_hash_map.items():
                if symbol in line and rv not in rv_nums:
                    rv_nums.append(rv)
            
    return rv_nums

def genomic_location_to_rv_extraction(file):
# Extracts the Rv number from a sample but from its genomic location
    rv_nums = []
    genomic_location_hash_map = genomic_location_hash_rv_num()
    
    with open(file) as f:
        for line in f:
            genomic_location = re.findall("[c0-9]+-[0-9]+", line)
            try:
                # If 'c' is in the genomic location, the beginning and end numbers are flipped for some reason
                if 'c' in genomic_location[0]:
                    c_index = genomic_location[0].find('c')
                    dash_index = genomic_location[0].find('-')
                    
                    begin_num = genomic_location[0][dash_index+1:]
                    end_num = genomic_location[0][c_index+1:dash_index]
                    
                    accession = f"NC_000962.3:{begin_num}-{end_num}"
                    if genomic_location_hash_map[accession] not in rv_nums:
                        rv_nums.append(genomic_location_hash_map[accession])
                        continue
                # You don't have to flip the numbers if there's no 'c'
                else:
                    accession = f"NC_000962.3:{genomic_location[0]}"
                    if genomic_location_hash_map[accession] not in rv_nums:
                        rv_nums.append(genomic_location_hash_map[accession])
                        continue
            except:
                pass
            
    return rv_nums
            
def rv_num_hash_mtub_num_GPL7477():
# Hashes Rv num to MTUB num in the series GPL7477
    hash_map = {}
            
    with open("RVintersection/GSE13246_family.xml/GPL7477-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    
                    hash_map[tag] = rv[0]
                except:
                    pass
    return hash_map

def rv_num_hash_qmt_num_GPL5773():
# Hashes Rv num to 4QMT num in the series GPL5773
    hash_map = {}
            
    with open("RVintersection/GSE23498_family.xml/GPL5773-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    
                    hash_map[tag] = rv[0]
                except:
                    pass
    return hash_map
    
def rv_num_hash_cust_num_GPL15565():
# Hashes Rv num to CUST num in the series GPL15565
    hash_map = {}
            
    with open("RVintersection/GSE37973_family.xml/GPL15565-tbl-1.txt") as file:
         for line in file:
            columns = line.split('\t')
            tag = columns[0]
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    
                    hash_map[tag] = rv[0]
                except:
                    pass
    return hash_map

def rv_num_hash_cust_num_GPL14824():
# Hashes Rv num to CUST num in the series GPL14824
    hash_map = {}
            
    with open("RVintersection/GSE40917_family.xml/GPL14824-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    
                    hash_map[tag] = rv[0]
                except:
                    pass
    return hash_map

def rv_num_hash_sym_num_GPL4519():
# Hashes Rv num to symbol in the series GPL4519
    hash_map = {}
            
    with open("RVintersection/GSE43749_family.xml/GPL4519-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    
                    hash_map[tag] = rv[0]
                except:
                    pass
    return hash_map

def rv_num_hash_tbid_num_GPL17082():
# Hashes Rv num to TBID in the series GPL17082
    hash_map = {}
            
    with open("RVintersection/GSE46432_family.xml/GPL17082-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    
                    hash_map[tag] = rv[0]
                except:
                    pass
    return hash_map

def rv_num_hash_jiehe_num_GPL18119():
# Hashes Rv num to jiehe num in the series GPL18119
    hash_map = {}
    gene_id_hash_map = gene_id_hash_rv_num()
            
    with open("RVintersection/GSE53843_family.xml/GPL18119-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]

            try:
                gene_id = columns[3]
                rv = gene_id_hash_map[gene_id]
                hash_map[tag] = rv
            except:
                pass
            
    return hash_map
    
def rv_num_hash_sym_num_GPL18055():
# Hashes Rv num to symbol in the series GPL18055
    hash_map = {}
    symbols_hash = symbols_hash_rv_num()         
    with open("RVintersection/GSE60376_family.xml/GPL18055-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    
                    hash_map[tag] = rv[0]
                except:
                    pass
            else:
                try:
                    hash_map[tag] = symbols_hash[tag]
                except:
                    pass
    return hash_map

def rv_num_hash_tb_num_rv_GPL22695():
# Hashes Rv num to TB num in the series GPL22695
    hash_map = {}
    genomic_location_hash = genomic_location_hash_rv_num()        
    with open("RVintersection/GSE90858_family.xml/GPL22695-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]
            accession = "NC_000962.3"
            begin_and_end = tag[tag.find(":")+1:]
            if "c" in begin_and_end:
                begin_num = begin_and_end[begin_and_end.find("-")+1:]
                end_num = begin_and_end[begin_and_end.find("c")+1:begin_and_end.find("-")]
                begin_and_end = f"{begin_num}-{end_num}"
            else:
                pass
            genomic_location = f"{accession}:{begin_and_end}"
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    hash_map[tag] = rv[0]
                except:
                    pass
            else:
                try:
                    hash_map[tag] = genomic_location_hash[genomic_location]
                except:
                    pass
                
    return hash_map

def rv_num_hash_sym_num_rv_GPL13237():
# Hashes Rv num to symbol in the series GPL13237
    hash_map = {}
            
    with open("RVintersection/GSE182749_family.xml/GPL13237-tbl-1.txt") as file:
        for line in file:
            columns = line.split('\t')
            tag = columns[0]
            
            rv = re.findall("Rv[a-zA-Z0-9]+", line)
            
            if rv:
                try:
                    
                    hash_map[tag] = rv[0]
                except:
                    pass
    return hash_map


def mtub_to_rv_GSE13246_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_mtub_num_GPL7477()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE13246_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split(' ')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE23498_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_qmt_num_GPL5773()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE23498_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[2]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE35484_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_mtub_num_GPL7477()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE35484_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE37973_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_cust_num_GPL15565()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE37973_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE40917_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_cust_num_GPL14824()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE40917_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE43466_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_cust_num_GPL14824()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE43466_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE43749_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_sym_num_GPL4519()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE43749_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE46432_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_tbid_num_GPL17082()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE46432_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names
    

# this series has repeats in the GPL file and lots of genes are missing
def mtub_to_rv_GSE53843_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_jiehe_num_GPL18119()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE53843_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE55979_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_cust_num_GPL14824()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE55979_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE60376_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_sym_num_GPL18055()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE60376_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE67445_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_mtub_num_GPL7477()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE67445_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE90858_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_tb_num_rv_GPL22695()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE90858_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE96639_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_mtub_num_GPL7477()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE96639_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def mtub_to_rv_GSE182749_extraction(intersection_set):
    #goes through sample files and splits into key and value
    hash_map = rv_num_hash_sym_num_rv_GPL13237()
    sample_names = []
    data_matrix = {}
    directory = 'RVintersection/GSE182749_family.xml'
    for filename in os.listdir(directory):
        if filename.startswith('log2') or filename.startswith('GSM'):
            sample_name = filename[filename.find("GSM"):filename.find("-")]
            sample_names.append(sample_name)
            sample_path = f"{directory}/{filename}"
            with open(sample_path) as file:
                for line in file:
                    columns = line.strip().split('\t')
                    key = columns[0]
                    data = columns[1]
                    # matches key to rv tag + add with data value to matrix
                    if key in hash_map:
                        rv_tag = hash_map[key]
                        if rv_tag in intersection_set:
                            if rv_tag not in data_matrix:
                                data_matrix[rv_tag] = {}
                            try:
                                data_matrix[rv_tag][sample_name].append(data)
                            except:
                                data_matrix[rv_tag][sample_name] = []
                                data_matrix[rv_tag][sample_name].append(data)
                                
                for rv_tag in data_matrix:
                    for num in range(len(data_matrix[rv_tag][sample_name])):
                        data_matrix[rv_tag][sample_name][num] = float(data_matrix[rv_tag][sample_name][num])

                    data_matrix[rv_tag][sample_name] = st.median(data_matrix[rv_tag][sample_name])
        
    # print(data_matrix)    
    # for rv_tag, values in data_matrix.items():
        # print(f"{rv_tag}: {values}")
    # print(f"Sample names: {sample_names}")
    # print()
    return data_matrix, sample_names

def create_2d_array(data_matrix, sample_names):
    """Create a 2D array with sample names, data, and Rv tags and print it."""
    sample_names = sorted(sample_names)
    array = [['RvTag'] + sample_names]
    for rv_tag, values in data_matrix.items():
        row = [rv_tag] + [str(values.get(sample, '')) for sample in sample_names]
        array.append(row)
    print("2D Array:")
    for row in array:
        print('\t'.join(row))
    print()
    dfarray = pd.DataFrame(array[1:], columns=array[0])
    return dfarray

if __name__ == "__main__":
    """
    Key
    GSE13246 - GPL7477 - 8 samples
    GSE23498 = GPL5773 - 18 samples
    GSE25484 - GPL7477 - 6 samples
    GSE40917 - GPL14824 - 9 samples
    GSE43466 - GPL14824 - 202 samples
    GSE46432 - GPL17082 - 9 samples
    GSE53843 - GPL18119 - 12 samples
    GSE55979 - GPL14824 - 12 samples
    GSE67445 - GPL7477 - 8 samples
    GSE90858 - GPL22696 - 4 samples
    GSE96639 - GPL7477 - 6 samples
    GSE182749 - GPL13237 - 8 samples
    
    GPL7477 - 28 total samples
    GPL5773 - 18 total samples
    GPL15565 - 5 total samples
    GPL14824 - 223 total samples
    GPL17082 - 9 total samples
    GPL18119 - 12 total samples
    GPL22696 - 4 total samples
    GPL13237 - 8 total samples
    Total: 302 samples
    """

    # Compiled the Rv genes from all GPL tables
    rv_GPL7477 = rv_extraction_complete("RVintersection/GSE13246_family.xml/GPL7477-tbl-1.txt") # 3847
    rv_GPL5773 = rv_extraction_complete("RVintersection/GSE23498_family.xml/GPL5773-tbl-1.txt") # 3838
    rv_GPL14824 = rv_extraction_complete("RVintersection/GSE40917_family.xml/GPL14824-tbl-1.txt") # 3931
    rv_GPL17082 = rv_extraction_complete("RVintersection/GSE46432_family.xml/GPL17082-tbl-1.txt") # 3948
    rv_GPL18119 = rv_extraction_complete("RVintersection/GSE53843_family.xml/GPL18119-tbl-1.txt") # 3863
    rv_GPL22695 = rv_extraction_complete("RVintersection/GSE90858_family.xml/GPL22695-tbl-1.txt") # 3986
    rv_GPL13237 = rv_extraction_complete("RVintersection/GSE182749_family.xml/GPL13237-tbl-1.txt") # 3940
    
    """
    # Find the number of duplicates and unique duplicates for each GPL table
    duplicates_rv_GPL7477 = duplicates("RVintersection/GSE13246_family.xml/GPL7477-tbl-1.txt") # 93 duplicates, 57 unique duplicates
    duplicates_rv_GPL5773 = duplicates("RVintersection/GSE23498_family.xml/GPL5773-tbl-1.txt") # 291 duplicates, 221 unique duplicates
    duplicates_rv_GPL14824 = duplicates("RVintersection/GSE40917_family.xml/GPL14824-tbl-1.txt") # 12 duplicates, 6 unique duplicates
    duplicates_rv_GPL17082 = duplicates("RVintersection/GSE46432_family.xml/GPL17082-tbl-1.txt") # 8849 duplicates, 2202 unique duplicates
    duplicates_rv_GPL18119 = duplicates("RVintersection/GSE53843_family.xml/GPL18119-tbl-1.txt") # 11039 duplicates, 3851 unique duplicates
    duplicates_rv_GPL22695 = duplicates("RVintersection/GSE90858_family.xml/GPL22695-tbl-1.txt") # 3224 duplicates, 2881 unique duplicates
    duplicates_rv_GPL13237 = duplicates("RVintersection/GSE182749_family.xml/GPL13237-tbl-1.txt") #  14 duplicates, 13 unique duplicates
    print([duplicates_rv_GPL7477, duplicates_rv_GPL5773, duplicates_rv_GPL14824, duplicates_rv_GPL17082, duplicates_rv_GPL18119, duplicates_rv_GPL22695, duplicates_rv_GPL13237])
    """
    
    intersection_set = intersection(rv_GPL7477, rv_GPL5773, rv_GPL14824, rv_GPL17082, rv_GPL18119, rv_GPL22695, rv_GPL13237) # 3652 # removed GPL15565 - 3719 # removed GPL18055 - 3730 # remove both GPL15565 and GPL18055 - 3805
    union_set = union(rv_GPL7477, rv_GPL5773, rv_GPL14824, rv_GPL17082, rv_GPL18119, rv_GPL22695, rv_GPL13237) # 4007 # removed GPL15565 - 4007 # removed GPL18055 - 4007 # remove both GPL15565 and GPL18055 - 4007
    difference = union_set - intersection_set
    
    # print(len(difference))
    # GPL_counts = np.array([len(rv_GPL7477), len(rv_GPL5773), len(rv_GPL15565), len(rv_GPL14824), len(rv_GPL17082), len(rv_GPL18119), len(rv_GPL18055), len(rv_GPL22695), len(rv_GPL13237)])
    
    # gpl_dict = {"rv_GPL7477": rv_GPL7477, "rv_GPL5773": rv_GPL5773, "rv_GPL15565": rv_GPL15565, "rv_GPL14824": rv_GPL14824, "rv_GPL17082": rv_GPL17082, "rv_GPL18119": rv_GPL18119, "rv_GPL18055": rv_GPL18055, "rv_GPL22695": rv_GPL22695, "rv_GPL13237": rv_GPL13237}
    # print(missing_gene_in_gpls(difference, gpl_dict))



    series = [
        #list of series functions
        mtub_to_rv_GSE13246_extraction(intersection_set),
        mtub_to_rv_GSE23498_extraction(intersection_set),
        mtub_to_rv_GSE35484_extraction(intersection_set),
        mtub_to_rv_GSE40917_extraction(intersection_set),
        mtub_to_rv_GSE46432_extraction(intersection_set),
        mtub_to_rv_GSE53843_extraction(intersection_set),
        mtub_to_rv_GSE55979_extraction(intersection_set),
        mtub_to_rv_GSE67445_extraction(intersection_set),
        mtub_to_rv_GSE90858_extraction(intersection_set),
        mtub_to_rv_GSE182749_extraction(intersection_set),
        
    ]
    # combines data from each series extraction
    combined_metadata = {}
    all_sample_names = set()
    for serie in series:
        metadata, sample_names = serie
        all_sample_names.update(sample_names)
        for rv_tag, values in metadata.items():
            if rv_tag not in combined_metadata:
                combined_metadata[rv_tag] = {}
            combined_metadata[rv_tag].update(values)
    
    # Prints the 2D array
    dfarray = create_2d_array(combined_metadata, list(all_sample_names))
    # saves to .tsv file
    output_filename = "samples_genes_matrix.tsv"
    dfarray.to_csv(output_filename, sep='\t', index=False)
    print(f"Data saved to {output_filename}")
    
    
    """
    # Compiled all the Rv gene counts from all samples based on their respective GPL table
    counts = np.array(total_sample_gene_count())
    print(counts)
    
    # Plots a histogram of the number of Rv genes in each set
    plt.hist(counts, bins=10)
    plt.xlim(3000,4000)
    plt.xlabel('Rv Genes')
    plt.ylabel('Counts')
    plt.title('Number of Rv genes in each sample')
    plt.grid()
    plt.show()
    """
    
    
   
