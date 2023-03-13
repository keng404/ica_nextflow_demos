#!/usr/bin/env python3

import sys
import glob
import os
import os.path
import re
from decimal import *
getcontext().prec = 4
import argparse

##Makes a summary Excel file when given a run folder from PhoeNiX
##Usage: >python Phoenix_Summary_Line_06-10-22.py -n Sequence_Name -t Trimmed_QC_Data_File -x Taxa_file -r Ratio_File -m MLST_File -u mutation_file -q Quast_File -a AR_GAMMA_File -v Hypervirulence_GAMMA_File -k trimd_kraken -s synopsis_file -o Out_File
## Written by Rich Stanton (njr5@cdc.gov) and Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary line')
    parser.add_argument('-n', '--name', required=True, help='sequence name')
    parser.add_argument('-t', '--trimmed', required=False, help='QC data file for trimmed reads')
    parser.add_argument('-x', '--taxa', dest="taxa", required=False, help='Taxa file')
    parser.add_argument('-r', '--ratio', required=False, help='assembly ratio file')
    parser.add_argument('-m', '--mlst', required=False, help='MLST file')
    parser.add_argument('-u', '--mutation', dest="mutations", required=False, help='Mutation file from AMRFinder')
    parser.add_argument('-q', '--quast', required=False, help='QUAST file')
    parser.add_argument('-a', '--ar', required=False, help='AR GAMMA file')
    parser.add_argument('-p', '--pf', required=False, help='PF GAMMA file')
    parser.add_argument('-v', '--vir', required=False, help='hypervirulence GAMMA file')
    parser.add_argument('-k', '--kraken_trim', dest="trimd_kraken", required=False, help='trimd_summary.txt from kraken2')
    parser.add_argument('-s', '--stats', dest="stats", required=False, help='Pipeline Stats file synopsis file')
    parser.add_argument('-o', '--out', required=True, help='output file name')
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'
CEND = '\033[0m'

# def MLST_Info_Only(MLST_file):
#     """Pulls MLST info from *_Scheme.mlst file"""
#     f = open(MLST_file, 'r')
#     f.readline()
#     String1 = f.readline()
#     ST = String1.split()[2]
#     Scheme = String1.split()[1]
#     Out = ST + ' (' + Scheme + ')'
#     return Out

# def MLST_ST(MLST_file):
#     """Pulls MLST info from *_Scheme.mlst file"""
#     ST_list = []
#     with open(MLST_file, 'r') as f:
#         lines = f.readlines()
#         lines.pop(0)
#         for line in lines:
#             ST = "ST" + line.split()[2]
#             ST_list.append(ST)
#         if "ST-" in ST_list:
#             ST_list = [re.sub("ST-", "-", i) for i in ST_list]
#     return ST_list

def MLST_Scheme(MLST_file):
    """Pulls MLST info from *_Scheme.mlst file"""
    Scheme_list = [[],[],[],[],[]]
    with open(MLST_file, 'r') as f:
        lines = f.readlines()
        lines.pop(0)
        #print(len(lines), lines)
        for rawline in lines:
            line=rawline.strip()
            split_line = line.split("\t")
            #print("\n".join(split_line))
            source = split_line[1]
            date = split_line[2]
            DB_ID = split_line[3]
            Scheme = str(split_line[4])
            alleles = "-".join(split_line[5:])
            if DB_ID in Scheme_list[0]:
                print("In Scheme_list[0]")
                print(Scheme_list[0])
                for i in range(0,len(Scheme_list[0])):
                    if DB_ID == Scheme_list[0][i]:
                        print("Adding to", Scheme_list[0][i], i)
                        if Scheme != "-" and "Novel" not in Scheme: #if Scheme != "-" and Scheme != "Novel_allele" and Scheme != "Novel_profile":
                            Scheme_list[1][i].append("ST"+str(Scheme))
                        else:
                            Scheme_list[1][i].append(Scheme)
                        Scheme_list[2][i].append(alleles)
                        Scheme_list[3][i].append(source)
                        Scheme_list[4][i].append(date)
                        print(Scheme_list)
            else:
                print("NOT in Scheme_list[0]")
                print(Scheme_list[0], Scheme_list[1], Scheme_list[2], Scheme_list[3], Scheme_list[4])
                Scheme_list[0].append(DB_ID)
                if Scheme != "-" and Scheme != "Novel_allele" and Scheme != "Novel_profile":
                    Scheme_list[1].append(["ST"+Scheme])
                else:
                    Scheme_list[1].append([Scheme])
                Scheme_list[2].append([alleles])
                Scheme_list[3].append([source])
                Scheme_list[4].append([date])
                print(Scheme_list[0], Scheme_list[1], Scheme_list[2], Scheme_list[3], Scheme_list[4])
            print("Yay")
            for i in Scheme_list:
                for j in i:
                    print(j)
            #print("\n".join(Scheme_list))
    return Scheme_list

def Contig_Count(input_quast):
    Contigs = '0'
    f = open(input_quast, 'r')
    String1 = f.readline()
    while String1 != '':
        if ('# contigs (>= 0 bp)' in String1):
            Contigs = String1.split()[-1]
            break
        String1 = f.readline()
    return Contigs

def Genome_Size(input_quast):
    Size = '0'
    f = open(input_quast, 'r')
    String1 = f.readline()
    while String1 != '':
        if ('Total length (>= 0 bp)' in String1):
            Size = String1.split()[-1]
            break
        String1 = f.readline()
    return Size

def N50_Length(input_quast):
    N50 = '0'
    f = open(input_quast, 'r')
    String1 = f.readline()
    while String1 != '':
        if ('N50' in String1):
            N50 = String1.split()[-1]
            break
        String1 = f.readline()
    return N50

def GC_Content(input_quast):
    NGC = '0'
    f = open(input_quast, 'r')
    String1 = f.readline()
    while String1 != '':
        if ('GC' in String1):
            GC = String1.split()[-1]
            break
        String1 = f.readline()
    return GC

def Assembly_Ratio(ratio_file):
    f = open(ratio_file, 'r')
    String1 = f.readline()
    Ratio = 'Unknown'
    SD = 'Unknown'
    while String1 != '':
        if ('Isolate_St.Devs' in String1):
            SD = String1.split()[1]
        elif ('Ratio' in String1):
            Ratio = String1.split()[1]
        String1 = f.readline()
    f.close()
    Out = Ratio + ' (' + SD + ')'
    return Out

def Assembly_Ratio_Length(ratio_file):
    f = open(ratio_file, 'r')
    String1 = f.readline()
    Length = 0
    while String1 != '':
        if ('Actual_length' in String1):
            Length = String1.split()[1]
        String1 = f.readline()
    f.close()
    Out = int(Length)
    return Out

def Assembly_Ratio_Species(ratio_file):
    f = open(ratio_file, 'r')
    String1 = f.readline()
    Species = 'Unknown'
    while String1 != '':
        if ('Tax:' in String1):
            Species = String1.split()[1:]
            Species = ' '.join(Species)
        String1 = f.readline()
    f.close()
    return Species

def Trimmed_BP(trimmed_counts_file):
    f = open(trimmed_counts_file, 'r')
    String1 = f.readline()
    String1 = f.readline()
    BP = String1.split()[-2]
    BP = int(BP)
    return BP

def Trim_Coverage(trimmed_counts_file, ratio_file):
    Length = Assembly_Ratio_Length(ratio_file)
    Trimmed = Trimmed_BP(trimmed_counts_file)
    Coverage = str(Decimal(Trimmed) / Decimal(Length))
    return Coverage

def Bla_Genes(input_gamma):
    with open(input_gamma, 'r') as f:
        header=next(f) # just use to skip first line
        Bla = []
        String1 = f.readline()
        for line in f:
            Cat = line.split('\t')[0].split('__')[4] # Drug category
            Gene = line.split('\t')[0].split('__')[2] # Gene Name
            Percent_Length = float(line.split('\t')[11])*100
            Codon_Percent = float(line.split('\t')[9])*100
            if ('LACTAM' in Cat.upper()):
                # Minimum 90% length required to be included in report, otherwise not reported
                if int(round(Percent_Length)) >= 90:
                    # Minimum 98% identity required to be included in report, otherwise not reported
                    if int(round(Codon_Percent)) >= 98:
                        Bla.append(Gene)
    Bla.sort()
    return Bla

def Non_Bla_Genes(input_gamma):
    with open(input_gamma, 'r') as f:
        header=next(f) # just use to skip first line
        Non_Bla = []
        String1 = f.readline()
        for line in f:
            Cat = line.split('\t')[0].split('__')[4] # Drug category
            Gene = line.split('\t')[0].split('__')[2] # Gene Name
            Percent_Length = float(line.split('\t')[11])*100
            Codon_Percent = float(line.split('\t')[9])*100
            if ('LACTAM' in Cat.upper()) == False:
                # Minimum 90% length required to be included in report, otherwise not reported
                if int(round(Percent_Length)) >= 90:
                    # Minimum 98% identity required to be included in report, otherwise not reported
                    if int(round(Codon_Percent)) >= 98:
                        Non_Bla.append(Gene)
    Non_Bla.sort()
    return Non_Bla

def HV_Genes(input_gamma):
    with open(input_gamma, 'r') as f:
        header=next(f) # just use to skip first line
        HV = []
        String1 = f.readline()
        for line in f:
            Gene = line.split('\t')[0]
            HV.append(Gene)
    HV.sort()
    num_lines = sum(1 for line in open(input_gamma))
    if num_lines==1:
        HV=""
    return HV

def WT_kraken_stats(stats):
    with open(stats, 'r') as f:
        for line in f:
            if line.startswith("KRAKEN2_CLASSIFY_WEIGHTED"):
                scaffold_match = re.findall(r': .*? with', line)[0]
                scaffold_match = re.sub( ": SUCCESS  : | with", '', scaffold_match)
    return scaffold_match

def QC_Pass(stats):
    status = []
    reason = []
    warning_count = 0
    with open(stats, 'r') as f:
        for line in f:
            if ": WARNING  :" in line:
                warning_count = warning_count + 1
            if line.startswith("Auto Pass/FAIL"):
                line_status = line.split(":")[1]
                line_reason = line.split(":")[2]
                status.append(line_status.strip())
                reason.append(line_reason.strip())
                status_end = str(status[0])
                if status_end == "PASS":
                    reason_end = ""
                else:
                    reason_end = str(reason[0])
    warning_count = str(warning_count)
    return status_end, reason_end, warning_count

def Get_Kraken_reads(stats, trimd_kraken):
#    if stats is not None:
#        with open(stats, 'r') as f:
#            for line in f:
#                if line.startswith("KRAKEN2_CLASSIFY_READS"):
#                    read_match = re.findall(r': .*? with', line)[0]
#                    read_match = re.sub( ": SUCCESS  : | with", '', read_match)
#    else:
    with open(trimd_kraken, "r") as f:
        for line in f:
            if line.startswith("G:"):
                genus_match = line.split(": ")[1]
                genus_percent = line.split(": ")[1]
                genus_match = re.sub( "\d+|\n|\s|\.", '', genus_match)
                genus_percent = re.sub( "[a-zA-Z]*|\n|\s", '', genus_percent)
            if line.startswith("s:"):
                species_match = line.split(": ")[1]
                species_percent = line.split(": ")[1]
                species_match = re.sub( "\d+|\n|\.", '', species_match)
                species_percent = re.sub( "[a-zA-Z]*|\n|\s", '', species_percent)
        read_match = genus_match + "(" + genus_percent + "%)" + species_match + "(" + species_percent + "%)"
    return read_match

def Get_Taxa_Source(taxa_file):
    with open(taxa_file, 'r') as f:
        first_line = f.readline()
        fline=first_line.strip().split("\t")
        #taxa_source = re.findall(r'\(.*?\)', first_line)[0]
        #taxa_source = re.sub( "\(|\)", '', taxa_source)
        taxa_source=fline[0]
        percent_match=fline[1]
        if (taxa_source == "ANI_REFSEQ"):
            #percent_match = re.findall(r'-.*?%ID', first_line)[0]
            #percent_match = re.sub( "-|ID", '', percent_match)
            percent_match = percent_match + " ANI_match"
        if (taxa_source == "kraken2_trimmed"):
            #percent_match = re.findall(r'-.*?-', first_line)[0]
            #percent_match = re.sub( "-", '', percent_match)
            percent_match = percent_match + "% Reads_assigned"
        if (taxa_source == "kraken2_wtasmbld"):
            #percent_match = re.findall(r'-.*?-', first_line)[0]
            #percent_match = re.sub( "-", '', percent_match)
            percent_match = percent_match + "% Scaffolds_assigned"
    return taxa_source, percent_match

def Get_Mutations(amr_file):
    point_mutations_list = []
    with open(amr_file, 'r') as f:
        #tsv_file = csv.reader(amr_file, delimiter="\t")
        lines = f.readlines()[1:]
        for line in lines:
            if "POINT" in line:
                point_mutations = line.split("\t")[5]
                point_mutations_list.append(point_mutations)
        if len(point_mutations_list) == 0:
            point_mutations_list = ""
        else:
            point_mutations_list = ','.join(point_mutations_list)
    return point_mutations_list

def Get_Plasmids(pf_file):
    plasmid_marker_list = []
    with open(pf_file, 'r') as f:
        header=next(f) # just use to skip first line
        for line in f:
            Gene = line.split('\t')[0]
            Percent_Length = float(line.split('\t')[14])*100
            Match_Percent = float(line.split('\t')[13])*100
            # Minimum 60% length required to be included in report, otherwise not reported
            if int(round(Percent_Length)) >= 60:
                # Minimum 98% identity required to be included in report, otherwise not reported
                if int(round(Match_Percent)) >= 98:
                    plasmid_marker_list.append(Gene)
    if len(plasmid_marker_list) == 0:
        plasmid_marker_list = ""
    else:
        plasmid_marker_list = ','.join(plasmid_marker_list)
    return plasmid_marker_list

def Get_BUSCO_Gene_Count(stats):
    with open(stats, 'r') as f:
        matched_line = [line for line in f if "BUSCO" in line]
        split_list = matched_line[0].split(':')
        #lineage = re.sub( "BUSCO_", '', split_list[0])
        #percent = split_list[2].split(' ')[1]
        #ratio = split_list[2].split(' ')[9].rstrip('\n')
        lineage="_".join(split_list[0].split("_")[-2:]).strip()
        percent=str(split_list[2].split("%")[0].strip())+"%"
        ratio="("+str(split_list[2].split("(")[1].strip())
        busco_line = percent + ' ' + ratio
    busco_file = True
    return busco_line, lineage, busco_file

def Isolate_Line(Taxa, ID, trimmed_counts, ratio_file, MLST_file, quast_file, gamma_ar, gamma_hv, stats, trimd_kraken, amr_file, pf_file):
    try:
        plasmid_marker_list = Get_Plasmids(pf_file)
    except:
        plasmid_marker_list = 'Unknown'
    try:
        point_mutations_list = Get_Mutations(amr_file)
    except:
        point_mutations_list = 'Unknown'
    try:
        taxa_source, percent_match = Get_Taxa_Source(Taxa)
    except:
        taxa_source = 'Unknown'
        percent_match = 'Unknown'
    try:
        Coverage = Trim_Coverage(trimmed_counts, ratio_file)
    except:
        Coverage = 'Unknown'
    try:
        Genome_Length = Genome_Size(quast_file)
    except:
        Genome_Length = 'Unknown'
    try:
        Ratio = Assembly_Ratio(ratio_file)
    except:
        Ratio = 'Unknown'
    try:
        Contigs = Contig_Count(quast_file)
    except:
        Contigs = 'Unknown'
    try:
        busco_line, lineage, busco_file = Get_BUSCO_Gene_Count(stats)
    except:
        busco_file = None
        busco_line = 'Unknown'
        lineage = 'Unknown'
    try:
        GC = GC_Content(quast_file)
    except:
        GC = 'Unknown'
    try:
        Species = Assembly_Ratio_Species(ratio_file)
    except:
        Species = 'Unknown'
    # try:
    #     ST = MLST_ST(MLST_file)
    #     if len(ST) > 1:
    #         ST_1 = ST[0]
    #         ST_2 = ST[1]
    #     else:
    #         ST_1 = ST[0]
    #         ST_2 = "-"
    # except:
    #     ST_1 = 'Unknown'
    #     ST_2 = 'Unknown'
    try:
        Scheme = MLST_Scheme(MLST_file)
        if len(Scheme[0]) > 1:
            if Scheme[0][0].lower() < Scheme[0][1].lower():
                MLST_scheme_1 = Scheme[0][0]
                #print("1,1-before sort", Scheme[1][0])
                mlst_types_1=sorted(Scheme[1][0])[::-1]
                MLST_type_1 = ",".join(mlst_types_1)
                if Scheme[3][0] == "srst2":
                    MLST_type_1 += '^'
                #print("1,1-after sort", MLST_type_1)
                #MLST_alleles_1 = ",".join(Scheme[2][0])
                MLST_scheme_2 = Scheme[0][1]
                #print("2,1-before sort", Scheme[1][1])
                mlst_types_2=sorted(Scheme[1][1])[::-1]
                MLST_type_2 = ",".join(mlst_types_2)
                if Scheme[3][1] == "srst2":
                    MLST_type_2 += '^'
                #print("2,1-after sort", MLST_type_2)
                #MLST_alleles_2 = ",".join(Scheme[2][1])
            else:
                MLST_scheme_1 = Scheme[0][1]
                #print("1,2-before sort", Scheme[1][1])
                mlst_types_1=sorted(Scheme[1][1])[::-1]
                MLST_type_1 = ",".join(mlst_types_1)
                if Scheme[3][1] == "srst2":
                    MLST_type_1 += '^'
                #print("1,2-after sort", MLST_type_1)
                #MLST_alleles_1 = ",".join(Scheme[2][1])
                MLST_scheme_2 = Scheme[0][0]
                #print("2,2-before sort", Scheme[1][0])
                mlst_types_2=sorted(Scheme[1][0])[::-1]
                MLST_type_2 = ",".join(mlst_types_2)
                if Scheme[3][0] == "srst2":
                    MLST_type_2 += '^'
                #print("2,2-after sort", MLST_type_2)
                #MLST_alleles_2 = ",".join(Scheme[2][0])
        else:
            MLST_scheme_1 = Scheme[0][0]
            MLST_type_1 = ",".join(Scheme[1][0])
            #MLST_alleles_1 = ",".join(Scheme[2][0])
            MLST_scheme_2 = "-"
            MLST_type_2 = "-"
            #MLST_alleles_2 = "-"
    except:
        MLST_scheme_1 = 'Unknown'
        MLST_scheme_2 = 'Unknown'
        MLST_type_1 = 'Unknown'
        MLST_type_2 = 'Unknown'
        #MLST_alleles_1 = 'Unknown'
        #MLST_alleles_2 = 'Unknown'
    try:
        Bla = Bla_Genes(gamma_ar)
        Bla = ','.join(Bla)
    except:
        Bla = 'Unknown'
    try:
        Non_Bla = Non_Bla_Genes(gamma_ar)
        Non_Bla = ','.join(Non_Bla)
    except:
        Non_Bla = 'Unknown'
    try:
        HV = HV_Genes(gamma_hv)
        if HV=="":
            pass
        else:
            HV = ','.join(HV)
    except:
        HV = 'Unknown'
    try:
        scaffold_match = WT_kraken_stats(stats)
    except:
        scaffold_match = "Unknown"
    try:
        QC_Outcome, Reason, warning_count = QC_Pass(stats)
    except:
        QC_Outcome = 'Unknown'
        warning_count = 'Unknown'
        Reason = ""
    try:
        read_match = Get_Kraken_reads(stats, trimd_kraken)
    except:
        read_match = "Unknown"
    if busco_file is None:
        Line = ID + '\t' + QC_Outcome + '\t' + warning_count + '\t'  + Coverage + '\t' + Genome_Length + '\t' + Ratio + '\t' + Contigs + '\t' + GC + '\t' + Species + '\t' + percent_match + '\t' + taxa_source + '\t' + read_match + '\t' + scaffold_match + '\t' + MLST_scheme_1 + '\t' + MLST_type_1 + '\t' + MLST_scheme_2 + '\t' + MLST_type_2 + '\t' + Bla + '\t' + Non_Bla + '\t' + point_mutations_list + '\t' + HV + '\t' + plasmid_marker_list + '\t' + Reason
        busco = False
    elif busco_file is not None:
        Line = ID + '\t' + QC_Outcome + '\t' + warning_count + '\t'  + Coverage + '\t' + Genome_Length + '\t' + Ratio + '\t' + Contigs + '\t' + GC + '\t' + busco_line + '\t' + lineage + '\t' + Species + '\t' + percent_match + '\t' + taxa_source + '\t' + read_match + '\t' + scaffold_match + '\t' + MLST_scheme_1 + '\t' + MLST_type_1 + '\t' + MLST_scheme_2 + '\t' + MLST_type_2 + '\t' + Bla + '\t' + Non_Bla + '\t' + point_mutations_list + '\t' + HV + '\t' + plasmid_marker_list + '\t' + Reason
        busco = True
    return Line, busco

def Isolate_Line_File(Taxa, ID, trimmed_counts, ratio_file, MLST_file, quast_file, gamma_ar, gamma_hv, out_file, stats, trimd_kraken, mutations, pf_file):
    with open(out_file, 'w') as f:
        Line, busco = Isolate_Line(Taxa, ID, trimmed_counts, ratio_file, MLST_file, quast_file, gamma_ar, gamma_hv, stats, trimd_kraken, mutations, pf_file)
        if busco == True:
            f.write('ID\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tBUSCO\tBUSCO_DB\tSpecies\tTaxa_Confidence\tTaxa_Source\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason\n')
        else:
            f.write('ID\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tSpecies\tTaxa_Confidence\tTaxa_Source\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason\n')
        f.write(Line)

def main():
    args = parseArgs()
    # if the output file already exists remove it
    Isolate_Line_File(args.taxa, args.name, args.trimmed, args.ratio, args.mlst, args.quast, args.ar, args.vir, args.out, args.stats, args.trimd_kraken, args.mutations, args.pf)

if __name__ == '__main__':
    main()
