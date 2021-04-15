import sys
import pandas as pd

def parse(vcf_reader):
    variants = []
    nr_of_vars = 0
    nr_filtered = 0
    ref_version = "Unknown"
    header = False
    #min_nr_progs = 2    # Define the number of programs that need to call the variant for it to continue
    
    for line in vcf_reader:
        line = line.rstrip().split()
       
        if not header:
            # Try finding a reference genome file
            if "##reference" in line[0]:
                ref_version = "\t".join(line)
                
            # Parse in header line for BSC, Charite, Truth file tsv
            if line[0] in ["CHR1", "#chrom1", "new_id"]:
                header = line
                print("[INFO] Header of TSV file found.")
                #print(header)
    
            # If variants are reached without column names been found, stop execution
            if line[0] in list(map(str,list(range(1,23)))) \
                or line[0] in ["X", "Y"] \
                    or line[0].startswith("chr") \
                        or line[0].startswith("truthset"):
                if header == False:
                    sys.exit("[ERROR] Column names not found. Exiting.")
    
            continue
    
        nr_of_vars += 1
        
        # FILTERING (Option to separate this in a different function)
        # BSC support
        if 'CHR1' in header: 
            if line[header.index('FILTER')] != "PASS":
                #print("[DEBUG] Found filter: " + str(line[header.index('CUSTOM_FILTER')]))
                nr_filtered += 1
                continue
            # Maybe we could recover these discarded variants in another dataframe for checking purposes
            '''
            # Filter out variants that do not get called by min_nr_progs
            if int(line[header.index('NPROGS')]) < min_nr_progs:
                nr_filtered += 1
                continue
            '''
        # Charite support
        elif '#chrom1' in header: # Charite does not have FILTER column so it does not apply. We need to ask them.
            pass 
        # COLO829 truth file (no need to filter)
        elif "new_id" in header:
            pass
    
    
        # FIND VARIANT INFO (option to separate in a different function)
        # BSC TSV support
        if 'CHR1' in header:
            start_chrom = line[header.index('CHR1')]
            start = line[header.index('POS1')]
            end_chrom = line[header.index('CHR2')]
            end = line[header.index('POS2')]
            ref = None
            alt = None
            if start_chrom == end_chrom:
                length = int(end) - int(start)
            else:
                length = None #Only for SV involving different chr.. maybe use sth more specific.
    
        # Charite TSV support
        elif '#chrom1' in header:   
            start_chrom = line[header.index('#chrom1')]
            start = line[header.index('start1')]
            end_chrom = line[header.index('chrom2')]
            end = line[header.index('start2')]
            ref = None
            alt = None
            if start_chrom == end_chrom:
                length = int(end) - int(start)
            else:
                length = None #Only for SV involving different chr.. maybe use sth more specific.
    
        #Truth COLO829 support
        elif "new_id" in header: 
            start_chrom = line[header.index('chr1')]
            start = line[header.index('pos1')]
            end_chrom = line[header.index('chr2')]
            #For insertions we assume that the chrom is the same and the end pos is start +1
            if end_chrom == "<INS>":
                end_chrom = start_chrom
            end = line[header.index('pos2')]
            if end == ".":
                end = int(start) +1
            ref = line[header.index('ref')]
            alt = None
            length = int(line[header.index('size')])
            if length in [0, "0"]:
                length = None #Only for SV involving different chr.. maybe use sth more specific.
    
            
        # FIND SV TYPE
        # BSC TSV support
        if 'CHR1' in header:
            svtypes = []
            progs = line[header.index('PROGS')]
            for prog in progs.split(","):
                if "SVTYPE_" + prog in header:
                    svtypes.append(line[header.index("SVTYPE_" + prog)])
            svtypes = list(dict.fromkeys(svtypes))
            if 'deletion' in svtypes and 'DEL' in svtypes:
                svtypes.remove('deletion')
            if 'DUP' in svtypes and 'tandem-duplication' in svtypes:
                svtypes.remove('tandem-duplication')
            if 'inversion' in svtypes and 'INV' in svtypes:
                svtypes.remove('inversion')
            if len(svtypes) == 1:
                sv_type = svtypes[0]
            elif len(svtypes) == 0:
                sv_type = None
            else:
                if "BND" in svtypes:
                    svtypes.remove('BND')
                    if len(svtypes) == 1:
                        sv_type = svtypes[0]
                    else:
                        print("Unknown SV combination! Please check!")
                        print(svtypes)
                        sv_type = None
                else:
                    print("Unknown SV combination! Please check!")
                    print(svtypes)
                    sv_type = None
    
        # Charite support
        elif '#chrom1' in header:
            sv_type = line[header.index('svtype')]
        #Truth COLO829 support
        elif "new_id" in header:
            sv_type = line[header.index('type')]
            
        # Reformat types
        sv_type = sv_type.replace("deletion", "DEL")
        sv_type = sv_type.replace("inversion", "INV")
        sv_type = sv_type.replace("translocation", "TRA")
        sv_type = sv_type.replace("insertion", "INS")
        sv_type = sv_type.replace("duplication", "DUP")
    
        #variants["length"] = pd.to_numeric(variants["length"])
    
        #Gather variant info
        variants.append([start_chrom, start, end_chrom, end, ref, alt, length, sv_type])
    
    print("[INFO] Reference version: " + ref_version)

    return(variants, nr_of_vars, nr_filtered)