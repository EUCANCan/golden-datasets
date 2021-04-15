import sys

def parse(vcf_reader, samplename):
    variants = []
    sample = None
    nr_of_vars = 0
    nr_filtered = 0

    # Get reference genome version
    if 'reference' in vcf_reader.metadata:
        ref_version = vcf_reader.metadata['reference']
    else:
        ref_version = "Unknown"
    print("[INFO] Reference version: " + ref_version)

    # Get VCF format version
    if 'fileformat' in vcf_reader.metadata:
        vcf_format = vcf_reader.metadata['fileformat']
    else:
        vcf_format = "Unknown"
    print("[INFO] File format: " + vcf_format)

    # Start reading in records and add to dataframe
    for record in vcf_reader:
        nr_of_vars += 1
        # Only use 'PASS' calls
        if record.FILTER:

            # If we want to keep vars with certain filters, list here
            if record.FILTER != ['INFERRED']:
                #print('[DEBUG] Found filter: ' + str(record.FILTER))
                nr_filtered += 1
                continue

        # Select only the relevant sample and skip this step if performed once
        if sample == None:
            if len(record.samples) > 1:

                # Need to filter, so need optional argument
                if samplename == None:
                    sys.exit("[ERROR] Multisample VCF found and no samplename given. Please input the tumor samplename with "
                             "-samplename to make filtering possible.")
                else:
                    sample = samplename
                    print("[INFO] Filtering for the following sample in the VFC: " + sample)
            else:
                sample = record.samples[0]
                print("[INFO] Using only sample in the VCF: " + sample)

            # First get index of sample
            index_sample = None
            for ind, call in enumerate(record.samples):
                if call.sample == sample:
                    index_sample = ind
            if not index_sample:
                sys.exit("[ERROR] Cannot find sample in the VCF file. Exiting.")

        start_chrom = record.CHROM.replace("CHR", "").replace("chr", "")
        start = record.POS

        # Check if variant is SV
        if record.is_sv:
            if len(record.ALT) > 1:
                sys.exit("[DEBUG] MORE THAN ONE ALT ALLELE!")

            # Get end coordinate from file if registered in INFO field
            if 'END' in record.INFO:
                end = record.sv_end
                end_chrom = start_chrom
            else:
                # Parse end from ALT
                end_pos = str(record.ALT[0]).replace("CHR", "").replace("chr", "").replace("A", "").replace("G", "").replace("T", "").replace("C", "").replace("N", "").replace("]", "").replace("[", "")
                if "." in end_pos:
                    '''
                    # SV is probably BND?, possible alts: ".str", "str.", "N.", ".N"
                    print("[DEBUG] IS SV BND?")
                    print(record)
                    print(record.samples)
                    print(record.INFO)
                    '''
                    end_pos = start
                    end_chrom = start_chrom
                else:
                    #print(record)
                    end_chrom, end = end_pos.split(":")

            if 'SVLEN' in record.INFO:
                length = record.INFO['SVLEN']
                if isinstance(length, list):
                    # Length is list
                    length = length[0]
            
            elif 'LEFT_SVINSSEQ' in record.INFO:
                length = len(str(record.INFO['LEFT_SVINSSEQ'])) + len(str(record.INFO['RIGHT_SVINSSEQ']))                                    
            else:
                if end_chrom == start_chrom:
                    length = int(end) - int(start)
                else:
                    # SV starts and ends on different chromosomes
                    length = None
            ref = record.REF
            alt = record.ALT[0]

            sv_type = None
            if 'EVENTTYPE' in record.INFO: # Support Hartwig
                sv_type = record.INFO['EVENTTYPE']
            elif 'SVTYPE' in record.INFO: # Support Curie
                sv_type = record.INFO['SVTYPE']
            else:
                sv_type = record.var_subtype


        elif record.is_snv:
            sys.exit("[DEBUG] SNV in SV file.")


        # Convert into table with: start_chrom, start, end_chrom, end, ref, alt, length, sv_type, genotype
        #print([start_chrom, start, end_chrom, end, ref, alt, length, sv_type])
        variants.append([start_chrom, start, end_chrom, end, ref, alt, length, sv_type])

    return(variants, nr_of_vars, nr_filtered)