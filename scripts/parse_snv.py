 #!/usr/bin/env python3
 
import sys

def parse(vcf_reader, vcf_writer, samplename):
    """
    Function to parse SNV vcf file (call or truth file)
    require the vcf path and the samplename if multisample vcf
    """

    variants_snv = []
    variants_indel = []
    #variants_sv = []
    sample = None
    nr_of_vars = 0
    nr_filtered = 0

    for record in vcf_reader:
        nr_of_vars += 1

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
                sample = str(record.samples[0])
                print("[INFO] Using only sample in the VCF: " + sample)

        chrom = record.CHROM.replace("CHR", "").replace("chr", "")
        pos = record.POS
        ref = record.REF
        alt = record.ALT[0]

        if record.var_type == "indel":
            len_ref = len(ref)
            len_alt = len(alt)
            if len_ref < 50 and len_alt < 50:
                vcf_writer.write_record(record)
                #print(record)
            if len_ref > 50 or len_alt > 50:
                nr_filtered += 1
                #print(record)
                #variants_sv.append([chrom, pos, ref, alt,len_ref,len_alt])
            #else:
            variants_indel.append([chrom, pos, ref, alt])
            #print("indel: {} {} {}".format(pos,ref,alt))
            #print(len(ref),len(alt))

        elif record.var_type == "snp":
            #print ("snp: {} {} {}".format(pos,ref,alt))
            variants_snv.append([chrom, pos, ref, alt])
            vcf_writer.write_record(record)
            #print(record)
        else:
            print(record)
            print("unknown: {} {} {} {}".format(record.var_type,pos,ref,alt))

    #return(variants_snv, variants_indel, variants_sv, nr_of_vars, nr_filtered)
    return(vcf_writer,variants_snv, variants_indel, nr_of_vars, nr_filtered)
