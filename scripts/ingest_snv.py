 #!/usr/bin/env python3
 
import sys
import argparse
import vcf
import pandas as pd
import parse_snv

def main():
    """
    Main function to parse SNV & Indel vcf files and export an indel dataframe & filtered vcf without indels bigger thant 50bp.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="Path to SNV vcf file")
    parser.add_argument("-samplename", help="Sample name if multisampleVCF")
    parser.add_argument("-outputfile", help="If dataframe  & filtered vcf should be stored")

    args = parser.parse_args()

    print("[INFO] ### Start ingesting files")
    print("[INFO] file: " + str(args.vcf))
    variants = []

    # Account for: .vcf, .vcf.gz
    if args.vcf[-3:] == ".gz":
        vcf_reader = vcf.Reader(filename=args.vcf)
    elif args.vcf[-4:] == ".vcf":
        vcf_reader = vcf.Reader(open(args.vcf, 'r'))

    else:
        sys.exit("[ERROR] Do not recognize type of input vcf " + args.vcf + ". Exiting.")

    # Input file type is VCF. Start parsing.
    vcf_writer = vcf.Writer(open(args.outputfile + ".filtered.vcf", 'w'), vcf_reader)

    vcf_writer, variants_snv, variants_indel, nr_of_vars, nr_filtered = parse_snv.parse(vcf_reader, vcf_writer, args.samplename)

    vcf_writer.flush()

    # Create dataframe from all variant data
    if variants_snv or variants_indel:
        # if variants_snv:
        #     columns = ["chrom", "pos", "ref", "alt"]
        #     data_snv = pd.DataFrame(variants_snv, columns=columns)
        #     print("\n[INFO] Dataframe for SNV variants")
        #     print(data_snv)
        #     if args.outputfile:
        #         data_snv.to_csv(args.outputfile + "_snv.csv",index=False)
        # else:
        #     print("\n[INFO] Dataframe for SNV variants")
        #     print("Empty Dataframe ")
        if variants_indel:
            columns = ["chrom", "pos", "ref", "alt"]
            data_indel = pd.DataFrame(variants_indel, columns=columns)
            print("\n[INFO] Dataframe for indel variants")
            print(data_indel)

            if args.outputfile:
                data_indel.to_csv(args.outputfile + "_indel.csv", index=False)
        else:
            print("\n[INFO] Dataframe for indel variants")
            print("Empty Dataframe ")
        # if variants_sv:
        #     columns = ["chrom", "pos", "ref", "alt", "len_ref", "len_alt"]
        #     data_sv = pd.DataFrame(variants_sv, columns=columns)
        #     print("\n[INFO] Dataframe for structural variants")
        #     print(data_sv)
        #
        #     if args.outputfile:
        #         data_sv.to_csv(args.outputfile + "_sv.csv",index=False)
        # else:
        #     print("\n[INFO] Dataframe for structural variants")
        #     print("Empty Dataframe ")
    else:
        sys.exit("No data parsed, please check your input data.")

    # Save dataframe to file
    #if args.outputfile:
        #print("saving file")
        #data_snv.to_csv(args.outputfile + "_snv.csv")
        #data_indel.to_csv(args.outputfile + "_indel.csv")
        #data_sv.to_csv(args.outputfile + "_sv.csv")
    #print("\n")
    print("\nTotal number of variants processed: " + str(nr_of_vars))
    print("Filtered variants of total: " + str(nr_filtered))

if __name__ == "__main__":
    main()
