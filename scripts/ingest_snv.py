import sys
import argparse
import vcf
import pandas as pd
import parse_snv

def main():
    """
    Main function to parse SNV vcf file and export it as a dataframe in the terminal or as a standalone csv file
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="Path to SNV vcf file")
    parser.add_argument("-samplename", help="Sample name if multisampleVCF")
    parser.add_argument("-outputfile", help="If dataframe should be stored")

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
    variants, nr_of_vars, nr_filtered = parse_snv.parse(vcf_reader, args.samplename)

    # Create dataframe from all variant data
    if variants:
        columns = ["chrom", "pos", "ref", "alt"]
        data = pd.DataFrame(variants, columns=columns)
        print(data)
    else:
        sys.exit("No data parsed, please check your input data.")

    # Save dataframe to file
    if args.outputfile:
        data.to_csv(args.outputfile)

    print("Total number of variants processed: " + str(nr_of_vars))
    print("Filtered variants of total: " + str(nr_filtered))

if __name__ == "__main__":
    main()
