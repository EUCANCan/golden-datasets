import sys
import argparse
import vcf
import pandas as pd
import parse_sv_vcf
import parse_sv_tsv
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="Path to SV vcf file")
    parser.add_argument("-samplename", help="Sample name if multisampleVCF")
    parser.add_argument("-outputfile", help="If dataframe should be saved to CSV, give output file path here")

    args = parser.parse_args()

    print("[INFO] ### Start ingesting files")
    print("[INFO] file: " + str(args.vcf))
    variants = []

    # Account for: .vcf, .vcf.gz, .tsv
    if args.vcf[-3:] == ".gz":
        vcf_reader = vcf.Reader(filename=args.vcf)
    elif args.vcf[-4:] == ".vcf":
        vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    elif args.vcf[-4:] == ".tsv":
        vcf_reader = open(args.vcf, 'r')
    else:
        sys.exit("[ERROR] Do not recognize type of input vcf " + args.vcf + ". Exiting.")

    if hasattr(vcf_reader, 'metadata'):
        # Input file type is VCF. Start parsing.
        print("[INFO] File format: VCF file")
        variants, nr_of_vars, nr_filtered = parse_sv_vcf.parse(vcf_reader, args.samplename)

    else:
        # Input file type is TSV. Start parsing.
        print("[INFO] File format: TSV file")
        variants, nr_of_vars, nr_filtered = parse_sv_tsv.parse(vcf_reader)

    # Create dataframe from all variant data
    if variants:
        columns = ["start_chrom", "start", "end_chrom", "end", "ref", "alt" ,"length", "type"]
        data = pd.DataFrame(variants, columns=columns)
        data["length"] = pd.to_numeric(data["length"])
        #print(data.dtypes)
        

        #Make all values in length positive (Curie had all DEL values in negative, and BND within the same chr)
        for index, row in data.iterrows():
            if row['length'] ==None:
                pass
            elif row['length'] < 0:
                #print("AAA")
                data['length'] = data['length'].replace(row['length'], abs(int(row['length'])))  

   #elif isinstance(row['length'], (np.floating, float, str)):
        #Delete length values for BND, we can evaluate this later on
        
        data.loc[data['type'] == 'BND', 'length'] = ''
        print("[INFO] All lengths are considered in their absolute value")
        print("[INFO] Breakends are considered to have no length")
        
        # Check for duplicate entries, to prevent penalizing the nodes when they have duplicate entries
        
        data = data.astype(str)  #this line is for drop.duplicates() to work properly 
        
        # Check for duplicate entries, to prevent penalizing the nodes when they have duplicate entries
        dups = data[data.duplicated(keep=False)]
        if not dups.empty:
            print("[WARNING] " + str(dups.shape[0]) + " duplicates found. Only keeping the first line of each duplicate entry. List of duplicate entries: ")
            print(dups)
            data = data.drop_duplicates(keep='first').reset_index(drop=True)   #Only keeping first of the duplicate rows
        #print(data)
        
       
        '''
        if data.shape() != data_no_duplicates.shape():
            print("[WARNING] duplicates found. Only keeping the first line of each duplicate entry. List of duplicate entries: ")
            print("Number of duplicates found:" + str( int(data.shape()) - int(data_no_duplicates.shape())  ))
            data = data_no_duplicates
        
        
        if not dups.empty:
            print("[WARNING] " + str(dups.shape[0]) + " duplicates found. Only keeping the first line of each duplicate entry. List of duplicate entries: ")
            print(data[data.duplicated(keep=False)])
            data = data.drop_duplicates(keep='first').reset_index(drop=True)   #Only keeping first of the duplicate rows
        print(data)
        '''
        
       
        
        # Save dataframe to file
        if args.outputfile:
            data.to_csv(args.outputfile)

        return(data)

    else:
        sys.exit("No data parsed, please check your input data.")

    # Todo: what about fusions?

    print("Total number of variants processed: " + str(nr_of_vars))
    print("Filtered variants of total: " + str(nr_filtered))

if __name__ == "__main__":
    main()
