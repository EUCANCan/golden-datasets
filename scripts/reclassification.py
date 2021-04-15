import sys
import argparse
import pandas as pd

#TO DO
#Make sure both indel and SV file are normalized in the same way beforehand
#Stablish proper formats (Columns), indel calling do not yields end chr and end position and sv type etc.
#Deal with real indel file, now it is tested with dummy files
#Make sure to re-evaluate the events with ref length == alt length
'''
# COLUMNS SNV: chrom,position,ref,alt
# COLUMNS SV: index,start_chrom,start,end_chrom,end,ref,alt,length,type

From INDEL/SNV to SV
# 1: end chrom = start chrom
# 2: position end = postion start + length
# 3: length = max(length(ref), length(alt))
# 4: type = del/ins/delins/other?

ref = AGC
alt = CTAATA
type = delins
length = ? > length(alt) - length(ref)

From SV to INDEL/SNV
1: chrom = start_chrom
2: position = start
3: ref = ref
4: alt = alt
'''


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("df_indels", help="Path to indels dataframe file")
    parser.add_argument("df_sv", help="Path to SV dataframe file")
    parser.add_argument("new_df_indels", help="Path to new indels dataframe file variants <50 bp")
    parser.add_argument("new_df_sv", help="Path to new sv dataframe file variants >=50 bp")   

    args = parser.parse_args()
   
    print("[INFO] ### Reclassificating variants <50bp as indels and variants >=50bp as structural variants")


    df_indels = pd.read_csv(args.df_indels, dtype={'chrom': str, 'pos':int, 'ref':str, 'alt':str, })
    df_sv = pd.read_csv(args.df_sv, index_col=0,dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': int, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})


    if len(df_indels.index)>0 and len(df_sv.index)>0 : #If files are not empy

        for index, row in df_indels.iterrows(): #Add required columns to go from indel format to SV format
            df_indels.at[index,'length'] = max(len(df_indels.at[index,'ref']), len(df_indels.at[index,'alt']))
            df_indels.at[index,'type'] = 'DEL' if len(df_indels.at[index,'ref']) > len(df_indels.at[index,'alt']) else 'INS'
        

        df_indels =  df_indels.astype( {'length' : int})
        df_sv = df_sv.astype({'length': int},errors='ignore')

    
    
        indels_longer_than_50 =  df_indels[df_indels['length']>=50]
        indels_smaller_than_50 = df_indels.append(indels_longer_than_50).drop_duplicates(keep = False)
        
        indels_longer_than_50 = indels_longer_than_50.rename(columns={'chrom':'start_chrom','pos':'start'})             
        indels_longer_than_50.insert(loc=2, column='end_chrom', value=indels_longer_than_50['start_chrom'])
        indels_longer_than_50.insert(loc=3, column='end', value=None  )
        indels_longer_than_50['end'] = indels_longer_than_50['start'] + indels_longer_than_50['length']
        indels_longer_than_50 = indels_longer_than_50.astype({"end": int})
        

        sv_longer_than_50_and_BND = df_sv.loc[(df_sv['length']>= 50) | (df_sv.length.isnull()), :]
       
        sv_smaller_than_50 = df_sv.append(sv_longer_than_50_and_BND).drop_duplicates(keep = False) 
       
        
        final_svs = sv_longer_than_50_and_BND.append(indels_longer_than_50, ignore_index=True)
        
        #Rename columns and drop some of them to go from SV format to indel dataframe format
        sv_smaller_than_50 = sv_smaller_than_50.rename(columns={"start_chrom": "chrom", "start": "pos"})
        sv_smaller_than_50 = sv_smaller_than_50.drop(['end_chrom','end'], axis=1)
        
        final_indels = indels_smaller_than_50.append(sv_smaller_than_50,ignore_index=True)
        final_indels = final_indels.drop(['length', 'type'], axis=1)
    
        
    
        if args.new_df_indels and args.new_df_sv:
                final_indels.to_csv(args.new_df_indels)
                final_svs.to_csv(args.new_df_sv)
    
        return final_indels, final_svs

    else:
        sys.exit("No data parsed, please check your input data.")



if __name__ == "__main__":
    main()
