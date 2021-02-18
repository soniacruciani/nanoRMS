### this script is only necessary when epinano-rms successfully creates a temporary folder with *_freq files but fails to merge them into the final .csv file. 
### run it as: df_proc.py $yourtempfolder


#!/usr/bin/env python
import sys
import numpy as np
import numpy as np
#import df_proc
import dask
import dask.dataframe as dd
import pandas as pd
import os

def df_is_not_empty(df):
        '''
        input df is a df filtred on reference id
        if is is empty: next (df.iterrows()) does not work
        otherwise it returns a row of df
        '''
        try:
                next (df.iterrows())
                return True
        except:
                return False


def df_proc (df, outfn):
        '''
        input is a dataframe for either forward or reverse strand
        '''
        if not df_is_not_empty (df):
                print ("empty dataframe for {}".format(outfn), file=sys.stderr)
                return None
        outfh = open (outfn, 'w')
        header = "#Ref,pos,base,strand,cov,q_mean,q_median,q_std,mis,ins,del,ACGT_freq"
        print (header, file=outfh)
        gb = df.groupby(['#Ref','pos','base','strand']).agg({
               'cov':['sum'],
               'mis':['sum'],
               'ins':['sum'],
               'del':['sum'],
               'qual':['sum'],
                'bases':['sum']})
        for row in gb.itertuples():
                index = ",".join (map (str,row[0]))
                cov, mis, ins, _del, qual,bases = row[1:]
                mis = '%0.5f' % (mis/cov)
                ins = '%0.5f' % (ins/cov)
                _de = '%0.5f' % (_del/cov)
                q = np.array ([x for x in qual.split(':') if x ]).astype(int)
                qmn,qme,qst = '%0.5f' % np.mean(q), '%0.5f' % np.median(q), '%0.5f' % np.std(q)
                ACGT_freq = []
                acgt = ''
                try:
                        bases = np.array(  [ ele for ele in  bases.split(':') if ele] ).astype(int)
                        acgt = np.array([0,0,0,0])
                        for x in range (0, len(bases), 4):
                                acgt = acgt + np.array (bases[x:x+4])
                        outfh.write ("{},{},{},{},{},{},{},{},{}\n".format(index,cov,qmn,qme,qst,mis,ins,_de,":".join(map (str,acgt))))
                except:
                        print ('warning:',row,file=sys.stderr)
                        raise

def main():
        tmp_folder=sys.argv[1] # os.path.dirname (sys.argv[1])
        df = dd.read_csv ("{}/small_*freq".format(tmp_folder))
        df = df.compute()
        tmp_folder = os.path.split(os.path.realpath(tmp_folder))[1]
        out_prefix = tmp_folder.replace ('.tmp_splitted_','')
        df_fr = df[df['strand']=="+"]
        df_rev = df[df['strand'] == "-"]
        df_proc(df_fr,'{}.fw.txt'.format(out_prefix))
        df_proc(df_rev,'{}.re.txt'.format(out_prefix))

if __name__ == "__main__":
        main()
