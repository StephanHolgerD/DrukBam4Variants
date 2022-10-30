from symbol import comp_if
import pysam
import numpy as np
import pandas as pd
import sys
from DataCollector import DataCollector_cls
from VariantLooker import VariantLooker_cls




class ReadFiles():
    def __init__(self,ReferenceFasta,BamFile,VcfFile):
        self.ReferenceFasta=ReferenceFasta
        self.BamFile = BamFile
        self.VcfFile = VcfFile
        self.DataCollector = DataCollector_cls(self.ReferenceFasta,self.BamFile)
        a = self.IterVariants()


    def IterVariants(self):
        with pysam.VariantFile(self.VcfFile) as f:
            for variant in f:
                c_all=0
                c_var=0
                AlignmentData = self.DataCollector.GetReads(str(variant.contig),variant.pos)
                for alignment in AlignmentData:
                    c_all=c_all+1
                    x = VariantLooker_cls(variant,alignment)
                    c_var = c_var + x.EvaluateAlignment()
                print(c_var/c_all)


            

if __name__ == "__main__":
    
    b = sys.argv[1]
    v = sys.argv[2]
    f = sys.argv[3]
    xx = ReadFiles(f,b,v)