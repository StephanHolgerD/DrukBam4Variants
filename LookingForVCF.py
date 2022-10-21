import pysam
import numpy as np
import pandas as pd
import sys






class ReadFiles():
    def __init__(self,ReferenceFasta,BamFile,VcfFile) -> None:
        self.ReferenceFasta=ReferenceFasta
        self.BamFile = BamFile
        self.VcfFile = VcfFile
        a = self.IterVariants()
        print(a)


    def IterVariants(self):
        with pysam.VariantFile(self.VcfFile) as f:
            for variant in f:
                print(variant)
                ReadData = self.GetReads(variant.contig,variant.pos)
                print(ReadData)


    def GetReads(self,chrom,pos):
        ResultArray = []
        with pysam.AlignmentFile(self.BamFile) as f:
            for alignment in f.fetch(chrom,pos):
                ResultArray.append(alignment.reference_start)
        return ResultArray
    


if __name__ == "__main__":
    
    b = sys.argv[1]
    v = sys.argv[2]
    
    xx = ReadFiles('foo',b,v)