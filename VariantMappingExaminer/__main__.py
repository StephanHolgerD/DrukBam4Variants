from concurrent.futures import process
from VariantMappingExaminer.ReadFiles.VcfReader import ReadVcfFile
from VariantMappingExaminer.WriteFiles.VcfWriter import WriteVcf_cls

from VariantMappingExaminer.ParseAlignments.VariantCounter import VariantCounter_cls
import sys
from time import time

if __name__ == "__main__":
    
    v = sys.argv[1]
    f = sys.argv[2]
    bs = sys.argv[3:]
    for b in bs:
        tt= time()
        x = ReadVcfFile(f,b,v)
        VC_object = VariantCounter_cls(x,processes=6)
        Vcf_Writer_object = WriteVcf_cls(VC_object,x)
        print(tt-time())    