from typing import Counter
from VariantMappingExaminer.ParseAlignments.VariantLooker import VariantLooker_cls
from VariantMappingExaminer.ReadFiles.DataCollector import DataCollector_cls
from multiprocessing import Pool, process


class VariantCounter_cls():
    def __init__(self,VcfReaderObject,processes=1):
        self.variants=VcfReaderObject.variants
        self.ReferenceFasta=VcfReaderObject.ReferenceFasta
        self.BamFile = VcfReaderObject.BamFile
        
        
        self.processes = processes
        self.VariantsInAlignments = self.CallMultiProcess()



    def CallMultiProcess(self):
        mp = [(k,v) for k,v in self.variants.items()]
        with Pool(processes=self.processes) as p:
            r = p.starmap(self.CountVariant,mp)
        return r
            

    def CountVariant(self,v1,v2):
        variant_key = v1
        variant = v2
        c_all=0
        c_var=0
        DataCollector = DataCollector_cls(self.ReferenceFasta,self.BamFile)
        AlignmentData = DataCollector.GetReads(str(variant.contig),variant.pos)
        for alignment in AlignmentData:
            c_all=c_all+1
            x = VariantLooker_cls(variant,alignment)
            c_var = c_var + x.EvaluateAlignment()
        print(variant_key,(c_var,c_all))
        
        return (variant_key,(c_var,c_all))
