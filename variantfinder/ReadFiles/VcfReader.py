import pysam
import sys
from variantfinder.ReadFiles.DataCollector import DataCollector_cls


class ReadVcfFile():
    def __init__(self,ReferenceFasta,BamFile,VcfFile):
        self.ReferenceFasta=ReferenceFasta
        self.BamFile = BamFile
        self.VcfFile = VcfFile
        self.DataCollector = DataCollector_cls(self.ReferenceFasta,self.BamFile)
        self.ReadVariants()
    
    def ReadVariants(self):
        variants = {}
        with pysam.VariantFile(self.VcfFile) as f:
            for variant in f:
                x=Variant_cls(variant)
                variants[self.VarKey(variant)]=x
        self.variants = variants
        
    
    
    def VarKey(self,variant):
        return f'{variant.chrom}-{variant.pos}-{variant.ref}'
    
    



class Variant_cls():
    def __init__(self,variant):
        self.pos=variant.pos
        self.contig=variant.contig
        self.alts=variant.alts
        self.ref=variant.ref
        self.ref=variant.ref
        
        
