import pysam
import sys
from variantfinder.ReadFiles.DataCollector import DataCollector_cls


class ReadVcfFile():
    def __init__(self,ReferenceFasta,BamFile,VcfFile,tmpdir=None,mode='shortread'):
        self.ReferenceFasta=ReferenceFasta
        self.BamFile = BamFile
        self.VcfFile = VcfFile
        self.mode =mode
        self.tempdir = tmpdir
        if self.mode =='longread':
            self.tempdir = tempfile.mkdtemp()

        self.DataCollector = DataCollector_cls(self.ReferenceFasta,self.BamFile,tmpdir=self.tempdir,mode=self.mode)
        self.ReadVariants()
        if self.tempdir != None:
            shutil.rmtree(self.tempdir)
    
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
        
        
