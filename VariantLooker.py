from re import S
import pysam

class VariantLooker_cls():
    def __init__(self,Variant,AlignmnentData):
        self.Variant = Variant
        self.GenomicPositions=AlignmnentData[0]
        self.fastaChunk=AlignmnentData[1]
        self.query_alignment_sequence=AlignmnentData[2]
        self.Alignment_cigarstring=AlignmnentData[3]
        self.InsertionPositions=AlignmnentData[4]
    
    def EvaluateAlignment(self):
        for AltAlelle in self.Variant.alts:
            if len(AltAlelle) == len(self.Variant.ref):
                SupportingAlignments = self.SameLengthSubstitution(AltAlelle)
                return SupportingAlignments
            
            
    def SameLengthSubstitution(self,AltAlelle):
        c = 0
        pos = self.Variant.pos -1
        GenomicStart = pos
        GenomicEnd = pos + len(self.Variant.ref)
        FastaRef = [self.fastaChunk[n] for n,x in enumerate(self.GenomicPositions) if x >= GenomicStart and x < GenomicEnd]
        ReadSeq = [self.query_alignment_sequence[n] for n,x in enumerate(self.GenomicPositions) if x >= GenomicStart and x < GenomicEnd]
        
        if len(FastaRef) > 0:
            print(FastaRef,AltAlelle,ReadSeq)