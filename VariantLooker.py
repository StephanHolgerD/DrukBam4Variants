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
        self.Insertions=AlignmnentData[5]
        
    
    def EvaluateAlignment(self):
        for AltAlelle in self.Variant.alts:
            if len(AltAlelle) == len(self.Variant.ref):
                SupportingAlignments = self.SameLengthSubstitution(AltAlelle)
                return SupportingAlignments
            
            if len(AltAlelle)>len(self.Variant.ref):
                SupportingAlignments = self.Insertion(AltAlelle)
                return SupportingAlignments
                
            if len(AltAlelle)<len(self.Variant.ref):
                SupportingAlignments = self.Deletion(AltAlelle)
                return SupportingAlignments
                
            
            
    def SameLengthSubstitution(self,AltAlelle):
        c = 0
        pos = self.Variant.pos -1
        GenomicStart = pos
        GenomicEnd = pos + len(self.Variant.ref) + 5
        FastaRef = [self.fastaChunk[n] for n,x in enumerate(self.GenomicPositions) if x >= GenomicStart and x < GenomicEnd]
        ReadSeq = [self.query_alignment_sequence[n] for n,x in enumerate(self.GenomicPositions) if x >= GenomicStart and x < GenomicEnd]
        
        if len(FastaRef) > 0:
            print(FastaRef,AltAlelle,ReadSeq)
            if ReadSeq==AltAlelle:
                return 1
        return 0
            
    def Insertion(self,AltAlelle):
        c = 0
        pos = self.Variant.pos -1
        GenomicStart = pos
        GenomicEnd = pos + len(self.Variant.ref)
        
        FastaRef = [self.fastaChunk[n] for n,x in enumerate(self.GenomicPositions) if x >= GenomicStart and x < GenomicEnd]
        ReadSeq = [self.query_alignment_sequence[n] for n,x in enumerate(self.GenomicPositions) if x >= GenomicStart and x < AlignmentEnd]
        for i in self.Insertions:
            start= i[0][0]
            seq = i[1]
            
            if seq == AltAlelle[1:] and self.Variant.pos == self.GenomicPositions[0]+start:
                
                print(self.Variant.pos,FastaRef,AltAlelle,self.Insertions,self.GenomicPositions[0]+start)
                return 1
        return 0
    
    
    def Deletion(self,AltAlelle):
        c = 0
        pos = self.Variant.pos -1
        GenomicStart = pos
        GenomicEnd = pos + len(self.Variant.ref) 
#        AlignmentEnd = pos + len(self.Variant.ref)

        FastaRef = [self.fastaChunk[n] for n,x in enumerate(self.GenomicPositions) if x >= GenomicStart and x < GenomicEnd]
        ReadSeq = [self.query_alignment_sequence[n] for n,x in enumerate(self.GenomicPositions) if x >= GenomicStart and x < GenomicEnd]
        seq = ''.join([x for x in ReadSeq if x != '-'])
        if seq == AltAlelle:
            print(FastaRef,ReadSeq,AltAlelle,self.Variant.ref)
            return 1
        return 0
        #if len(FastaRef) > 0:
         #   print(FastaRef,AltAlelle,ReadSeq)