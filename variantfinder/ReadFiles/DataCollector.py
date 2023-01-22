import pysam

class DataCollector_cls():
    def __init__(self,OpenReferenceFasta,OpenBamFile):
        self.OpenReferenceFasta=OpenReferenceFasta
        self.OpenBamFile = OpenBamFile
    
    def GetReads(self,chrom,pos):
        fetchPos = pos -1
        ResultArray_temp = []
        ResultArray = []
        #with pysam.AlignmentFile(self.BamFile) as f:
        for alignment in self.OpenBamFile.fetch(chrom,fetchPos,pos,until_eof=False):
            if alignment.is_unmapped:
                continue
            parsedCigarString = self.CigChunker(alignment.cigarstring)
            ResultArray_temp.append((alignment.reference_name,
                                    alignment.reference_start,
                                    alignment.reference_end,
                                    alignment.query_alignment_sequence,
                                    parsedCigarString,
                                    alignment.query_name,
                                    alignment.is_read1))
        for r in ResultArray_temp:
            ResultArray.append(self.ChunkFastaForAlignment(r))
        
        return ResultArray
    
    def CigChunker(self,cigarstring):
        cigL=[]
        def parseCig(cig):
            for p,l in enumerate(cig):
                if l.isalpha():
                    liste=([cig[p]]*int(cig[:p]))
                    cig=cig[p+1:]
                    return cig,liste
        while cigarstring:
            cigL=cigL+parseCig(cigarstring)[1]
            cigarstring=parseCig(cigarstring)[0]
        return cigL

    def listConsec(self,l):
        retL=[]
        c=0
        for e,x in enumerate(l):
            if x+1 not in l:
                retL.append(l[c:e+1])
                c=e+1
        return(retL)
    
    def ChunkFastaForAlignment(self,AlignmentTuple):
        chrom = AlignmentTuple[0]
        start = AlignmentTuple[1]
        end = AlignmentTuple[2]
        query_alignment_sequence=AlignmentTuple[3]
        cigarstring = AlignmentTuple[4]
        query_name = AlignmentTuple[5]
        is_read1 = AlignmentTuple[6]
        
        Alignment_cigarstring=[x for x in cigarstring if x !='S' and x !='H']
        GenomicPositions = list(range(start,end))
        #with pysam.FastaFile(self.ReferenceFasta) as fa:
        if chrom not in self.OpenReferenceFasta.references:
            return
        fastaChunk=str(self.OpenReferenceFasta.fetch(str(chrom),start,end)).upper()
        InsertionPositions=[x for x,y in enumerate(Alignment_cigarstring) if y=='I']
        InsertionPositions=self.listConsec(InsertionPositions)
        iposCounter=0
        insertions = []
        for i in InsertionPositions:
            insertion_sequence="".join([x for _,x in enumerate(query_alignment_sequence) if _+iposCounter in i])
            insertions.append((i,insertion_sequence))
            query_alignment_sequence="".join([x for _,x in enumerate(query_alignment_sequence) if _+iposCounter not in i])
            
            iposCounter=iposCounter+len(i)
        for p,N in enumerate(Alignment_cigarstring):
            if N=='D' or N=='N':
                qs1=query_alignment_sequence[:p]
                qs2=query_alignment_sequence[p:]
                query_alignment_sequence=qs1+'-'+qs2
        return (GenomicPositions,fastaChunk,query_alignment_sequence,Alignment_cigarstring,InsertionPositions,insertions,query_name,is_read1)