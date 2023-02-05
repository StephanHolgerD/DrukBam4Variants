import pysam

class DataCollector_cls():
    def __init__(self,OpenReferenceFasta,OpenBamFile):
        self.OpenReferenceFasta=OpenReferenceFasta
        self.OpenBamFile = OpenBamFile
    
    def GetReads(self,chrom,pos,share_d):
        fetchPos = pos -1
        ResultArray = []
        cov =0
        #with pysam.AlignmentFile(self.BamFile) as f:
        for alignment in self.OpenBamFile.fetch(chrom,fetchPos,pos,until_eof=False):

            if alignment.is_unmapped or alignment.query_alignment_sequence ==None:
                continue
            
            if cov >=200:
                return ResultArray
                
            cov = cov + 1
            
            if (alignment.query_name,alignment.is_read1) in share_d:
                try:
                    ResultArray.append(share_d[(alignment.query_name,alignment.is_read1)])
                    continue
                except KeyError:
                    pass
            
            parsedCigarString = self.CigChunker(alignment.cigarstring)
            r=(alignment.reference_name,
                                    alignment.reference_start,
                                    alignment.reference_end,
                                    alignment.query_alignment_sequence,
                                    parsedCigarString,
                                    alignment.query_name,
                                    alignment.is_read1)
            
            rr=self.ChunkFastaForAlignment(r)
            if len(list(share_d.keys())) > 1000:
                k = list(share_d.keys())[0]
                try:
                    share_d.pop(k)
                    gc.collect()
                except KeyError:
                    pass
            share_d[(alignment.query_name,alignment.is_read1)]=rr
            
            ResultArray.append(rr)
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
        
        Alignment_cigarstring=[x for x in cigarstring if x !='S' and x !='H' and x!='D']

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
            xx=[x for _,x in enumerate(query_alignment_sequence) if _+iposCounter in i]
            if len(xx)>0:
                
                insertion_sequence="".join(xx)
            else:
                insertion_sequence=""
            insertions.append((i,insertion_sequence))
            query_alignment_sequence="".join([x for _,x in enumerate(query_alignment_sequence) if _+iposCounter not in i])
            
            iposCounter=iposCounter+len(i)

        Alignment_cigarstring=[x for x in cigarstring if x !='S' and x !='H' and x!='I']
        for p,N in enumerate(Alignment_cigarstring):
            if N=='D' or N=='N':
                qs1=query_alignment_sequence[:p]
                qs2=query_alignment_sequence[p:]
                query_alignment_sequence=qs1+'-'+qs2
        return (GenomicPositions,fastaChunk,query_alignment_sequence,Alignment_cigarstring,InsertionPositions,insertions,query_name,is_read1)