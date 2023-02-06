import pysam
import os
import json
class DataCollector_cls():
    def __init__(self,OpenReferenceFasta,OpenBamFile,tmpdir=None,mode='shortread'):
        self.OpenReferenceFasta=OpenReferenceFasta
        self.OpenBamFile = OpenBamFile
        self.mode = mode
        self.tmpdir = tmpdir
    
    def GetReads(self,chrom,pos):
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
            
            if self.mode=='longread':
                try:
                    if os.path.isfile(f'{self.tmpdir}/{alignment.query_name}.json'):
                        print(f'{self.tmpdir}/{alignment.query_name}.json')
                        with open(f'{self.tmpdir}/{alignment.query_name}.json') as json_file:
                            data = json.load(json_file)
                            rr =(data['GenomicPositions'],
                             data['fastaChunk'],
                             data['query_alignment_sequence'],
                             data['Alignment_cigarstring'],
                             data['InsertionPositions'],
                             data['insertions'],
                             data['query_name'],
                             data['is_read1'])
                            ResultArray.append(rr)
                            continue
                except:
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
            jsondumpDict = {'GenomicPositions':rr[0],
                'fastaChunk':rr[1],
                'query_alignment_sequence':rr[2],
                'Alignment_cigarstring':rr[3],
                'InsertionPositions':rr[4],
                'insertions':rr[5],
                'query_name':rr[6],
                'is_read1':rr[7]}
            if self.mode=='longread':
                
                with open(f'{self.tmpdir}/{alignment.query_name}.json','w') as outfile:
                    json.dump(jsondumpDict, outfile)
            
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