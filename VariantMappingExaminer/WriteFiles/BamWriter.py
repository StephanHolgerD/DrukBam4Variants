import pysam


class WriteBam_cls():
    def __init__(self,VariantCounterObject,VcfReaderObject,Outfile):
        self.OutPath = Outfile.split('.')[0]
        self.ReferenceFasta=VcfReaderObject.ReferenceFasta
        self.BamFile = VcfReaderObject.BamFile
        self.OpenBamfile = pysam.AlignmentFile(self.BamFile, "rb")
        
         
        self.VariantsInAlignments = VariantCounterObject.VariantsInAlignments
        self.WriteBamFile()
        self.OpenBamfile.close()
    
    def VarKey(self,variant):
        return f'{variant.chrom}-{variant.pos}-{variant.ref}'
    
    def WriteBamFile(self):
        for variant in self.VariantsInAlignments:
            chrom = variant[0].split('-')[0]
            pos = int(variant[0].split('-')[1])
            fetchPos = pos-1
            
            for allele in variant[2]:
                o = f'{self.OutPath}-{variant[0]}-AltAllele-{allele}.bam'
                oo = f'{self.OutPath}-{variant[0]}-NotAltAllele-{allele}.bam'
                with pysam.AlignmentFile(o,"wb",template=self.OpenBamfile) as OutBam:
                    with pysam.AlignmentFile(oo,"wb",template=self.OpenBamfile) as OutBam2:
                
                        for r in self.OpenBamfile.fetch(chrom,fetchPos,pos,until_eof=False):
                            testSet=set()
                            testSet.add((r.query_name,r.is_read1))
                            if len(testSet.intersection(variant[2][allele])) > 0:
                                OutBam.write(r)
                            else:
                                OutBam2.write(r)
        