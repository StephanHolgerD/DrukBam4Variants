import pysam


class WriteTSV_cls():
    def __init__(self,VariantCounterObject,VcfReaderObject,Outfile):
        self.ReferenceFasta=VcfReaderObject.ReferenceFasta
        self.BamFile = VcfReaderObject.BamFile
        self.VcfFile = VcfReaderObject.VcfFile
        self.VariantsInAlignments = {x[0]:x[1] for x in  VariantCounterObject.VariantsInAlignments}
        self.outfile = Outfile
        self.WriteVcfFile()
    
    
    def VarKey(self,variant):
        return f'{variant.chrom}-{variant.pos}-{variant.ref}'
    
    def WriteVcfFile(self):
        with pysam.VariantFile(self.VcfFile) as InVcf:
            with open(self.outfile,'w') as outfile:
                outfile.write(f'chrom\tpos\tref\tdepth\talts\n')
                for r in InVcf:
                    varkey = self.VarKey(r)
                    AF = self.VariantsInAlignments[varkey]
                    outfile.write(f'{r.chrom}\t{r.pos}\t{r.ref}\t{AF[1]}\t{AF[0]}\n')