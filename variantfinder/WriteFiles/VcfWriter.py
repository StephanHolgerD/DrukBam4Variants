import pysam


class WriteVcf_cls():
    def __init__(self,VariantCounterObject,VcfReaderObject,Outfile):
        self.ReferenceFasta=VcfReaderObject.ReferenceFasta
        self.BamFile = VcfReaderObject.BamFile
        self.InfoFieldShort = self.BamFile.split('/')[-1]
        self.InfoFieldLong = 'AF in ' + self.BamFile.split('/')[-1] + ' Alignments with Alt Allele Count / Total Allignemnts Count'
         
        self.VcfFile = VcfReaderObject.VcfFile
        self.VcfFileOut = Outfile
        self.VariantsInAlignments = {x[0]:x[1] for x in  VariantCounterObject.VariantsInAlignments}
        self.WriteVcfFile()
    
    
    def VarKey(self,variant):
        return f'{variant.chrom}-{variant.pos}-{variant.ref}'
    
    def WriteVcfFile(self):
        with pysam.VariantFile(self.VcfFile) as InVcf:
            InVcf.header.add_meta(key="INFO", items=[('ID',self.InfoFieldShort), 
                                                  ('Number','.'),
                                                  ('Type','String'), 
                                                  ('Description',self.InfoFieldLong)])
            with pysam.VariantFile(self.VcfFileOut,'w',header=InVcf.header) as OutVcf:
                for r in InVcf:
                    varkey = self.VarKey(r)
                    AF = self.VariantsInAlignments[varkey]
                    AF = str(AF)
                    r.info.__setitem__(self.InfoFieldShort, AF)
                    OutVcf.write(r)
                
        
        