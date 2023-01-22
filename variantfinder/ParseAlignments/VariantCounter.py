from unittest import result
from variantfinder.ParseAlignments.VariantLooker import VariantLooker_cls
from variantfinder.ReadFiles.DataCollector import DataCollector_cls
from multiprocessing import Pool
import numpy as np 
import pysam
from tqdm import tqdm
class VariantCounter_cls():
    def __init__(self,VcfReaderObject,processes=1):
        self.variants=VcfReaderObject.variants
        self.ReferenceFasta=VcfReaderObject.ReferenceFasta
        self.BamFile = VcfReaderObject.BamFile
        
        
        self.processes = processes
        self.VariantsInAlignments = self.CallMultiProcess()



    def CallMultiProcess(self):
        mp = [(k,v) for k,v in self.variants.items()]
        
        mp_split = np.array_split(mp, self.processes)
        print(f'....start looking for variant supporting aligments in {self.BamFile}')
        
        with Pool(processes=self.processes) as p:
            r = p.map(self.ParseMultiList,mp_split)
        rr = [xx for x in r for xx in x]
        return rr
    
    def ParseMultiList(self,l):
        results = []
        with pysam.AlignmentFile(self.BamFile) as OpenBamFile:
            with pysam.FastaFile(self.ReferenceFasta) as OpenReferenceFasta:
                for entry in tqdm(l):
                    results.append(self.CountVariant(entry[0],entry[1],OpenReferenceFasta,OpenBamFile))
        return results

    def CountVariant(self,v1,v2,OpenReferenceFasta,OpenBamFile):
        variant_key = v1
        variant = v2
        alignmentsVariant  = {}
        c_all=0
        c_var={}
        DataCollector = DataCollector_cls(OpenReferenceFasta,OpenBamFile)
        AlignmentData = DataCollector.GetReads(str(variant.contig),variant.pos)
        for alignment in AlignmentData:
            c_all=c_all+1
            x = VariantLooker_cls(variant,alignment)
            alt,res = x.EvaluateAlignment()
            if alt in c_var:
                c_var[alt]=c_var[alt]+res
            else:
                c_var[alt]=res
            if res == 0:
                continue
            
            if alt in alignmentsVariant:
                x = alignmentsVariant[alt]
                x.add((alignment[6],alignment[7]))
                alignmentsVariant[alt] = x
            else:
                alignmentsVariant[alt] = set()
                x = alignmentsVariant[alt]
                x.add((alignment[6],alignment[7]))
                alignmentsVariant[alt] = x
                
                
            #c_var = c_var + x.EvaluateAlignment()
        
        return (variant_key,(c_var,c_all),alignmentsVariant)
