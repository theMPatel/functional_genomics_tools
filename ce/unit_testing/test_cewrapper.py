# -*- coding: utf-8 -*-

import sys
import os
import shutil

flPath = os.path.realpath(__file__)
myPath = os.path.split(flPath)[0]
basePath = os.path.split(os.path.split(myPath)[0])[0]
if basePath not in sys.path: sys.path.append(basePath)

import importlib

from mock import patch

def test_allelecalling():
    cewrapper = importlib.import_module('cewrapper.cewrapper')

    testargs = [
        "prog", 
        '--localdir', r'c:\users\hannes\temp\test_allelecalling\local',
        '--resultsdir', r'c:\users\hannes\temp\test_allelecalling',
        '--shareddir', r'c:\users\hannes\temp\test_allelecalling\shared',
        '--toolsdir', r'C:\tfs\CppCode\MAIN\Standalones\CEExecutables\Python\extern',
        '--tempdir', r'\users\hannes\temp\test_allelecalling\shared\temp',         
        '--algorithm', 'allelecalling',
        '--server', 'http://dev-ce',
        '--project', 'test_hannes',
        '--login', '601-597-463',
        '--password', 'test123',
        '--organism', 'CAMP',
        '--alleleCallsFlName', r'C:\Users\hannes\temp\test_allelecalling\ERR1046033_70.0_loci.xml',
        '--minSequenceIdentity', '70',
        '--nGaps', '1000',
        '--requireStartStop', 'True',
        '--avoidRepeats', 'True',
        '--sampleId', 'test',
        '--labId', 'hannes',
        '--isCacheOwner', 'True',
        '--submitCachedAlleles', 'False'
    ]
    
    with patch.object(sys, 'argv', testargs):
        retval = cewrapper.main()
        print retval

   
def test_genotyping():

    cewrapper = importlib.import_module('cewrapper.cewrapper')

    testargs = [
        "prog", 
        '--localdir', r'c:\users\hannes\temp\test_genotyping\local',
        '--resultsdir', r'c:\users\hannes\temp\test_genotyping',
        '--shareddir', r'C:\tfs\CppCode\MAIN\standalones\ceexecutables\python\genotyping\shared',
        '--toolsdir', r'C:\tfs\CppCode\MAIN\Standalones\CEExecutables\Python\extern',
        '--tempdir', r'C:\tfs\CppCode\MAIN\standalones\ceexecutables\python\genotyping\shared\temp',         
        '--configFile', r'[SHAREDDIR]\config.json',
        '--algorithm', 'genotyping',      
        '--organism', 'ecoli',        
        '--query', r'c:\users\hannes\temp\test_genotyping\Escherichia_coli_O157H7_Sakai.fasta',
        '--readFile', r'C:\Users\hannes\temp\test_genotyping\samples_ecoli\XH43223C_S6_L001_R1_001.fastq.gz',
        '--readFile', r'C:\Users\hannes\temp\test_genotyping\samples_ecoli\XH43223C_S6_L001_R2_001.fastq.gz',
        #'--readFile', r'c:\users\hannes\temp\test_genotyping\reads_16S.fastq',        
        #'--readFile', r'c:\users\hannes\temp\test_genotyping\reads_16S_2.fastq',
        #'--genotyper', 'plasmids',
        #'--genotyper', 'resistance', 
        #'--genotyper', 'ecoli.serotype', 
        #'--genotyper', 'ecoli.pathotype',
        #'--genotyper', 'virulence',
        #'--genotyper', 'insilicopcr',
        '--genotyper', 'ecoli.stxtype',
        #'--genotyper', 'mtbc.spoligotyping',
        #'--genotyper', 'mtbc.speciesdetermination',
        #'--genotyper', 'mtbc.resistancemapper',
        #'--readFile', r'c:\users\hannes\temp\test_genotyping\samples_mtbc\reads\BCG-GTAGAGGA-CTAAGCCT_L002_R1_001.fastq.gz',
        #'--readFile', r'c:\users\hannes\temp\test_genotyping\samples_mtbc\reads\BCG-GTAGAGGA-CTAAGCCT_L002_R2_001.fastq.gz',
    ]
    
    #from tools.tools import ReadSeqFile, WriteSeqFile
    #seqs=  ReadSeqFile( r'c:\users\hannes\temp\test_genotyping\SCP03-11.fasta', rename=True)
    #WriteSeqFile(seqs, r'c:\users\hannes\temp\test_genotyping\SCP03-11.fasta')
    with patch.object(sys, 'argv', testargs):
        retval = cewrapper.main()
        print retval

        
def test_genotyping_folder_2():

    cewrapper = importlib.import_module('cewrapper.cewrapper')

    import shutil
    from tools.tools import EnsureExists
    
    cewrapper = importlib.import_module('cewrapper.cewrapper')

    resultsDir = r'c:\users\hannes\temp\test_genotyping'
    rawResultsDir = os.path.join(resultsDir, 'results', 'raw')
    outputDir = r'c:\Users\hannes\temp\test_genotyping\mtbc_results'
    samples = {}
    
    EnsureExists(resultsDir)
    
    baseDir = r'D:\Anne.Holmes_nhs.uk'
    
    def GetKey(root, fl):
        return fl.split('_')[0] 
    def GetExt(fl):
        if 'fastq.gz' in fl: return 'fastq.gz'
        return '?'
        
    #for root, dirs, files in os.walk(baseDir):
    #   for file in files:
    for file in os.listdir(baseDir):
            root=baseDir
            key = GetKey(root, file)        
            ext = GetExt(file)
            
            if (ext!='fastq.gz'):
                continue
            
            if key not in samples: samples[key] = []
            samples[key].append(os.path.join(root, file))
     
    for key, files in samples.iteritems():
        print key, " " , len(files)

    for key, files in samples.iteritems():
        print "Analyzing sample {} ...".format(key)

        testargs = [
            "prog", 
            '--localdir', r'c:\users\hannes\temp\test_genotyping\local',
            '--resultsdir', r'c:\users\hannes\temp\test_genotyping',
            '--shareddir', r'C:\tfs\CppCode\MAIN\standalones\ceexecutables\python\genotyping\shared',
            '--toolsdir', r'C:\tfs\CppCode\MAIN\Standalones\CEExecutables\Python\extern',
            '--tempdir', r'C:\tfs\CppCode\MAIN\standalones\ceexecutables\python\genotyping\shared\temp',         
            '--configFile', r'[SHAREDDIR]\config.json',
            '--algorithm', 'genotyping',      
            '--organism', 'ecoli',        
            '--query', r'c:\users\hannes\temp\test_genotyping\Escherichia_coli_O157H7_Sakai.fasta',
            '--readFile', r'C:\Users\hannes\temp\test_genotyping\samples_ecoli\XH43223C_S6_L001_R1_001.fastq.gz',
            '--readFile', r'C:\Users\hannes\temp\test_genotyping\samples_ecoli\XH43223C_S6_L001_R2_001.fastq.gz',
            #'--readFile', r'c:\users\hannes\temp\test_genotyping\reads_16S.fastq',        
            #'--readFile', r'c:\users\hannes\temp\test_genotyping\reads_16S_2.fastq',
            #'--genotyper', 'plasmids',
            #'--genotyper', 'resistance', 
            #'--genotyper', 'ecoli.serotype', 
            #'--genotyper', 'ecoli.pathotype',
            #'--genotyper', 'virulence',
            #'--genotyper', 'insilicopcr',
            '--genotyper', 'ecoli.stxtype',
            #'--genotyper', 'mtbc.spoligotyping',
            #'--genotyper', 'mtbc.speciesdetermination',
            #'--genotyper', 'mtbc.resistancemapper',
            #'--readFile', r'c:\users\hannes\temp\test_genotyping\samples_mtbc\reads\BCG-GTAGAGGA-CTAAGCCT_L002_R1_001.fastq.gz',
            #'--readFile', r'c:\users\hannes\temp\test_genotyping\samples_mtbc\reads\BCG-GTAGAGGA-CTAAGCCT_L002_R2_001.fastq.gz',
        ]
        for fl in files: 
            target = os.path.join(r'c:\users\hannes\temp\test_genotyping\local', os.path.split(fl)[-1])
            shutil.copyfile(fl, target)
            testargs+= [ '--readFile', '"{}"'.format(target) ]
    
        #from tools.tools import ReadSeqFile, WriteSeqFile
        #seqs=  ReadSeqFile( r'c:\users\hannes\temp\test_genotyping\SCP03-11.fasta', rename=True)
        #WriteSeqFile(seqs, r'c:\users\hannes\temp\test_genotyping\SCP03-11.fasta')
        with patch.object(sys, 'argv', testargs):
            retval = cewrapper.main()
            print retval

        for fl in files: 
            target = os.path.join(r'd:\users\hannes\temp\test_genotyping\local', os.path.split(fl)[-1])
            os.remove(target)
            
        myOutputDir = os.path.join(outputDir, key)
        EnsureExists(myOutputDir)
        for fl in os.listdir(rawResultsDir):
            shutil.copyfile(
                os.path.join(rawResultsDir, fl),
                os.path.join(myOutputDir, fl),        
            )
        

   
def test_genotyping_folder():

    import shutil
    from tools.tools import EnsureExists
    
    cewrapper = importlib.import_module('cewrapper.cewrapper')

    resultsDir = r'c:\users\hannes\temp\test_genotyping'
    rawResultsDir = os.path.join(resultsDir, 'results', 'raw')
    outputDir = r'c:\Users\hannes\temp\test_genotyping\mtbc_results'
    samples = {}
    
    EnsureExists(resultsDir)
    
    baseDir = r'C:\Users\hannes\temp\test_genotyping\samples_mtbc'
    
    def GetKey(root, fl):
        return fl.split('_')[0] 
    def GetExt(fl):
        if 'fastq.gz' in fl: return 'fastq.gz'
        return '?'
        
    #for root, dirs, files in os.walk(baseDir):
    #   for file in files:
    for file in os.listdir(baseDir):
            root=baseDir
            key = GetKey(root, file)        
            ext = GetExt(file)
            
            if (ext!='fastq.gz'):
                continue
            
            if key not in samples: samples[key] = []
            samples[key].append(os.path.join(root, file))
     
    for key, files in samples.iteritems():
    
        print key, " " , len(files)
    
    for key, files in samples.iteritems():
    
        testargs = [     
            "prog", 
            '--localdir', r'c:\users\hannes\temp\test_genotyping\local',
            '--resultsdir', r'c:\users\hannes\temp\test_genotyping',
            '--shareddir', r'C:\tfs\CppCode\MAIN\standalones\ceexecutables\python\genotyping\shared',
            '--toolsdir', r'C:\tfs\CppCode\MAIN\Standalones\CEExecutables\Python\extern',
            '--tempdir', r'C:\tfs\CppCode\MAIN\standalones\ceexecutables\python\genotyping\shared\temp',         
            '--configFile', r'[SHAREDDIR]\config.json',
            '--algorithm', 'genotyping',      
            '--organism', 'mtbc',        
            '--query', r'',
            #'--readFile', r'c:\users\hannes\temp\test_genotyping\reads_16S.fastq',        
            #'--readFile', r'c:\users\hannes\temp\test_genotyping\reads_16S_2.fastq',
            #'--genotyper', 'plasmids',
            #'--genotyper', 'resistance', 
            #'--genotyper', 'ecoli.serotype', 
            #'--genotyper', 'ecoli.pathotype',
            #'--genotyper', 'virulence',
            #'--genotyper', 'insilicopcr',
            #'--genotyper', 'ecoli.stxtype',
            '--genotyper', 'mtbc.spoligotyping',
            '--genotyper', 'mtbc.speciesdetermination',
            '--genotyper', 'mtbc.resistancemapper',
        ]
        for fl in files: 
            #target = os.path.join(r'd:\users\hannes\temp\test_genotyping\local', os.path.split(fl)[-1])
            #shutil.copyfile(fl, target)
            testargs+= [ '--readFile', '"{}"'.format(fl) ]
        
        #from tools.tools import ReadSeqFile, WriteSeqFile
        #seqs=  ReadSeqFile( r'c:\users\hannes\temp\test_genotyping\SCP03-11.fasta', rename=True)
        #WriteSeqFile(seqs, r'c:\users\hannes\temp\test_genotyping\SCP03-11.fasta')A
        with patch.object(sys, 'argv', testargs):
            retval = cewrapper.main()
            print retval   
            
        myOutputDir = os.path.join(outputDir, key)
        EnsureExists(myOutputDir)
        for fl in os.listdir(rawResultsDir):
            shutil.copyfile(
                os.path.join(rawResultsDir, fl),
                os.path.join(myOutputDir, fl),        
            )
        
        
        
if __name__ =='__main__':     
    test_genotyping_folder_2()