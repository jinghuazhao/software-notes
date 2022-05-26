#!/usr/bin/env python

########################################################
# plink2metasoft.py                                    
#   Convert Plink .assoc files to Metasoft input file  
#   Free license -- you are free to use it in any ways 
#   Buhm Han (2012)                                    
########################################################

import sys, subprocess, os

comple = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

# PROCESS ARGUMENTS
if len(sys.argv) < 3:
    print_and_exit()
out=sys.argv[1]
files=sys.argv[2:]

# READ FILES
studies=[]
for f in files:
    study={}
    fin=open(f)
    colnames=fin.next().split()
    for line in fin:
        snp={}
        for (x,y) in zip(colnames,line.split()):
            snp[x]=y
        rsid=snp['SNP']
        study[rsid]=snp
    fin.close()
    studies.append(study)
    
# UNION OF SNPS
allsnps={}
for study in studies:
    for rsid, snp in study.iteritems():
        allsnps[rsid]=snp

# SORT SNPS
rsids=sorted(allsnps.keys(), 
             key=lambda x:(allsnps[x]['CHR'],int(allsnps[x]['BP'])))

# MERGE STUDIES
fout=open(out+'.meta.tmp','w') 
fmap=open(out+'.mmap','w')
flog=open(out+'.log','w')  
for rsid in rsids:
    output=rsid+'\t'
    pivot=-1
    pivotstudyindex=-1
    numstudy=0
    for study in studies:
        studyindex=studies.index(study)+1
        if rsid in study:
            snp=study[rsid]
            if 'OR' in snp:
                beta='log(%s)'%snp['OR']
            elif 'BETA' in snp:
                beta=snp['BETA']
            else:
                assert 0, 'OR or BETA must be in columns'
            if 'P' in snp:
                p=snp['P']
            else:
                assert 0, 'P must be in columns'
            stderr='abs(%s/qnorm(%s/2))'%(beta,p)
            if pivot == -1:
                pivot=snp # 1ST STUDY's SNP INFO IS PIVOT
                pivotstudyindex=studyindex
            else:
                # CHECK ALLELE TO PIVOT
                if 'A2' in pivot and 'A2' in snp:
                    if pivot['A1'] == snp['A1'] and \
                       pivot['A2'] == snp['A2']:
                        # GOOD
                        pass
                    elif pivot['A1'] == snp['A2'] and \
                         pivot['A2'] == snp['A1']:
                        # SIMPLE FLIP
                        beta='-(%s)'%beta
                    elif pivot['A1'] == comple[snp['A1']] and \
                         pivot['A2'] == comple[snp['A2']]:
                        # STRAND INCONSIS., BUT GOOD
                        flog.write('FLIP_STRAND %s in study %d\n'%(rsid,studyindex))
                    elif pivot['A1'] == comple[snp['A2']] and \
                         pivot['A2'] == comple[snp['A1']]:
                        # STRAND INCONSIS., SIMPLE FLIP
                        flog.write('FLIP_STRAND %s in study %d\n'%(rsid,studyindex))
                        beta='-(%s)'%beta
                    else:
                        flog.write('EXCLUDE %s due to allele inconsistency: A1:%s A2:%s in study %d but A1:%s A2:%s in study %d\n'
                                   %(rsid, pivot['A1'], pivot['A2'], pivotstudyindex,
                                     snp['A1'], snp['A2'], studyindex))
                else: 
                    if pivot['A1'] == snp['A1']:
                        # GOOD
                        pass
                    else:
                        flog.write('EXCLUDE %s due to allele inconsistency: A1:%s in study %d but A1:%s in study %d\n'
                                   %(rsid, pivot['A1'], pivotstudyindex,
                                     snp['A1'], studyindex))
                # CHECK CHR & BP TO PIVOT
                if pivot['CHR'] != snp['CHR']:
                    flog.write('WARNING %s has chr inconsistency\n'%rsid)
                if pivot['BP'] != snp['BP']:
                    flog.write('WARNING %s has basepair inconsistency\n'%rsid)
            output+=beta+' '+stderr+' '
            numstudy+=1
        else:
            output+='NA NA '
    if numstudy == 1:
        flog.write('EXCLUDE %s due to being in single study\n'%rsid)
    else:
        fout.write(output+'\n')
        fmap.write('%s\t%s\t%s\t%s\t%d\n'%
                   (rsid, pivot['CHR'], pivot['BP'], pivot['A1'], numstudy))
fout.close()
fmap.close()
flog.close()

# CALL R TO EVALUATE MATH EXPRESSION
subprocess.call(['R --vanilla --slave "--args '+out+'.meta.tmp '+out+'.meta" < plink2metasoft_subroutine.R'],shell=True)
subprocess.call(['rm '+out+'.meta.tmp'], shell=True)


