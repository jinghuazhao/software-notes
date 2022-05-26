


import os, sys
import numpy as np


args = sys.argv

if (len(args) != 8):
    print("Usage: python pmplot.py <infile> <outfile> <study_names_file> <study_order_files> rsid genename <pmplot_file(.pdf)>")
    sys.exit()


infile = args[1]
outfile = args[2]
studynamesfile = args[3]
studyorderfile = args[4]
rsid = args[5]
genename = args[6]
pmplot_file = args[7]

rsids = []
for line in open(outfile):
    rsids.append(line.split('\t')[0])
rsids = np.array(rsids)

outtab = []
for line in open(outfile):
    outtab.append(line)
outtab = np.array(outtab)
intab = np.loadtxt(infile, dtype='str')
studies = np.loadtxt(studynamesfile, dtype='str')
total_studies = len(studies)

rs_index_list = np.where(rsids == rsid)[0]
if (len(rs_index_list) == 0):
    print("rsid: " + rsid + " not found!")
    print("exiting ...")
    exit()

rs_index = np.where(rsids == rsid)[0][0]
def get_valid_index(intab, rs_index):
    intab = intab[rs_index-1,:] 
    return(np.where(intab[range(1,len(intab), 2)] != 'NA')[0])


def cleanup(outname):
    cmd = "rm tmpoutput.txt"
    os.system(cmd)
    cmd = "rm tmpinput.txt"
    os.system(cmd)
    cmd = "rm tmpgenenames.txt"
    os.system(cmd)
    cmd = "rm tmpstudies.txt"
    os.system(cmd)
    cmd = "rm tmpstudyorder.txt"
    os.system(cmd)
    cmd = "rm " + outname + ".Rout"
    os.system(cmd)


def cleanup_err():
    cmd = "rm tmpoutput.txt"
    os.system(cmd)
    cmd = "rm tmpinput.txt"
    os.system(cmd)
    cmd = "rm tmpgenenames.txt"
    os.system(cmd)
    cmd = "rm tmpstudies.txt"
    os.system(cmd)
    cmd = "rm tmpstudyorder.txt"
    os.system(cmd)



def const_outtab(outtab, rs_index, total_studies, sels):
    outtab = np.array(outtab[rs_index].split('\t'))
    outtab = outtab[:-1]
    header = outtab[0:16]
    pvals = outtab[16:(16+total_studies)]
    mvals = outtab[(16+total_studies):]
    newout = np.hstack((header, pvals[sels], mvals[sels]))
    return(newout)


def const_intab(intab, rs_index):
    intab = intab[rs_index-1,:] 
    return(intab[np.where(intab != 'NA')[0]])
    


if (len(intab.shape) == 1):
    intab.shape = (1, len(intab))

sels = get_valid_index(intab, rs_index)
intab = const_intab(intab, rs_index)
outtab = const_outtab(outtab, rs_index, total_studies, sels)
nstudy = int(outtab[1])
height = 4 + 0.23*nstudy


np.savetxt("tmpoutput.txt", outtab, fmt='%s')
np.savetxt("tmpinput.txt", intab, fmt='%s')
np.savetxt("tmpgenenames.txt", np.array([genename]), fmt='%s')
np.savetxt("tmpstudies.txt", studies[sels], fmt='%s')
np.savetxt("tmpstudyorder.txt", range(1, (len(sels)+1)), fmt='%s')
outname = str(np.random.randint(100000))


cmd = "R CMD BATCH --no-save --no-restore '--args tmpoutput.txt tmpinput.txt tmpgenenames.txt tmpstudies.txt tmpstudyorder.txt " + pmplot_file + " " + str(height) + "' forestpmplot.R "+outname+ ".Rout"
os.system(cmd)




fp = open(outname+".Rout")
lines = fp.readlines()
if ("proc.time" in lines[-3]):
    cleanup(outname)
    print(pmplot_file + " is successfully generated.")
else:
    cleanup_err()
    print("Problem occurred, while generating " + pmplot_file)
    ##os.system("cat outname.Rout")
    ## added by yurang.park
    os.system("cat "+str(outname)+".Rout")






