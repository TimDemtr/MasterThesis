##### REASSES THIS IN THE FUTURE ->>>> DEF. NOT WORKING PROPERLY... SOIMETIMES ITS +1 -
# BUT AIN'T NOBODY GOT TIME FOR THAT RIGHT NOW

import os
from Bio import pairwise2
import multiprocessing

folder_path = "/Work2/Tim_D_Sequencing/2021_11_25_Epidermal_Immuno_EcoliPolyIC/08-piRNAClusters/"
bedfile_files = sorted([f"{folder_path}{x}" for x in os.listdir(folder_path) if x.endswith(".bed")])

## Formats explanation:
# getfasta:
# dd_Smes_g4_1:6-36(-)	tgtaaacgagaacaatgctcgagttaaatt

# bed2
# dd_Smes_g4_1	6	36	1	173	-	TGTAAACGAGAACAATGCTCGAGTTAAATT

# bed
# dd_Smes_g4_1	6	36	1825143	255	-

# desired .map output
# Chr1  238  TGTTACGGCTAGCTCAGTACGGC  23  TGTTACGGCTAGCTCAGTACGAA  2  +

## alternative .map ELAND3 output - seems to be valid as well...
# Chr1  LEFTMOSTSTART  -  NAME  SEQUENCE_Mapped  -  +
# -> so we can ditch the getfasta entirely... so unnecessary...

def multiprocess_wrapper(bedfile):
    total_size = os.path.getsize(bedfile)
    progress = 0
    stepstone = 0
    with open(bedfile, "r") as bed_h, \
        open(f"./bed2/{bedfile.split('.')[0]}.bed2", "r") as bed2_h, \
        open(f"{bedfile.split('/')[-1].replace('.bed', '.map')}", "w") as map_out:
        for n, bed_l in enumerate(bed_h):
            bed2_l = bed2_h.readline().split()
            chr = bed2_l[0]
            start = int(bed2_l[1])+1 # Checked ELAND3 output of "proTRAC" and all entries were off by +1... maybe that's why it's hanging at ~70% no matter what. 
            name = bed_l[-3]
            strand = bed_l[-1]
            mapped_seq = bed2_l[-1].upper()
            write_l = "\t".join([str(item) for item in [chr, start, ref_seq, name, mapped_seq, mismatches, strand]])
            map_out.write(f'{write_l}\n')
            progress += len(str.encode(bed_l))
            if (progress/total_size) > stepstone:
                with open(f"{bedfile.split('/')[-1].replace('.bed', '.log')}", "a+") as progress_h:
                    progress_h.write(f"{int((progress/total_size)*100)}%\n")
                stepstone += 0.01

Processes = []
for bedfile in bedfile_files:
    new_thread = multiprocessing.Process(target=multiprocess_wrapper, args=(bedfile, ))
    new_thread.start()
    Processes.append(new_thread)

for process in Processes:
    process.join()


