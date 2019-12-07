#!/usr/bin/env python

import re
import argparse


#output file
#output_file = "/Users/GioiTran/Documents/shell/Bi624/Deduper/output.sam"
output_file = "/home/ntran2/bgmp/Bi624/Deduper/output.sam"


def get_arguements():
    ```This function handles all argparse agruements i.e, specifying input sam files, UMI file, and pair-ended files (optional)```
    parser = argparse.ArgumentParser(description="flags are set up to specify the desired (sorted) input sam file, UMI file, and an optional flag for pair-ended reads ")
    parser.add_argument("-f", "--file", help="this flag will specify the desired sorted, input sam file. Please specify the absolute path", type =str, required=True)
    parser.add_argument("-p", "--paired", help="optional arg, designates file is paired end (not single-end)", type =str, required=False)
    parser.add_argument("-u", "--umi", help="optional arg, designates file containing the list of UMIs (unset if randomers instead of UMIs)", type =str, required=True)
    return parser.parse_args()

args=get_arguements()
f=args.file
p=args.paired
u=args.umi


def get_UMI(file):
    ```This function will take in a file specifying the list of known UMIs. This function will then capture the UMIs and return all known UMIs into a list ```
    #reading in file
    with open(u, "rt") as fh:
        #creating empty list to hold UMIs
        known_UMI = []
        #iterate through each line
        for line in fh:
            #appending UMIs to the known_UMI list each time, remove the \n
            known_UMI.append(line.rstrip('\n'))
    #return known_UMI list when function is called
    return(known_UMI)

def parsing_QNAME(qname):
    ```This function takes in the QNAME, which is the first element in the SAM format, as a parameter. It will split the qname by the ":" and return the last element (which contains the )```
    #.split() will automatically store sections into a list, delimited by each ":"
    qname_list = qname.split(":")
    #print(qname_list)
    #return the last element in the qname_list list
    return qname_list[-1] #get the last element in the list (where barcode is)

def forward(CIGAR_str, POS):
    ```This function will take in the cigar string and original position taken from the SAM read as its parameters if the reads are forward. It will return the adjusted position by subtracting the original POS if the leftmost cigar string contains a "S"```
    #using the regex, findall method to store the number and letter of the cigar string into a tuple
    cigar_list = re.findall('(\d+)([MNS])', CIGAR_str )

    #if there's an "S" in the first element of the cigar tuple, the adjusted cigar is: (original position - number in the leftmost "S" in cigar string)
    if "S" in cigar_list[0]:
        fw_softclip = (cigar_list[0])

        adj_POS = (POS) - int(fw_softclip[0])
        return(adj_POS)
    # if there's no "S" in the first element of the cigar tuple, the original POS is uneffected
    if "S" not in cigar_list[0]:
        adj_POS = POS
        return(adj_POS)




def reverse(CIGAR_str, POS):
    ```This function will take in the cigar string and original position taken from the SAM read as its parameters if the reads are reverse. It will return the adjusted position by adding numbers from N's, M's, D's and S's to the original POS if the rightmost cigar string contains a "S. Return the adjusted POS"```
    #create a tuple by spliting the store the number and character of the cigar
    cigar_list = re.findall('(\d+)([MNSD])', CIGAR_str )
    #print(cigar_list)
    #if there's an "S" in the leftmost position, delete that tuple from our list
    if "S" in cigar_list[0]:
        del cigar_list[0]
        #print(cigar_list[0])

        #add N's, D's, M's and "S" to the original cigar string
        for char in cigar_list:
            if "N" or "D" or "M" or "S" in char:
                POS = POS + int(char[0])

    #if there isn't a "S" in the leftmost position, add N's, D's, M's and "S" to the original cigar string
    else:
        for char in cigar_list:
            if "N" or "D" or "M" or "S" in char:
                POS = POS + int(char[0])
                #print(char[0])
    return(POS)



#open the output file, once
file_out=open(f+"_deduped", "a")

#open the input file, store as a file handle
with open(f, "rt") as file_handle:

    #create empty cets for the forward and reverse records
    record_set_fw = set()
    record_set_rv = set()

    #call on the get_UMI to return the list of known UMIs
    umi_list = get_UMI(u)
    previous_chromosome = 0

    #iterate through each line of file without allocating to everything to memory
    for line in file_handle:

        #if there's an "@" in header, write it to file
        if "@" in line:
            file_out.write(line)
            #pass

        #split record into a list, deimited by tabs
        else:
            ref_rec = line.split()

            #extract chromosome/RNAME
            next_chromosome = ref_rec[2]
            #group sets by chromosome

            #clear the sets with each chromosome, method of not slowing down the algorithm
            if next_chromosome != previous_chromosome:
                previous_chromosome = next_chromosome
                record_set_fw.clear()
                record_set_rv.clear()


            #execute this condiiton if the UMI is known and if the bitwise flag is forward
            if parsing_QNAME(ref_rec[0]) in umi_list and ((int(ref_rec[1]) & 16) != 16):
                #call on the forward function()
                adj_5prime_pos_fw = forward(ref_rec[5], int(ref_rec[3]))
                #add qname, chromosome, and adjusted 5' position into a tuple
                record_tuple_fw=  (parsing_QNAME(ref_rec[0]) ,ref_rec[2],adj_5prime_pos_fw)


                #write to file if that record tuple is unique to the set
                if record_tuple_fw not in record_set_fw:
                    record_set_fw.add(record_tuple_fw)

                    file_out.write(line)


            #account for adjusted position if bitwise flag is reverse
            if parsing_QNAME(ref_rec[0]) in umi_list and ((int(ref_rec[1]) & 16) == 16):

                #call on the reverse function
                adj_5prime_pos_rv = reverse(ref_rec[5], int(ref_rec[3]))

                record_tuple_rv = (parsing_QNAME(ref_rec[0]), ref_rec[2], adj_5prime_pos_rv)

                #write to file if the record tuple is unique to the set
                if record_tuple_rv not in record_set_rv:
                    record_set_rv.add(record_tuple_rv)
                    file_out.write(line)



file_out.close()
