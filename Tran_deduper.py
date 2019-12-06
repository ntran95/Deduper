#!/usr/bin/env python

import re
import argparse


#test file
#sam_file = "/Users/GioiTran/Documents/shell/Bi624/Deduper/input_sorted_unit_test.sam"
#sam_file = "/projects/bgmp/shared/deduper/test.sam"
sam_file = "/projects/bgmp/ntran2/test_sorted.sam"

#known barcode file
#UMI_file = "/Users/GioiTran/Documents/shell/Bi624/Deduper/STL96.txt"
UMI_file = "/home/ntran2/bgmp/Bi624/Deduper/STL96.txt"

#output file 
#output_file = "/Users/GioiTran/Documents/shell/Bi624/Deduper/output.sam"
output_file = "/home/ntran2/bgmp/Bi624/Deduper/output.sam"


def get_arguements():
    ```This function handles all argparse agruements i.e, specifying input sam files, UMI file, and pair-ended files (optional)```
    parser = argparse.ArgumentParser(description="flags are set up to specify the desired (sorted) input sam file, UMI file, and an optional flag for pair-ended reads ")
    parser.add_argument("-f", "--R1_file", help="this argument specifies the R1 file", type =str, required=True)
    parser.add_argument("-p", "--R2_file", help="this argument specifies the R2 file", type =str, required=True)
    parser.add_argument("-u", "--R3_file", help="this argument specifies the R3 file", type =str, required=True)
    parser.add_argument("-h", "--R4_file", help="this argument specifies the R4 file", type =str, required=True)
   
    return parser.parse_args()

args=get_arguements()
a=args.R1_file
b=args.R2_file
c=args.R3_file
d=args.R4_file
q=args.avg_qscore

def get_UMI(file):
    with open(UMI_file, "rt") as fh:
        known_UMI = []
        
        for line in fh:
            known_UMI.append(line.rstrip('\n'))
    return(known_UMI)

def parsing_QNAME(qname):
    qname_list = qname.split(":")
    #print(qname_list)
    return qname_list[-1] #get the last element in the list (where barcode is)

def forward(CIGAR_str, POS):
    cigar_list = re.findall('(\d+)([MNS])', CIGAR_str )
    
    if "S" in cigar_list[0]:
        fw_softclip = (cigar_list[0])
        
        adj_POS = (POS) - int(fw_softclip[0])
        return(adj_POS)
    if "S" not in cigar_list[0]:
        adj_POS = POS
        return(adj_POS)
         
    


def reverse(CIGAR_str, POS):
    cigar_list = re.findall('(\d+)([MNSD])', CIGAR_str )
    #print(cigar_list)
    if "S" in cigar_list[0]:
        del cigar_list[0]
        #print(cigar_list[0])
        for char in cigar_list:
            if "N" or "D" or "M" or "S" in char:
                POS = POS + int(char[0])
                #print(char[0])                
    else:
        for char in cigar_list:
            if "N" or "D" or "M" or "S" in char:
                POS = POS + int(char[0])
                #print(char[0])
    return(POS)
            
    


file_out=open(output_file, "a")

with open(sam_file, "rt") as file_handle:
    
    record_set_fw = set()
    record_set_rv = set()
    
    umi_list = get_UMI(UMI_file)
    previous_chromosome = 0
    
    #print(len(umi_list))
    
    for line in file_handle:
    
        if "@" in line:
            file_out.write(line)
            #pass
        else:
            ref_rec = line.split()
            
            next_chromosome = ref_rec[2]
            #group sets by chromosome
            
            if next_chromosome != previous_chromosome:
                previous_chromosome = next_chromosome
                record_set_fw.clear()
                record_set_rv.clear()
                
            
            
            if parsing_QNAME(ref_rec[0]) in umi_list and ((int(ref_rec[1]) & 16) != 16):
                
                adj_5prime_pos_fw = forward(ref_rec[5], int(ref_rec[3]))
                #add qname, chromosome, and adjusted 5' position into a tuple
                record_tuple_fw=  (parsing_QNAME(ref_rec[0]) ,ref_rec[2],adj_5prime_pos_fw, ref_rec[5])
                #print(record_tuple_fw)
                #record_set.add(record_tuple)
                #print(adj_5prime_pos)
                #print(record_set)
                
                if record_tuple_fw not in record_set_fw:
                    record_set_fw.add(record_tuple_fw)
                    
                    file_out.write(line)
                    
                                    
            #account for adjusted position if bitwise flag is reverse    
            if parsing_QNAME(ref_rec[0]) in umi_list and ((int(ref_rec[1]) & 16) == 16):
                adj_5prime_pos_rv = reverse(ref_rec[5], int(ref_rec[3]))
                
                record_tuple_rv = (parsing_QNAME(ref_rec[0]), ref_rec[2], adj_5prime_pos_rv, ref_rec[5])
                
                
                if record_tuple_rv not in record_set_rv:
                    record_set_rv.add(record_tuple_rv)
                    file_out.write(line)
                   
                

file_out.close()
            
           
            
        
                
            
                
      
        

        
    