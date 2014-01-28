import sys
import  os
import re
"""
##################### match mature and detected miRNAs ############ 
list_detected_mirnas = []
with open("intersect_pre-mirANDdpi/list_mirnas_detected_hsa_only.txt", "r") as fd:
    for line in fd:
        list_detected_mirnas.append(line.strip("\n"))
#print len(list_detected_mirnas)        


list_tuple = []
with open("intersect_pre-mirANDdpi/maturemirna_precursor_list.txt", "r") as fd2:
    fd2.readline()
    for line in fd2:
        res = re.search("^(hsa-.+)\t(hsa.+)$", line)
        if res:
            mir = res.group(1)
            pre = res.group(2)
            t = (mir, pre)
            list_tuple.append(t)
#print len(list_tuple)

list_tuple2 = [] # tuple(mature, prec)
for t in list_tuple:
    m_t = (t[0], t[1])
    for m in list_detected_mirnas:
        if m == t[1]+"-3p" and "-3p" in t[0]:
            m_t = (m, t[1])
        elif m == t[1]+"-5p" and "-5p" in t[0]:
            m_t = (m, t[1])
    list_tuple2.append(m_t)

#print len(list_tuple2)
for t in list_tuple2:
    print t[0]+"\t"+t[1]
    
    
"""
######################## match mature miRNA and precursor #################
dict_precursor = {}
with open("intersect_pre-mirANDdpi/intersect_pre-mirANDdpi_list.txt", "r") as fd:
    for line in fd:
        res = re.search("^(hsa-[mirlet]+-[0-9a-z]+)-?[0-9]?$", line)        
        if res:
            dict_precursor[line.strip("\n")] = res.group(1)

#print dict_precursor

print "mature_mirna\tprecursor"
with open("all_miRNAs_in_miRBase.txt", "r") as fd2:
    for line in fd2:
        res = re.search("^(hsa-[mirlet]+-[0-9a-z]+)-?[53]?p?$", line)        
        if res:
            for k,v in dict_precursor.items():
                if res.group(1) == v:
                    print line.strip("\n")+"\t"+k
