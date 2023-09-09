import csv
import orange
from Bio import Entrez, SeqIO
import os
import subprocess
import time


def method_2_prediction(sequence):
    pwm_diff = [[0.470003629245736,0.628608659422374,0.191055236762709,0.969400557188104,0.887303195000903,0.594707107746693,1.12846525181779,1.09861228866811,0.810930216216329,1.09861228866811,0,0,0.550046336919272,0,0.470003629245736,0.628608659422374,0.191055236762709,0.969400557188104,0.887303195000903,0.594707107746693,1.12846525181779,1.09861228866811,0.810930216216329,1.09861228866811], [0.287682072451781,0.693147180559945,0,0.231801614057324,1.04145387482816,0.371563556432483,0.961411167154625,0.587786664902119,0.204794412646013,0,0.559615787935423,0.0953101798043249,0.122602322092332,0.424883193965266,0.287682072451781,0.693147180559945,0,0.231801614057324,1.04145387482816,0.371563556432483,0.961411167154625,0.587786664902119,0.204794412646013,0], [0,0.405465108108164,0,0.476924072090309,0.481838086892738,0.594707107746693,0.390866308687012,1.01856958099457,0.300104592450338,0.916290731874155,0.287682072451781,0.146603474191875,0.424883193965266,0.367724780125317,0,0.405465108108164,0,0.476924072090309,0.481838086892738,0.594707107746693,0.390866308687012,1.01856958099457,0.300104592450338,0.916290731874155], [0,0,0.362905493689368,0,0,0,0,0,0,0.0339015516756814,0.559615787935423,0.0953101798043249,0,0.262364264467491,0,0,0.362905493689368,0,0,0,0,0,0,0.0339015516756814]];
    pwm_p90 = [[0.847297860387204,0.998528830111127,1.14513230430300,1.01160091167848,0.641853886172395,0.546543706368070,0.737598943130779,0.650587566141149,1.99243016469021,1.09861228866811,0,0.0606246218164348,0.348306694268216,0,0.847297860387204,0.998528830111127,1.14513230430300,1.01160091167848,0.641853886172395,0.546543706368070,0.737598943130779,0.650587566141149,1.99243016469021,1.09861228866811],[0.271933715483642,0,0.200670695462151,0.606135803570316,0.305381649551182,0.236388778064230,0.832909122935104,0.427444014826940,0.382992252256106,0,0.117783035656383,0,0,0.133531392624523,0.271933715483642,0,0.200670695462151,0.606135803570316,0.305381649551182,0.236388778064230,0.832909122935104,0.427444014826940,0.382992252256106,0],[0.336472236621213,0.0540672212702758,0,0.146603474191875,0.0540672212702758,0.171850256926659,0.302280871872934,0.737598943130779,0.0465200156348929,0.538996500732687,0.117783035656383,0.348306694268216,0.0606246218164348,0.0645385211375712,0.336472236621213,0.0540672212702758,0,0.146603474191875,0.0540672212702758,0.171850256926659,0.302280871872934,0.737598943130779,0.0465200156348929,0.538996500732687],[0,0.111225635110224,0.451985123743057,0,0,0,0,0,0,0.470003629245736,0.492476485097794,0.0606246218164348,0.0606246218164348,0,0,0.111225635110224,0.451985123743057,0,0,0,0,0,0,0.470003629245736]];
    temp=''
    b_ret=bayes_classifier(sequence)
    p90_score=pwm_2_pred(pwm_p90, sequence)
    diff_score=pwm_2_pred(pwm_diff,sequence)
    pwm_cls=pwm_3_classification(p90_score, diff_score)
    p_ret=[pwm_cls, p90_score, diff_score] #type, p90 score, diff score
    #return p_ret
    return [p_ret,b_ret]
    #return p_ret


def pwm_2_pred(pwm=[], seq=''): # to gain the p90/diff score for the seq by pwm (diff pwm or p90 pwm) provided
    #partial_seq=seq[seq.find('GACTC')-14:seq.find('GACTC')]+seq[seq.find('GACTC')+6:seq.find('GACTC')+16]
    partial_seq=seq[seq.find('GACTC')-14:seq.find('GACTC')]+seq[seq.find('GACTC')-14:seq.find('GACTC')-4]
    if len(partial_seq)<24:
        return 99
    matrix_1 = [[0 for col in range(24)] for row in range(4)]
    matrix_2 = [[0 for col in range(24)] for row in range(4)]
    matrix_score = [0 for x in range(24)]
    score =0
    for j in range(24):
        if partial_seq[j]=='A':
            matrix_1[0][j]=1
        if partial_seq[j]=='T':
            matrix_1[1][j]=1
        if partial_seq[j]=='G':
            matrix_1[2][j]=1
        if partial_seq[j]=='C':
            matrix_1[3][j]=1
    for i in range(24):
        for j in range(4):
            matrix_2[j][i]=matrix_1[j][i]*pwm[j][i]
    for i in range(24):
        for j in range(4):
        	matrix_score[i]=matrix_score[i]+matrix_2[j][i]
    score= sum(matrix_score)
    return score

def pwm_3_classification(P90_score, Diff_score): # to gain the predicted classification of sequence by the p90 and diff score
    if Diff_score > -0.8510*P90_score+15.8764 :
        return 'bad'
    else:
        return 'good'

def import_pwm_csv (infile,outfile,pwm_p90,pwm_diff): # input csv file with template id and template sequence, output the csv file with pwm mwthod predicted p90 score and diff score and predicted classification
    fout=open(outfile,'w', newline='')
    writer=csv.writer(fout)
    with open(infile, newline='') as fin:
        reader = csv.reader(fin)
        for row in reader:
            p90_score=pwm_2_pred(pwm_p90,row[1])
            diff_score=pwm_2_pred(pwm_diff,row[1])
            seq_class=pwm_3_classification(p90_score,diff_score)
            writer.writerow([row[0],row[1],p90_score,diff_score,seq_class])
    fout.close()

def bayes_classifier(sequence):
    f=open('header.tab','r')
    feature_header=f.read()
    f.close()
    f=open('new_seq.tab','w')
    f.write(feature_header)
    f.write('unknow\t')
    features=bayes_pre(sequence)
    for feature in features:
        f.write(str(feature))
        f.write('\t')
    f.close()
    data = orange.ExampleTable("data3_v5.tab")
    learner = orange.BayesLearner()
    classifier = learner(data)
    new_data=orange.ExampleTable("new_seq.tab")
    p = [str(classifier(new_data[0])),orange.GetProbabilities]
    return p

def bayes_pre(sequence):
    #sequence='CTACTTACTGGGATGACTCTCTACTTACTG'
    #feature=['AT-1', 'T-2', 'CC-3']
    #result=[0 for i in range(len(features))]
    features=['C-8',	'A-9',	'A-10',	'T-10',	'GA-8',	'C-7',	'A-1',	'G-8',	'C-6',	'GA-9',	'G-3',	'CC-7',	'AG-9',	'CT-9',	'C-10',	'GG-3',	'CG-10',	'CG-6',	'A-3',	'CC-8',	'AG-1',	'A-6',	'C-9',	'GA-5',	'AC-7',	'T-9',	'AC-10',	'AAC-11',	'A-8',	'CGA-7',	'TT-9',	'AG-8',	'TG-2',	'TC-6',	'T-2',	'AC-12',	'TA-10',	'TA-7',	'AGT-1',	'CCT-8',	'CC-1',	'CT-7',	'CA-2',	'GAC-9',	'CT-8',	'TGG-3',	'AA-11',	'AG-12',	'GG-1',	'CG-7',	'G-7',	'AA-10',	'GAC-5',	'GG-4',	'AG-6',	'A-7',	'G-2',	'CTC-6',	'GT-8',	'GAG-8',	'AC-3',	'AA-3',	'AA-1',	'CAA-9',	'GC-5',	'GT-7',	'TAC-7',	'ACT-10',	'GA-2',	'GGAC-4',	'CCG-5',	'T-1',	'CTG-9',	'ATC-4',	'GTT-8',	'A-12',	'AACA-11',	'CGG-5',	'GGA-1',	'A-2',	'TG-8',	'CAA-10',	'ACA-3',	'TGG-2',	'TGA-7',	'TAG-11',	'GAC-2',	'CC-9',	'AT-1',	'ACTC-10',	'GAG-5',	'CCT-7',	'CTC-9',	'TC-4',	'AA-8',	'CCTC-8',	'ACC-12',	'ACCG-12',	'C-11',	'T-8',	'CAG-7',	'CCGC-1',	'CGA-8',	'CCG-1',	'TGGA-3',	'GAC-6',	'CG-3',	'CGAC-5',	'ACC-6',	'TA-9',	'ACG-9',	'GACC-5',	'TAG-10',	'TTAG-9',	'GGT-3',	'GTTA-8',	'GG-2',	'GGA-7',	'TCT-6',	'TC-10',	'CT-11',	'CA-4',	'GT-4',	'TCG-4',	'G-11',	'A-4',	'CGT-4',	'ATGT-9',	'TCTG-8',	'TA-2',	'TCT-8',	'AAAC-10',	'TGC-4',	'CCT-4',	'G-9',	'CCC-8',	'C-4',	'T-5',	'GGA-4',	'GCC-3',	'ACGA-6',	'GAT-5',	'AC-9',	'CT-10',	'GAA-8',	'GT-1',	'TG-10',	'GTT-11',	'TTG-9',	'TG-6',	'CG-8',	'CTC-5',	'TAG-7',	'TA-1',	'CTA-10',	'TAA-10',	'CGC-6',	'GTA-1',	'AGTG-12',	'AAG-1',	'GGAC-1',	'CGCC-10',	'TC-12',	'GTA-6',	'AGT-12',	'CTG-8',	'CA-7',	'TC-9',	'ACG-5',	'AAC-6',	'AGA-10',	'AA-4',	'A-5',	'CT-1',	'GTG-7',	'GGT-2',	'GT-5',	'TAC-11',	'AC-6',	'TGA-8',	'GTG-4',	'CGAA-7',	'AT-6',	'ACG-6',	'C-13',	'GT-9',	'ATG-1',	'AAA-8',	'C-1',	'AG-2',	'AGA-2',	'TTA-9',	'CTC-11',	'TC-5',	'GAC-8',	'AT-4',	'CCC-7',	'ACA-1',	'AGA-9',	'AG-4',	'GCC-7',	'AGA-6',	'GAA-7',	'TCC-4',	'TAC-2',	'GTG-12',	'GGG-10',	'GC-7',	'GTA-9',	'TTC-12',	'AGA-8',	'GT-3',	'CTCC-9',	'GACT-9',	'TTCG-12',	'CCGC-5',	'TGGA-6',	'GTGG-12',	'GTAA-9',	'G-5',	'CGC-2',	'CT-5',	'CTA-5',	'CTG-7',	'TAAT-10',	'CTAC-4',	'CGC-10',	'AAC-10',	'AGAA-8',	'TT-4',	'TA-11',	'GCTA-9',	'TTT-4',	'GGG-6',	'CCTG-7',	'CGAG-7',	'TTGT-9',	'TGGG-2',	'CC-13',	'CG-5',	'CCGA-13',	'CCG-13',	'TGG-8',	'AAG-9',	'CGAC-14',	'CGTG-12',	'CGAC-8',	'C-14',	'CGT-12',	'CTGG-7',	'TGT-6',	'CGA-14',	'CG-14',	'TC-7',	'ATC-3',	'GCC-11',	'AG-10',	'GTG-5',	'GCT-11',	'TCA-6',	'TCA-5',	'GCA-4',	'TCC-10',	'GAAG-8',	'TCCT-7',	'ACT-3',	'CGG-10',	'GGCG-8',	'TAG-1',	'AAC-8',	'ACC-8',	'TT-10',	'ACAA-3',	'CG-2',	'GGT-7',	'AT-9',	'GGC-8',	'AT-10',	'CTT-11',	'TGA-4',	'TACG-5',	'CG-4',	'ACGA-7',	'CGA-1',	'AGG-6',	'AC-2',	'TGG-9',	'GCA-6',	'TACG-12',	'ACCT-7',	'GCA-8',	'TCTG-12',	'GCT-9',	'ATG-9',	'CCGC-9',	'TT-12',	'TCT-12',	'ACG-7',	'CTC-3',	'GCG-7',	'TTA-11',	'GCG-2',	'TAC-12',	'AGG-1',	'AC-1',	'GTA-7',	'CCAG-12',	'TAGC-11',	'G-10',	'CCTC-4',	'TAC-1',	'CCA-12',	'TG-7',	'GG-6',	'ACT-9',	'CGT-8',	'CCC-6',	'ATC-5',	'CTAG-10',	'GC-2',	'CGC-1',	'GGG-3',	'ATG-5',	'TACC-1',	'GGT-6',	'ATGA-1',	'CA-9',	'AAA-11',	'CC-5',	'GC-3',	'AA-13',	'TAGA-1',	'GCG-5',	'GCG-11',	'GA-11',	'GCG-9',	'TCAT-6',	'GCT-7',	'GCGA-5',	'CTG-6',	'ATA-4',	'TGC-6',	'AGAT-9',	'CGT-6',	'GGA-11',	'TAC-5',	'AAGA-13',	'ATC-11',	'GAAC-10',	'GCT-5',	'GAT-1',	'CAC-7',	'AAG-13',	'ACC-3',	'GATG-1',	'ATAC-1',	'ACA-9',	'ATG-2',	'GGG-1',	'CT-12',	'CTA-11',	'CAC-4',	'AGC-9',	'TTT-5',	'GG-8',	'TA-8',	'ACTA-10',	'AA-5',	'TGA-6',	'ATA-1',	'ACG-1',	'TCGC-1',	'CAT-2',	'CA-11',	'AC-8',	'TGCT-4',	'CAT-11',	'TCC-7',	'ACG-3',	'AAC-9',	'TACA-11',	'AGC-8',	'GAA-10',	'GTT-2',	'GAT-3',	'GTC-9',	'GCTG-6',	'ACC-2',	'ACC-7',	'CAA-4',	'CCT-10',	'CCG-9',	'GCC-9',	'CCG-8',	'AACC-11',	'CCCG-6',	'CTC-4',	'GT-6',	'CCG-12']
    result=[]
    for feature in features:
        ind=feature.index('-')
        simbol=feature[0:ind]
        start=feature[ind+1:]
        if cmp(sequence[int(start)-1:int(start)-1+len(simbol)],simbol)==0:
             result.append('y')
        else:
            result.append('n')
    return result

def trigen(temp_filename):
    trigger_filename=str(int(time.time()))+'.txt'
    outfile=open(trigger_filename,'w')
    for cur_record in SeqIO.parse(open(temp_filename), "fasta") :
        if cur_record.seq.find('GACTC')==-1:
            return 0 #'error, no GACTC in the templte sequence'
	    else:
	        trigger=cur_record.seq[:cur_record.seq.find('GACTC')]
                print(trigger)
	        print >>outfile, '>tri_',
	        print >>outfile, str(cur_record.id)
	        print >>outfile, trigger
    return trigger_filename

def unafold_pred(temp_filename,trig_filename):
    dic_therm={}
    name=[]
    for cur_record_o in SeqIO.parse(open(temp_filename), "fasta") :
        name.append(cur_record_o.id)
    num=0
    nc="D:/NAR_tool/UNAFold/bin/melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.00000005 "+trig_filename+" "+temp_filename+" >DGHS.txt\n"
    nc2="D:/NAR_tool/UNAFold/bin/melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.00000005 "+temp_filename+" "+temp_filename+" >DGHS_2.txt\n"
    helix_c1="D:/NAR_tool/UNAFold/bin/UNAFold.pl "+temp_filename+" -n DNA"
    helix_c2="D:/NAR_tool/UNAFold/bin/ct2rnaml "+temp_filename+".ct"
    subprocess.Popen(nc, shell=True)
    print nc,nc2,helix_c1,helix_c2
    os.system(nc)
    os.system(nc2)
    os.system(helix_c1)
    os.system(helix_c2)
    infile = open("DGHS.txt")
    infile2 = open("DGHS_2.txt")
    line = infile.readline()
    line2 = infile2.readline()
    print line, line2
    helix_all=cal_helix(temp_filename+'.rnaml')
    #print "<tr><td>rank</td><td>type</td><td>location</td><td>","trigger", "</td><td>","template","</td>","<td>DH</td>","<td>DS</td>","<td>tem-tem_TM</td>","<td>TM</td>","<td>number of helix</td>","<td>bonds</td>"
    i=0
    while 1:
    	line = infile.readline().rstrip()
    	line2 = infile2.readline().rstrip()
	print line, line2
    	if line !='' and line2 !='':
    	    line=line.replace('-','')
    	    ener=line.split("\t")
    	    ener2=line2.split("\t") #template-template TM
    	    TM=float(ener[1])/((float(ener[2])-(1.987*log(0.00000005)))/1000)-273.15 #trigger-template TM
    	    bondadd=0;
    	    one=helix_all[i]
    	    list=one[1]
    	    if one[0]>0:
    	        for j in range(len(list)):
    	              bondadd=bondadd+int(list[j])
            dic_therm[name[i]]=[TM,float(ener2[3]),bondadd] # Trigger-template TM, template-template TM,bonds
            num=num+1
            i=i+1
    	if line == "":
            # Reached the end of the record or end of the file
    	    break

    infile.close()
    infile2.close()
    return dic_therm

def cal_helix(file):
    input_file=open(file, 'r')
    site=-1
    leng=[]
    count=0
    all_cal_helix=[]
    line=input_file.readline()
    while (line!=''):
        line=input_file.readline()
        if (line.find('seq-data')!=-1):
            count=0
            leng=[]
            #print "new"
        #if (line.find('<model id')!=-1):
            #a=line
        elif (line.find('End of folding 1 for')!=-1):
            #print "end"
            count=count/2
            #print leng
            all_cal_helix.append([count,leng])
            count=0
            leng=[]
            continue
        elif (line.find('helix')!=-1): #count the number of helix
            count=count+1
        elif (line.find('length')!=-1): #count the length of bonds in this helix
            t1=line.find('>')
            t2=line.find('</')
            leng.append(line[t1+1:t2])
            
        
    return all_cal_helix
    input_file.close()


if __name__ == "__main__":
	seq1='CCTACGACTGAACAGACTCTCCTACGACTG' # good
	seq2='CCTACGACTTAACAGACTCTCCTACGACTT' # good
	seq3='CCTACTACTGAACAGACTCTCCTACTACTG' # good
	seq4='CCTGCGACTGAACAGACTCTCCTGCGACTG' # bad
	seq5='GGGGAAATAGGTGAGACTCTGGGGAAATAG' # bad
	seq6='TGGCGTGAAAAACGGACTCTTGGCGTGAAA' # bad
	print method_2_prediction(seq1)
	print method_2_prediction(seq2)
	print method_2_prediction(seq3)
	print method_2_prediction(seq4)
	print method_2_prediction(seq5)
	print method_2_prediction(seq6)

#import orange data = orange.ExampleTable("voting")
#classifier = orange.BayesLearner(data)
#print "Possible classes:", data.domain.classVar.values
#print "Probabilities for democrats:"
#for i in range(5):
    #p = classifier(data[i], orange.GetProbabilities)
    #print "%d: %5.3f (originally %s)" % (i+1, p[1], data[i].getclass())