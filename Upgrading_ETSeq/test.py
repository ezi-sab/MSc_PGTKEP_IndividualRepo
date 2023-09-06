from PyQt4 import QtCore, QtGui
import re
from Bio import Entrez, SeqIO
outfilename='1329730321.txt'
seq_list={}
def check_template(template_seq):
    template_seq=template_seq.upper()
    GACTC_start=[i - 1 for i in range(len(template_seq)) if template_seq.startswith('GACTC', i - 1)]
    if re.search('[^a^g^t^c^A^G^T^C]',template_seq):
        return "some leters other than ATGC/atgc shows in template sequences!"
    else:
        if len(GACTC_start)!=1:
            return "more than one complementary nicking enzyme site or no complementary nicking enzyme site!"
        else:
            trigger=template_seq[:GACTC_start[0]-4]
            if len(trigger)<10:
                return "template too short!"
    return 1






for cur_record_o in SeqIO.parse(open(outfilename), "fasta"):
    if seq_list.has_key(str(cur_record_o.id)):
	QtGui.QMessageBox.about(self, "Exist", str(cur_record_o.id)+" is existed before, so ignored!")
    else:
	if check_template(str(cur_record_o.seq))!=1:
		t=check_template(str(cur_record_o.seq))
	else:
		seq_list[str(cur_record_o.id)]=str(cur_record_o.seq)
