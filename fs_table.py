## manual script to calculate fsc, fsd

import pandas as pd
import os
import sys

sh=os.system

sample_prefix=sys.argv[1]
outpath=sys.argv[2]

out_fsc_window="{}/{}_fsc_window.bed.gz".format(outpath,sample_prefix)
# out_fsc_window = "./fsc/fsc_window/{}_fsc_window.bed.gz".format(sample_prefix2)
test_input = pd.read_csv(out_fsc_window,sep="\t",header=None)
sample_prefix2=sample_prefix
# fsd
test_input["region"] = test_input[4]+":"+test_input[5].map(str)+":"+test_input[6].map(str)
thres = [100,150,220,300]
thres_label = ["short","mid","long"]
test_input["len_thres"] = pd.cut(test_input[3],thres,labels=thres_label)
test_output=pd.DataFrame(test_input.groupby(["region","len_thres"]).aggregate("count").iloc[:,0])

test2 = test_output.reset_index(level=None).pivot(index="region",columns="len_thres",values=0)
fsc = pd.DataFrame((test2["long"]+test2["short"])/test2["mid"])
fsc.columns = [sample_prefix2]

thres = [i for i in range(100,301,5)]
thres_labels = ["{0}-{1}".format(i,i+5) for i in range(100,300,5)]


test_input['range']=pd.cut(test_input.iloc[:,3],thres,labels=thres_labels)
fsd =pd.DataFrame(test_input.groupby([0,"range"]).aggregate("count").iloc[:,0])
fsd.columns = [sample_prefix2]

fsd.to_csv("{}/fsd_table.txt".format(outpath),sep='\t',index=True,header=True)
fsc.to_csv("{}/fsc_table.txt".format(outpath),sep='\t',index=True,header=True)
