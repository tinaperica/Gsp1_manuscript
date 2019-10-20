email from Hannes on 20181105
"I put it together in 2009 based on unaveraged E-MAP files I got from Sean Collins and Mike Shales. They're unaveraged because I wanted to first scale them, then merge them, and finally average them (using the compute averaged scores function in the EMAP toolbox). For these a linear scaling was fine, but the Boone data complicated things...

SeanUberAryScreens = ubermap screens scores 07-01-08.mat  (I think this is a bunch of extra screens that we never published)

TF = TFunaveragedFromText.mat  (transcription factor emap published with Hao Li)

cF3 = chrom bio unavg corrected 06-21-07.mat (chromosome biology emap)

cc2I = cell cycle screens 04-09-08.mat  (again screens I'm not sure if they were published)

combRNAP = RNAPcombinedUnaveraged.mat  (RNA processing published with Christine Guthrie)

endoUnaveraged = endoUnaveraged.mat   (endocytosis emap - think this was published at some point... maybe with Bob Farese)

espT =  esp unavg corrected 06-19-07.mat   (early secretory pathway EMAP - the first EMAP together with chromosome biology)

kinaseJoinedUnAvGd = kinaseJoinedUnaveraged(d4).mat  (kinase emap with kevan shokat)

scorematQC1p = surmaQC1unavg.mat   (lipid emap we published in maybe 2010-2012.  first author: Michal Surma)

slp = mito unavg 04-09-09.mat   (mitochondrial emap. sean collins. not sure if/when this was published. Sean Collins website would have it listed)


The following were not 1536 format:
TFunaveragedFromText.mat
Endocytosis
Chromosome biology
ESP
Kinase
RNAPcombinedUnaveraged.mat

These were 1536 format:
Mitochondrial
Surma/Lipid"


Individual screens were averaged and exported from Matlab with:
load('Users/tina/Desktop/Uber_091216.mat');

scorematQC1p_avg = computeAveragedScores_hannes4(scorematQC1p);
exportTabwORFNames(scorematQC1p_avg, '~/Desktop/scorematQC1p.txt');
exportForCluster3_0(scorematQC1p_avg, '~/Desktop/scorematQC1p_gene_names.txt');


slp_avg = computeAveragedScores_hannes4(slp);
exportTabwORFNames(slp_avg, '~/Desktop/slp.txt');
exportForCluster3_0(slp_avg, '~/Desktop/slp_gene_names.txt');


TF_avg = computeAveragedScores_hannes4(TF);
exportTabwORFNames(TF_avg, '~/Desktop/TF.txt');
exportForCluster3_0(TF_avg, '~/Desktop/TF_gene_names.txt');


cF3_avg = computeAveragedScores_hannes4(cF3);
exportTabwORFNames(cF3_avg, '~/Desktop/cF3.txt');
exportForCluster3_0(cF3_avg, '~/Desktop/cF3_gene_names.txt');


cc2l_avg = computeAveragedScores_hannes4(cc2l);
exportTabwORFNames(cc2l_avg, '~/Desktop/cc2l.txt');
exportForCluster3_0(cc2l_avg, '~/Desktop/cc2l_gene_names.txt');


combRNAP_avg = computeAveragedScores_hannes4(combRNAP);
exportTabwORFNames(combRNAP_avg, '~/Desktop/combRNAP.txt');
exportForCluster3_0(combRNAP_avg, '~/Desktop/combRNAP_gene_names.txt');


endoUnaveraged_avg = computeAveragedScores_hannes4(endoUnaveraged);
exportTabwORFNames(endoUnaveraged_avg, '~/Desktop/endo.txt');
exportForCluster3_0(endoUnaveraged_avg, '~/Desktop/endo_gene_names.txt');


espT_avg = computeAveragedScores_hannes4(espT);
exportTabwORFNames(espT_avg, '~/Desktop/espT.txt');
exportForCluster3_0(espT_avg, '~/Desktop/espT_gene_names.txt');


kinaseJoinedUnAvGd_avg = computeAveragedScores_hannes4(kinaseJoinedUnAvGd);
exportTabwORFNames(kinaseJoinedUnAvGd_avg, '~/Desktop/kinase.txt');
exportForCluster3_0(kinaseJoinedUnAvGd_avg, '~/Desktop/kinase_gene_names.txt');


SeanUberAryScreens_avg = computeAveragedScores_hannes4(SeanUberAryScreens);
exportTabwORFNames(SeanUberAryScreens_avg, '~/Desktop/SeanUberAryScreens.txt');
exportForCluster3_0(SeanUberAryScreens_avg, '~/Desktop/SeanUberAryScreens_gene_names.txt');

Code for computeAveragedScores_hannes4 function is included in the folder as computeAveragedScores_hannes4.m

