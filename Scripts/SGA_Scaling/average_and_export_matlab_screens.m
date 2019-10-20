load('Data/SGA_Scaling/Uber_091216.mat');

cF3_avg = computeAveragedScores_hannes4(cF3);
exportTabwORFNames(cF3_avg, 'Data/SGA_Scaling/cF3.txt');
exportForCluster3_0(cF3_avg, 'Data/SGA_Scaling/cF3_gene_names.txt');

load('Data/SGA_Scaling/SGA_NxN_unaveraged.mat');

SGA_NxN_avg = computeAveragedScores_hannes4(SGA_NxN_unaveraged);
exportTabwORFNames(SGA_NxN_avg, 'Data/SGA_Scaling/SGA_NxN_avg.txt');

load('Data/SGA_Scaling/SGA_unaveraged.mat');

SGA_avg = computeAveragedScores_hannes4(SGA_unaveraged);
exportTabwORFNames(SGA_avg, 'Data/SGA_Scaling/SGA_avg.txt');
