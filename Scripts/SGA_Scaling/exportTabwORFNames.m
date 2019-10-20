function status = exportTabORFNames(scoremat,filename,annotateMut)
%exports the data in a tab-delimited form suitable for
%clustering with the Cluster 3.0 algorithm using ORFs instead of gene names. annotateMut is used as flag to
%indicate that the exported gene names should indicate the type of mutation
%used. A value of 1 for annotateMut indicates that all mutations should be
%annotated. A value greater than one indicates that all nondeletion
%mutations should be annotated.
%
%written by Sean Collins (2006) as part of the EMAP toolbox


if nargin<2
    filename='data for clustering.txt';
end
if nargin<3
    annotateMut=2;
end
fn=fieldnames(scoremat);
fid= fopen(filename,'w');
fprintf(fid,'Gene');
len1=length(scoremat.rowlabels);
len2=length(scoremat.collabels);
for i=1:len2
    orfname = char(scoremat.collabels(i));
    %genename = char(getGeneName(scoremat,orfname));
    genename = char(scoremat.collabels(i));
    if (strcmp(genename,''))
        genename = char(scoremat.collabels(i));
    end
    if annotateMut>0 && ismember('colMut',fn)
        if annotateMut==1 || ~strcmpi(scoremat.colMut(i),'DELETION')
            %genename=[genename ' - ' char(scoremat.colMut(i))];
            genename=[genename ' - ' orfname];
        end
    end
    fprintf(fid,'\t%s',genename);
end
fprintf(fid,'\n');
for i=1:len1
    orfname = char(scoremat.rowlabels(i));
    genename = char(scoremat.rowlabels(i));
    %genename = char(getGeneName(scoremat,orfname));
    if (strcmp(genename,''))
        genename = char(scoremat.rowlabels(i));
    end
    if annotateMut>0 && ismember('rowMut',fn)
        if annotateMut==1 || ~strcmpi(scoremat.rowMut(i),'DELETION')
            %genename=[genename ' - ' char(scoremat.rowMut(i))];
            genename=[genename ' - ' orfname];
        end
    end
    fprintf(fid,'%s',genename);
    for j=1:len2
        if (isnan(scoremat.data(i,j)))
            fprintf(fid,'\t');
        else
            fprintf(fid, '\t%f',scoremat.data(i,j));
        end
    end
    fprintf(fid,'\n');
end
status=fclose(fid);
