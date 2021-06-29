library(tidyverse)
library(cowplot)
library(factoextra)
library(dendextend)
library(ComplexHeatmap)
library(circlize)

source('ucsf_colors.R')

load('Data/spitzemapko_for_corr.rda')


# get the mutants and 278 alleles we are using for Figure 4 from filtered_v6_correlations.RData
load('Data/filtered_v6_correlations.RData')
gsp1_alleles <- filtered_correlations %>% pull(query_uniq1) %>% unique()
query_alleles <- filtered_correlations %>% pull(query_uniq2) %>% unique()

length(gsp1_alleles) # 22 Gsp1 mutants
length(query_alleles) # 278 yeast genes

# we only want to compute correlations for those queries, since we are assessing how much allelic
# dependencies between library genes matter
df_GI <- spitzemapko_for_corr %>% 
    filter(query_allele_name %in% c(gsp1_alleles, query_alleles)) %>% 
    select(-weight)

# this confirms that we are using the 1129 library genes that 
df_GI %>% pull(array_ORF) %>% unique() %>% length()

# why do we have multiple values for crm1_damp, yrb2_damp, and yrb2_damp???
df_GI %>% 
    group_by(query_allele_name, array_ORF) %>% 
    summarise(n=n()) %>% 
    filter(n>1) %>% 
    select(-array_ORF) %>% 
    unique()

spitzemapko_for_corr %>% 
    filter(query_allele_name %in% c('crm1_damp')) %>% 
    arrange(array_ORF)

# load the object SGA_scaled_to_EMAP
# load('Data/SGA_scaling/SGA_2016_full_scaled.rda')

# SGA_scaled_to_EMAP %>%
#     filter(query_strain_id == 'crm1_damp', array_ORF == 'YAL002W') %>% 
#     group_by(query_strain_id) %>% 
#     mutate(n=n()) %>% 
#     filter(n>1) %>%  
#     print(width=Inf)

# we see from the code above that the difference between the two measurements 
# are that some deletion mutant array strains were included in the TSA (the 
# array of TS alleles) to allow for calibration of the TSA. These array strains
# were of course also used in the deletion mutant array (DMA), and their
# scores are likely different because of the two different growth temps
# (DMA30 at 30C, TSA26 at 26C). if we use these strains, then we should
# only use the DMA30 score, but I think we should just not use DAmPs.
#
# this is not a problem for the rest of the spitzemapko, because I filtered
# out strains from the TSA (using !grepl('_tsa', array_strain_id)) in the 
# file spitzemap_construction.Rmd
df_GI <- filter(df_GI, !query_allele_name %in% c('crm1_damp', 'yrb2_damp'))



mat_GI <-
    df_GI %>% 
    pivot_wider(id_cols=query_allele_name, names_from=array_ORF, values_from=score) %>% 
    column_to_rownames('query_allele_name') %>% 
    as.matrix()


mat_cor <- cor(mat_GI, use='pairwise.complete.obs', method='pearson')

# plot the distribution of correlation values
corvals <- c(mat_cor)
qplot(corvals, geom="histogram", bins=1000) + 
    xlab('Correlations between library genes for GI with \n 22 Gsp1 mutants and 276 alleles') +
    ylab('Count') +
    theme_light() +
    theme(
        text = element_text(family = "Helvetica", size = 6),
        axis.title = element_text(size = 6), axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.05), axis.ticks.length = unit(0.05, 'cm'),
        strip.text.x = element_text(size = 6)
    )

ggsave('Revisions2/Fig_and_analysis_for_StatReviewer/corr_between_library_genes.png', height=2, width=3)



# Make df_cor is the square matrix, in table form, without the diagonal entries
# we count the number of high correlations (r > 0.5) for each ORF
df_cor <- mat_cor %>% 
    as_tibble(rownames='array_ORF1') %>% 
    pivot_longer(cols=-array_ORF1, names_to='array_ORF2', values_to='cor') %>% 
    filter(array_ORF1 != array_ORF2) %>% 
    group_by(array_ORF1) %>% 
    mutate(number_high_cor=sum(cor > 0.5))


# distribution of counts of correlations > 0.5
plot <- df_cor %>% 
    select(array_ORF1, number_high_cor) %>% 
    unique() %>% 
    ggplot(aes(x=number_high_cor)) +
    geom_bar(stat = 'count') +
    xlab('N. of correlations > 0.5\n(between library genes)') + ylab('Count') +
    theme_light() +
    theme(
        text = element_text(family = "Helvetica", size = 6),
        axis.title = element_text(size = 6), axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.05), axis.ticks.length = unit(0.05, 'cm'),
        strip.text.x = element_text(size = 6)
    )

inset <- df_cor %>% 
    select(array_ORF1, number_high_cor) %>% 
    unique() %>% 
    ggplot(aes(x=number_high_cor)) +
    geom_bar(stat = 'count') +
    xlab('N. of correlations > 0.5\n(between library genes)') + ylab('Count') +
    xlim(c(-1,25)) +
    theme_light() +
    theme(
        plot.background = element_rect(colour = "black",size = 0.2),
        text = element_text(family = "Helvetica", size = 6),
        axis.title = element_text(size = 6), axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.05), axis.ticks.length = unit(0.05, 'cm'),
        strip.text.x = element_text(size = 6)
    )

plot + cowplot::draw_plot(inset, x = 40, y = 100, width = 150, height = 500)
ggsave('Revisions2/Fig_and_analysis_for_StatReviewer/n_corr_per_library_gene.png', height=2, width=3)


# what are these two outliers?
df_cor %>% 
    select(array_ORF1, number_high_cor) %>% 
    unique() %>% 
    arrange(desc(number_high_cor))

# clearly we need to remove these top two outliers:
#   YBR156C - SL115, a subunit of the chromosomal passenger complex (CPC)
#   YFR036W - CDC26, a subunit of the Anaphase-Promoting Complex/Cyclosome (APC/C)

# remake the data structures in the exact same way after removing these library genes
df_GI <- filter(df_GI, !array_ORF %in% c('YBR156C', 'YFR036W'))

mat_GI <-
    df_GI %>% 
    pivot_wider(id_cols=query_allele_name, names_from=array_ORF, values_from=score) %>% 
    column_to_rownames('query_allele_name') %>% 
    as.matrix()

mat_cor <- cor(mat_GI, use='pairwise.complete.obs', method='pearson')

df_cor <- mat_cor %>% 
    as_tibble(rownames='array_ORF1') %>% 
    pivot_longer(cols=-array_ORF1, names_to='array_ORF2', values_to='cor') %>% 
    filter(array_ORF1 != array_ORF2) %>% 
    group_by(array_ORF1) %>% 
    mutate(number_high_cor=sum(cor > 0.5))


# now our worst offender has 20 high correlations, so we still need to pare down
df_cor %>% 
    select(array_ORF1, number_high_cor) %>% 
    unique() %>% 
    arrange(desc(number_high_cor))

# To do so, we can cluster the library genes using the correlations as a 
# distance metric, which will put the most highly correlated genes together.
# Then we can cut the tree into a certain number of clusters, and only pick 
# one gene from each cluster. This is effectively taking just one member from
# each functional cluster of genes.

pick_least_corr_ORF <- function(df) {
    # function to pick a less-dependent representative ORF
    # it takes in a dataframe like
    #     # A tibble: 12 x 2
    #        ORF     number_high_cor
    #        <chr>             <int>
    #      1 YBR023C              14
    #      2 YER149C               8
    #      3 YGR217W               9
    # and would in this case return "YER149C"
    # Note that we randomly sample which row to take, to resolve ties

    df %>% 
    filter(number_high_cor==min(number_high_cor)) %>% 
    sample_n(1) %>% 
    pull(ORF)
}

clustfn <- function(mat) {
  set.seed(3)
  cormat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  dissim <- as.dist((1 - cormat)/2)
  hc <- hclust(dissim, method = "average" )
  return(hc)
}
row.hc <- clustfn(t(mat_GI))
col.hc <- clustfn(mat_GI)

# in Figure 1D, we found cutting the tree of 1444 genes into 1200 clusters 
# (83.1% reduction) was an effective means of clustering genes together by 
# biological process so I started with this same factor of reduction for the
# 1127 genes leftover at this stage of the analysis (since we needed to use 
# genes that were present for the 276 alleles). This means we should cut into
# k=936 clusters. This did not get rid of all dependencies, so I wrote a loop
# to decrease the number k until there were no high correlations lets

# I tested a few values to narrow down the range
# I found that cutting into 830 clusters had removed all dependencies but 835 had not
# The correct number ends up being 834, as seen in the code output

# for (k in c(835)) {
for (k in seq(835, 830)) {

    print(paste0('testing ', k, ' clusters...'))
    cluster_assignments <- cutree(col.hc, k = k)

    least_dependent_ORFs <- 
        enframe(cluster_assignments, name='ORF', value='cluster_number') %>% 
        right_join(df_cor %>% 
                    select(array_ORF1, number_high_cor) %>% 
                    unique(),
                by=c('ORF'='array_ORF1')) %>% 
        nest(cluster=c(ORF, number_high_cor)) %>% 
        mutate(representative_ORF = map_chr(cluster, pick_least_corr_ORF)) %>% 
        select(cluster_number, representative_ORF) %>% 
        pull(representative_ORF)

    # now remake the df_cor, but with only this subset of ORFs
    # recalculate the number of high correlations and get the max
    # if max > 0, need smaller clusters
    
    df_cor_tmp <- mat_cor %>% 
        as_tibble(rownames='array_ORF1') %>% 
        pivot_longer(cols=-array_ORF1, names_to='array_ORF2', values_to='cor') %>% 
        filter(array_ORF1 != array_ORF2) %>% 
        filter(array_ORF1 %in% least_dependent_ORFs, array_ORF2 %in% least_dependent_ORFs) %>% 
        group_by(array_ORF1) %>% 
        mutate(number_high_cor=sum(cor > 0.5)) 
        
    max_n_high_cor <- df_cor_tmp %>% 
        select(array_ORF1, number_high_cor) %>% 
        unique() %>% 
        pull(number_high_cor) %>% 
        max()

    print(paste0('The gene with the most high correlations has ', max_n_high_cor, ' high correlation after pruning.'))

    if (max_n_high_cor == 0) {
        print('All high correlations have been removed')
        print(paste0('Correct number of clusters to make is k = ', k))
        break
    }

}

# we now have a table with no high correlations
df_nhc <- df_cor_tmp

#let's cluster it and see what it looks like
mat_nhc <-
    df_nhc %>% 
    pivot_wider(id_cols=array_ORF1, names_from=array_ORF2, values_from=cor) %>% 
    column_to_rownames('array_ORF1') %>% 
    as.matrix()

# compare it to the heatmap with high correlating genes arround
mat_whc <-
    df_cor %>% 
    pivot_wider(id_cols=array_ORF1, names_from=array_ORF2, values_from=cor) %>% 
    column_to_rownames('array_ORF1') %>% 
    as.matrix()


hm_nhc <- 
    Heatmap(mat_nhc, name = 'GI Profile correlations\nHigh correlating genes\nremoved',
            col = colorRamp2(c(-0.2, 0, 0.2), c('red', 'black', 'green')),
            show_row_dend = F, show_row_names = F,
            show_column_dend = F, show_column_names = F,
            heatmap_legend_param = list(direction='horizontal',position='bottom')
            )

hm_whc <- 
    Heatmap(mat_whc, name = 'GI Profile correlations\nHigh correlating genes\nkept',
            col = colorRamp2(c(-0.2, 0, 0.2), c('red', 'black', 'green')),
            show_row_dend = F, show_row_names = F,
            show_column_dend = F, show_column_names = F,
            heatmap_legend_param = list(direction='horizontal',position='bottom')
            )

dev.off()
pdf('Revisions2/Fig_and_analysis_for_StatReviewer/corr_between_library_genes_heatmap_no-hi-corr.pdf', 
    width=7, height=5)
draw(hm_nhc)
dev.off()

pdf('Revisions2/Fig_and_analysis_for_StatReviewer/corr_between_library_genes_heatmap_with-hi-corr.pdf', 
    width=7, height=5)
draw(hm_whc)
dev.off()


# Now let's remake Fig 1b using the smaller set of library genes


GI_cyan <- '#0BC3E8'
GI_yellow <- '#FDFB00'
GI_black <- '#000000'

# read in raw E-MAP data, process into a matrix
ORF2name <-
    data.frame('ORF'=colnames(read_tsv('Data/E-MAP/gsp1_pEMAP_avg_merged.txt')),
               'name'=colnames(read_tsv('Data/E-MAP/gsp1_pEMAP_avg_merged_gene_names.txt'))
    ) %>% 
    as_tibble()

emap <-
    read_tsv('Data/E-MAP/gsp1_pEMAP_avg_merged_gene_names.txt') %>%
    gather(-Gene, key = strain, value = score) %>%
    separate(strain, 'library_gene', sep = ' - ', remove = TRUE) %>%
    separate(Gene, c('Gene', 'mutant'), sep = ' - ') %>%
    left_join(ORF2name, by=c('library_gene'='name')) %>% 
    filter(ORF %in% least_dependent_ORFs) %>% 
    filter(! mutant %in% c('NTER3XFLAG WT','CTER3XFLAG WT', 'T34N')) %>%
    select(mutant, ORF, score) %>% 
    spread(ORF, score) %>%
    column_to_rownames(var="mutant") %>%
    as.matrix()


# cluster rows and columns
row_hc <- clustfn(t(emap))
col_hc <- clustfn(emap)

# use the first principal coordinate (classical multidimensional scaling) to flip leaves
row_pcoor <-
  cmdscale(dist(emap), eig = T, k = 1)$point %>%
  as_tibble(rownames = 'mutant') %>%
  arrange(`V1`) %>% pull(mutant)
row_hc <- rotate(row_hc, row_pcoor)

col_pcoor <-
  cmdscale(dist(t(emap)), eig = T, k = 1)$point %>%
  as_tibble(rownames = 'gene') %>%
  arrange(`V1`) %>% pull(gene)
col_hc <- rotate(col_hc, col_pcoor)

row_order <- row_hc$order
col_order <- col_hc$order


strong_mutants <- c('H141E','Y157A','D79A','D79S','T34Q', 'T34E',
                    'K101R','T34G','T34A','R108L','R108I','Q147E',
                    'H141I','G80A','H141R','R112A','R112S','R108Y',
                    'R108Q','Y148I','R108G','R108A','R78K')

# for the main figure, only show the storng mutants
# the clustering separates these mutants out as the first branch of the dendrogram
strong_dend <- as.dendrogram(row_hc)[[1]]
order.dendrogram(strong_dend) <- seq_along(get_nodes_attr(strong_dend, 'members'))
emap_strong <- emap[labels(strong_dend),]


# Make heatmap of full data, showing only strong mutants
hm <-
  Heatmap(emap_strong, name = 'GI S-score',
          col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
          # bottom_annotation = columnAnnotation(
          #   inset = colnames(emap_strong) %in% unlist(clusts),
          #   show_annotation_name = F, show_legend = F,
          #   col = list(inset = c('TRUE' = 'black', 'FALSE' = 'white')),
          #   simple_anno_size = unit(2, "mm")),
          # show_heatmap_legend = F,
          column_title = 'S. cerevisiae E-MAP deletion library, after highly correlated allele removal (n = 834)',
          row_title = 'Gsp1 point mutant',
          row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),
          cluster_rows = strong_dend, cluster_columns = as.dendrogram(col_hc),
          row_dend_reorder = FALSE, column_dend_reorder = FALSE,
          row_dend_width = unit(6, "mm"), row_dend_side = "left",
          show_column_dend = F, show_column_names = F,
          row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica')
          )

pdf('Revisions2/Fig_and_analysis_for_StatReviewer/1D_EMAP_no_dependent_alleles.pdf', height = 2.6, width = 5)
hm
dev.off()


#### Recompute the correlations between alleles and mutants
mat_GI_nhc <-
    df_GI %>% 
    filter(array_ORF %in% least_dependent_ORFs) %>% 
    pivot_wider(id_cols=query_allele_name, names_from=array_ORF, values_from=score) %>% 
    column_to_rownames('query_allele_name') %>% 
    as.matrix()


# vectorized matrix format, from https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
cor.test.p <- function(x){
    FUN <- function(x, y) cor.test(x, y, alternative = 'g', method='pearson')[["p.value"]]
    z <- outer(
      colnames(x), 
      colnames(x), 
      Vectorize(function(i,j) FUN(x[,i], x[,j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}

cor.test.cor <- function(x){
    FUN <- function(x, y) cor.test(x, y, alternative = 'g', method='pearson')[["estimate"]]
    z <- outer(
      colnames(x), 
      colnames(x), 
      Vectorize(function(i,j) FUN(x[,i], x[,j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}

mat_corr_square <- cor.test.cor(t(mat_GI_nhc))
mat_p_square <- cor.test.p(t(mat_GI_nhc))

# put alleles on rows, mutants on columns
mat_corr <-
    mat_corr_square %>% 
    as_tibble(rownames='query_1') %>% 
    pivot_longer(cols=-query_1, names_to='query_2', values_to='cor') %>% 
    filter(query_2 %in% gsp1_alleles, query_1 %in% query_alleles) %>% 
    separate(query_2, into=c(NA, 'query_2'), sep=' - ') %>% 
    pivot_wider(id_cols=query_1, names_from=query_2, values_from=cor) %>% 
    column_to_rownames('query_1') %>% 
    as.matrix()

mat_p <- 
    mat_p_square %>% 
    as_tibble(rownames='query_1') %>% 
    pivot_longer(cols=-query_1, names_to='query_2', values_to='p_raw') %>% 
    filter(query_2 %in% gsp1_alleles, query_1 %in% query_alleles) %>% 
    separate(query_2, into=c(NA, 'query_2'), sep=' - ') %>% 
    mutate('greater_fdr' = p.adjust(p_raw, method = 'fdr'),
           'greater_bonferroni' = p.adjust(p_raw, method = 'bonferroni')) %>% 
    ungroup() %>% 
    pivot_wider(id_cols=query_1, names_from=query_2, values_from=greater_fdr) %>% 
    column_to_rownames('query_1') %>% 
    as.matrix()


# load kinetics data for plotting as ratio of relative GAP, GEF efficiencies
kinetics <- read_delim('Data/kinetics_data_relative_to_WT.txt', delim='\t', col_types=cols())


##### CLUSTERING for heatmaps of mutants vs alleles

# clustering function for heatmaps
clustfn_hm <- function(mat, dist_method, clust_method) {
  
  # use euclidean clustering to get first principal coordinate ordering
  pcoor <-
    cmdscale(get_dist(mat, method = 'euclidean'), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'var') %>%
    arrange(`V1`) %>%
    pull(var)
  
  # use the distance methods and clustering method of choice to cluster
  hc <-
    mat %>%
    get_dist(method = dist_method) %>%
    hclust(method = clust_method) %>%
    rotate(pcoor)
  return(hc)
}


corr_mat <- mat_p 

# perform clustering
row.hc <- clustfn_hm(corr_mat, dist_method = 'euclidean', clust_method = 'ward.D2')
col.hc <- clustfn_hm(t(corr_mat), dist_method = 'pearson', clust_method = 'ward.D2')


##### ROTATION OF THE TREES

# rotate the yeast query gene tree so the density goes from upper left to lower right
row.hc <- rotate(row.hc, rev(row.hc$labels[row.hc$order]))

# rotate mutants based on kinetics ordering (ratio of GAP/GEF relative efficiencies)
kinetics_ordering <-
  kinetics %>%
  filter(measure == 'GAP/GEF kcat/Km') %>%
  arrange(rel_to_WT) %>%
  pull(mutant)

mutant_ordering <- colnames(corr_mat)
mutant_ordering <- mutant_ordering[order(match(mutant_ordering, kinetics_ordering))]
col.hc <- rotate(col.hc, mutant_ordering)



##### PREPARE DATA STRUCTURES FOR ANNOTATIONS

# KINETICS
# prepare dataframes of kinetic values with errors for barplots
# The ordering should be the same as the matrix columns, because the heatmap
# reorders the bars based on the dendrogram we provide
get_ordered_kinetics_df <- function(m) {
  col.hc$labels %>% 
    enframe(value = 'mutant', name = NULL) %>% 
    left_join(kinetics, by = 'mutant', ) %>%
    complete(measure, nesting(mutant), fill = list(rel_to_WT = 0, se = 0, ln_rel_to_WT = 0, ln_se = 0)) %>% 
    filter(measure == m) %>%
    mutate(mutant = factor(mutant, levels = col.hc$labels)) %>%
    arrange(mutant)
}

rel_effic <- get_ordered_kinetics_df('GAP/GEF kcat/Km')

# set color parameters for the kinetic barplots
rc_col_fn <- colorRamp2(c(-2, 0, 4.5), c(ucsf_colors$orange1, 'white', ucsf_colors$cyan1))
rc_colors <- rc_col_fn(rel_effic$ln_rel_to_WT)

# plotting heatmap
pval_heatmap <- Heatmap(
  corr_mat, show_heatmap_legend = F,
  col = colorRamp2(c(0, 0.1), c(ucsf_colors$gray1, 'white')),
  width = unit(4.3, 'cm'), gap = unit(1, "mm"),
  
  # columns (mutants)
  column_title = 'Gsp1 point mutant', column_title_gp = gpar(fontsize = 6),
  cluster_columns = as.dendrogram(col.hc), column_split = 3, column_dend_height = unit(6, "mm"),
  column_names_side = 'top', column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
  
  # mutant cluster annotations
  top_annotation = HeatmapAnnotation(
    which = 'column', height = unit(2, "mm"),
    annot = anno_block(
      labels = c('I','II','III'),
      labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica'),
      gp = gpar(fill = c('I' = ucsf_colors$navy1,
                         'II' = ucsf_colors$navy2,
                         'III' = ucsf_colors$navy3)
      ),
    )
  ),
  
  # mutant barplot annotations
  bottom_annotation = HeatmapAnnotation(
    which = 'col',
    height = unit(1.3, "cm"), 
    # gap = unit(c(4, 4), "mm"),
    annotation_height = unit(c(1.3), "cm"),

    `Ln ratio of\nGAP/GEF rel.\nefficiences` = anno_barplot(
      x = rel_effic$ln_rel_to_WT,
      ylim = c(-2.2,6.5),
      gp = gpar(fill = rc_colors, lex=0.5),
      axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'))),
    
    annotation_name_gp = gpar(fontsize = 6, fontfamily='Helvetica', lineheight = 0.9),
    annotation_name_side = 'left',
    annotation_name_rot = 0
  ),
  
  # rows (genes)
  row_title = 'Saccharomyces cerevisiae genes', row_title_gp = gpar(fontsize = 6),
  cluster_rows = as.dendrogram(row.hc), row_split = 7, row_dend_width = unit(6, "mm"),
  row_order = row.hc$order,
  show_row_names = F,
  # row_labels = label_only_gene_set_strains,
  # row_names_gp = gpar(fontsize = 2, fontfamily = 'Helvetica',
  #                     col = gene_set_label_colors),
  
  # yeast gene cluster annotation
  left_annotation = HeatmapAnnotation(
    which = 'row', width = unit(2, "mm"),
    annot = anno_block(
      labels = c('1','2','3','4','5','6','7'),
      labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica'),
      labels_rot = 0,
      gp = gpar(fill = c('1' = 'black','2' = 'black', '3' = 'black', '4' = 'black',
                         '5' = 'black','6' = 'black', '7' = 'black')
      )
    )
  ),

)

co <- unlist(column_order(pval_heatmap))
dev.off()

add_error_bars_and_stars <- function(annot, val, err) {
  decorate_annotation(annot, {
    # x positions
    err_w = 0.2  # width of error bars
    x = c(seq(1,9), c(10.54, 11.54), seq(13.08,22.08))
    x0 = c(seq(1,9), c(10.54, 11.54), seq(13.08,22.08)) - err_w
    x1 = c(seq(1,9), c(10.54, 11.54), seq(13.08,22.08)) + err_w
    
    # y positions
    y = val[co]
    y0 = val[co] - err[co]
    y1 = val[co] + err[co]
    
    # draw error bars
    grid.segments(x0=x, x1=x, y0=y0, y1=y1, default.units = 'native', gp = gpar(lex=0.5))
    grid.segments(x0=x0, x1=x1, y0=y0, y1=y0, default.units = 'native', gp = gpar(lex=0.5))
    grid.segments(x0=x0, x1=x1, y0=y1, y1=y1, default.units = 'native', gp = gpar(lex=0.5))
    
    # draw star
    x_star = c(11.54, 16.08, 20.08, 23.08)
    y_star = rep(0.6, 4)
    grid.points(x=x_star, y=y_star, pch=8, size = unit(0.75, 'mm'), gp = gpar(lex=0.5))
    
  }
  ) 
}

# save heatmap
pdf('Revisions2/Fig_and_analysis_for_StatReviewer/4A_Corr_Pvalue_Heatmap_no-hi-corr-library-genes.pdf', width = 2.7, height = 5)
draw(pval_heatmap)

add_error_bars_and_stars('Ln ratio of\nGAP/GEF rel.\nefficiences',
                         rel_effic$ln_rel_to_WT, rel_effic$ln_se)

dev.off()




### also remake the E-MAP heatmap for 1d