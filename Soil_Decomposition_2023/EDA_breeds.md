---
title: "EDA_breeds"
author: "Aleksei Zverev"
date: "2023-06-28"
output: 
  html_document: 
    toc: yes
    toc_float: yes
    toc_collapsed: yes
    toc_depth: 4
    number_sections: yes
    keep_md: yes
    code_folding: hide
    theme: lumen

---



# Soils and substrates


```r
format_data_for_intersection <- function(ps, group){
  physeq <- prune_taxa(taxa_sums(ps) > 0, ps)
  groups <- levels(sample_data(physeq)[[group]] %>% as.factor())
  data <- merge_samples(physeq, group) %>% 
    psmelt() %>% 
    group_by(Sample, OTU) %>% 
    summarise(ASVs_abund = list(paste(OTU, 1:sum(Abundance))), Abund = sum(Abundance), .groups='keep') %>% 
    filter(Abund > 0)
  
  asvs <- data %>% select(Sample, OTU) %>% group_by(Sample) %>% summarise(ASVs = list(OTU), .groups='keep') %>% as.list()
  d1 <- asvs[[2]]
  names(d1) <- asvs[[1]]
  d1
  
  weighted.asvs <- data %>% select(Sample, ASVs_abund) %>% group_by(Sample) %>% summarise(ASVs = list(unlist(ASVs_abund)))
  d2 <- weighted.asvs[[2]]
  names(d2) <- weighted.asvs[[1]]
  d2
  
  list(d1, d2)
}

vienn_from_intersects <- function(data){
  require(ggVennDiagram)
  
  list(ggVennDiagram(data[[1]]) + ggtitle("ASVs") + scale_fill_distiller(palette = "OrRd", trans = "reverse"),
       ggVennDiagram(data[[2]]) + ggtitle("Reads") + scale_fill_distiller(palette = "OrRd", trans = "reverse"))
}

find_intersections <- function(c0, c1, c2){
  require(dplyr)
  c0_unique <- setdiff(c0, c(c1, c2))
  intersect <- intersect(c0, intersect(c1, c2))
  c0_c1 <- setdiff(intersect(c0, c1), c2)
  c0_c2 <- setdiff(intersect(c0, c2), c1)
  c(c0_unique %>% length(),
    c0_c1 %>% length(),
    intersect %>% length(),
    c0_c2 %>% length())
}
```

## Soils


```r
p <- vienn_from_intersects(
  format_data_for_intersection(ps %>% subset_samples(
    Substrate %in% c("DF", "PS", "BG")), "Substrate"
    )
  )

ggarrange(p[[1]], p[[2]])
```

![](EDA_breeds_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

## Substrates


```r
p <- vienn_from_intersects(
  format_data_for_intersection(
    ps %>% subset_samples((Substrate %in% c("Cannb", "Wheat", "Oat")) & (Type == "substrate")), "Substrate"
    )
  )

ggarrange(p[[1]], p[[2]])
```

![](EDA_breeds_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

# Breeds

## Soils, Substrates and breeds

What share breeds takes from substrate and from soil?


```r
intersects_from_groups <- function(ps.int, groups){

  d <- format_data_for_intersection(ps.int, "Compare")
  list(
    qual = find_intersections(d[[1]][[paste0(groups[1], ".")]],
                              d[[1]][[paste0(".", groups[2])]],
                              d[[1]][[paste0(".", groups[3])]]),
    quan = find_intersections(d[[2]][[paste0(groups[1], ".")]], 
                              d[[2]][[paste0(".", groups[2])]], 
                              d[[2]][[paste0(".", groups[3])]])
    )
}


PS_Oat <- intersects_from_groups(
  ps %>% subset_samples((Substrate %in% c("PS", "Oat", "")) & 
                          (Primer %in% c("PS_Oat", "")) & 
                          (Type != "use")),
                       c("PS_Oat", "PS", "Oat"))

DF_Oat <- intersects_from_groups(
  ps %>% subset_samples((Substrate %in% c("DF", "Oat", "")) & 
                          (Primer %in% c("DF_Oat", "")) & 
                          (Type != "use")),
                       c("DF_Oat", "DF", "Oat"))

BG_Oat <- intersects_from_groups(
  ps %>% subset_samples((Substrate %in% c("BG", "Oat", "")) & 
                          (Primer %in% c("BG_Oat", "")) & 
                          (Type != "use")),
                       c("BG_Oat", "BG", "Oat"))

PS_Cannb <- intersects_from_groups(
  ps %>% subset_samples((Substrate %in% c("PS", "Cannb", "")) & 
                          (Primer %in% c("PS_Cannb", "")) & 
                          (Type != "use")),
                       c("PS_Cannb", "PS", "Cannb"))

DF_Cannb <- intersects_from_groups(
  ps %>% subset_samples((Substrate %in% c("DF", "Cannb", "")) & 
                          (Primer %in% c("DF_Cannb", "")) & 
                          (Type != "use")),
                       c("DF_Cannb", "DF", "Cannb"))

BG_Cannb <- intersects_from_groups(
  ps %>% subset_samples((Substrate %in% c("BG", "Cannb", "")) & 
                          (Primer %in% c("BG_Cannb", "")) & 
                          (Type != "use")),
                       c("BG_Cannb", "BG", "Cannb"))
```



```r
intersect.qual <- data.frame(PS_Oat$qual, PS_Cannb$qual,
                             DF_Oat$qual, DF_Cannb$qual,
                             BG_Oat$qual, BG_Cannb$qual) %>% t() %>% data.frame()
colnames(intersect.qual) <- c("Unique", "From soil", "From substrate", "Common")
intersect.qual$Name <- gsub(".qual", "", rownames(intersect.qual))
intersect.qual.long <- intersect.qual %>% tidyr::pivot_longer(!Name, names_to = "group", values_to = "count")

intersect.quan <- data.frame(PS_Oat$quan, PS_Cannb$quan,
                             DF_Oat$quan, DF_Cannb$quan,
                             BG_Oat$quan, BG_Cannb$quan) %>% t() %>% data.frame()
colnames(intersect.quan) <- c("Unique", "From soil", "From substrate", "Common")
intersect.quan$Name <- gsub(".quan", "", rownames(intersect.quan))
intersect.quan.long <- intersect.quan %>% tidyr::pivot_longer(!Name, names_to = "group", values_to = "count")



p1 <- ggplot(intersect.qual.long, aes(Name, count, fill=group)) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  theme_light() + ggtitle('ASVs')
p2 <- ggplot(intersect.quan.long, aes(Name, count, fill=group)) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  theme_light() + ggtitle('Reads')
ggarrange(p1, p2, nrow=2, common.legend = T)
```

![](EDA_breeds_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

## Soils, Substrates and primers for the major taxa (more than 1%)


```r
ps.2 <- filter_taxa(ps, function (x) {sum(x) > 60}, prune=TRUE)


PS_Oat <- intersects_from_groups(
  ps.2 %>% subset_samples((Substrate %in% c("PS", "Oat", "")) & 
                          (Primer %in% c("PS_Oat", "")) & 
                          (Type != "use")),
                       c("PS_Oat", "PS", "Oat"))

DF_Oat <- intersects_from_groups(
  ps.2 %>% subset_samples((Substrate %in% c("DF", "Oat", "")) & 
                          (Primer %in% c("DF_Oat", "")) & 
                          (Type != "use")),
                       c("DF_Oat", "DF", "Oat"))

BG_Oat <- intersects_from_groups(
  ps.2 %>% subset_samples((Substrate %in% c("BG", "Oat", "")) & 
                          (Primer %in% c("BG_Oat", "")) & 
                          (Type != "use")),
                       c("BG_Oat", "BG", "Oat"))

PS_Cannb <- intersects_from_groups(
  ps.2 %>% subset_samples((Substrate %in% c("PS", "Cannb", "")) & 
                          (Primer %in% c("PS_Cannb", "")) & 
                          (Type != "use")),
                       c("PS_Cannb", "PS", "Cannb"))

DF_Cannb <- intersects_from_groups(
  ps.2 %>% subset_samples((Substrate %in% c("DF", "Cannb", "")) & 
                          (Primer %in% c("DF_Cannb", "")) & 
                          (Type != "use")),
                       c("DF_Cannb", "DF", "Cannb"))

BG_Cannb <- intersects_from_groups(
  ps.2 %>% subset_samples((Substrate %in% c("BG", "Cannb", "")) & 
                          (Primer %in% c("BG_Cannb", "")) & 
                          (Type != "use")),
                       c("BG_Cannb", "BG", "Cannb"))
```



```r
intersect.qual <- data.frame(PS_Oat$qual, PS_Cannb$qual,
                             DF_Oat$qual, DF_Cannb$qual,
                             BG_Oat$qual, BG_Cannb$qual) %>% t() %>% data.frame()
colnames(intersect.qual) <- c("Unique", "From soil", "From substrate", "Common")
intersect.qual$Name <- gsub(".qual", "", rownames(intersect.qual))
intersect.qual.long <- intersect.qual %>% tidyr::pivot_longer(!Name, names_to = "group", values_to = "count")

intersect.quan <- data.frame(PS_Oat$quan, PS_Cannb$quan,
                             DF_Oat$quan, DF_Cannb$quan,
                             BG_Oat$quan, BG_Cannb$quan) %>% t() %>% data.frame()
colnames(intersect.quan) <- c("Unique", "From soil", "From substrate", "Common")
intersect.quan$Name <- gsub(".quan", "", rownames(intersect.quan))
intersect.quan.long <- intersect.quan %>% tidyr::pivot_longer(!Name, names_to = "group", values_to = "count")



p1 <- ggplot(intersect.qual.long, aes(Name, count, fill=group)) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  theme_light() + ggtitle('ASVs')
p2 <- ggplot(intersect.quan.long, aes(Name, count, fill=group)) + 
  geom_bar(aes(), stat="identity", position="stack") + 
  theme_light() + ggtitle('Reads')
ggarrange(p1, p2, nrow=2, common.legend = T)
```

![](EDA_breeds_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

As a result: most parts of the taxa are unique for the breeds, or common with soils. This differences are not in minor taxa, but in both major and minor. Soils themselves have small share of common reads, so, it's effect is directly the result of training of the communities.

## Interserction of breeds

### Oat breeds


```r
# p <- vienn_from_intersects(
#   format_data_for_intersection(ps %>% subset_samples(
#     (Primer %in% c("DF_Oat", "PS_Oat", "BG_Oat")) & (Type == "breed")), "Primer"
#     )
#   )
# 
# ggarrange(p[[1]], p[[2]])


p <- vienn_from_intersects(
  format_data_for_intersection(ps.2 %>% subset_samples(
    (Primer %in% c("DF_Oat", "PS_Oat", "BG_Oat")) & (Type == "breed")), "Primer"
    )
  )

ggarrange(p[[1]], p[[2]])
```

![](EDA_breeds_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

### Cannb breeds


```r
# p <- vienn_from_intersects(
#   format_data_for_intersection(ps %>% subset_samples(
#     (Primer %in% c("DF_Cannb", "PS_Cannb", "BG_Cannb")) & (Type == "breed")), "Primer"
#     )
#   )
# 
# ggarrange(p[[1]], p[[2]])

p <- vienn_from_intersects(
  format_data_for_intersection(ps.2 %>% subset_samples(
    (Primer %in% c("DF_Cannb", "PS_Cannb", "BG_Cannb")) & (Type == "breed")), "Primer"
    )
  )

ggarrange(p[[1]], p[[2]])
```

![](EDA_breeds_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

Breeds themselves are different (alongside a bit more similarity between PS and DF-based)

## Soils, substrates and breeds on beta-diversity


```r
beta_plot <- function(ps, method, distance, ...){
  require(phyloseq)
  require(ggplot2)
  
  ps.prop <- transform_sample_counts(ps, function(x) x/sum(x))
  ords <- ordinate(ps.prop, method=method, distance=distance)
  plot_ordination(ps.prop, ords, title=deparse(substitute(ps)), ...) +
    geom_point(size=3, alpha=0.7) + 
    theme_light()
}


# beta_plot(ps %>% subset_samples(Type != "use"), "PCoA", "bray", color="Compare", shape = "Type")
beta_plot(ps.2 %>% subset_samples(Type != "use"), "PCoA", "bray", color="Compare", shape = "Type")
```

![](EDA_breeds_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

BAGS and their breeds clearly forms own cluster, and soils are from the other one. In terms of taxa, what is the difference?

## Taxa differences between BG - DF&PS clusters


```r
sig_table <- function(ps, formula, threshold){
  require(DESeq2)
  require(dplyr)
  
  ds <- phyloseq_to_deseq2(ps, formula)
  ds = estimateSizeFactors(ds, type="poscounts")
  ds = estimateDispersions(ds, fitType = "local")
  ds = DESeq(ds)
  #mcols(ds, use.names=TRUE)
  res = results(ds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab <- sigtab %>% data.frame() %>% dplyr::filter(baseMean >= threshold[1],
                                                      log2FoldChange >= threshold[2] | log2FoldChange <= threshold[2])
  if (nrow(sigtab) == 0) {
    return(NA)
  }
  sigtab = cbind(as(sigtab, "data.frame"),
                 as(tax_table(ps)[rownames(sigtab), ], "matrix")
  )
  return(sigtab)
}

draw_sig_table <- function(sig_table, rank){
  require(ggplot2)
  
  sig_table$Class <-  ifelse(is.na(sig_table$Class), 
                             paste(sig_table$Phylum, "// NA"), 
                             paste(sig_table$Class))
  sig_table$Order <-  ifelse(is.na(sig_table$Order), 
                             paste(sig_table$Class, "// NA"), 
                             paste(sig_table$Order))
  sig_table$Family <- ifelse(is.na(sig_table$Family), 
                             paste(sig_table$Order, "// NA"), 
                             paste(sig_table$Family))
  sig_table$Genus <- ifelse(is.na(sig_table$Genus), 
                            paste(sig_table$Class, "// NA"), 
                            paste(sig_table$Genus))
  
  sig_table[sig_table == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
  sig_table[sig_table == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Pararhizobium"
  
  ggplot(sig_table, aes(y = sig_table[,rank], x=log2FoldChange, size = log(baseMean))) + 
    geom_point(aes(color=sig_table[,rank])) + guides(colour=FALSE) +
    theme_light() + facet_grid(rows = vars(Phylum), space = 'free_y', scales = 'free') +
    theme(strip.text.y = element_text(angle=0)) +
    theme(legend.position = "top") + ylab(rank)
}
```

### BG - DF


```r
ps.pair.bg.df <- ps.2 %>% subset_samples(
  ((Substrate %in% c("BG", "DF")) | ((Primer %in% c("BG_Cannb", "DF_Cannb", "BG_Oat", "DF_Oat")) & (Type == "breed")))
                                         )

ps.pair.bg.df@sam_data$Compare <- gsub("_Oat.", "", ps.pair.bg.df@sam_data$Compare)
ps.pair.bg.df@sam_data$Compare <- gsub("_Cannb.", "", ps.pair.bg.df@sam_data$Compare)
ps.pair.bg.df@sam_data$Compare <- gsub("\\.", "", ps.pair.bg.df@sam_data$Compare)
ps.pair.bg.df@sam_data$Compare <- as.factor(ps.pair.bg.df@sam_data$Compare)

sig.tabe.pair.bg.df <- sig_table(ps.pair.bg.df, ~Compare, c(30, 2.5))

draw_sig_table(sig.tabe.pair.bg.df, "Genus")
```

![](EDA_breeds_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

### BG - PS


```r
ps.pair.bg.ps <- ps.2 %>% subset_samples(
  ((Substrate %in% c("BG", "PS")) | ((Primer %in% c("BG_Cannb", "PS_Cannb", "BG_Oat", "PS_Oat")) & (Type == "breed")))
                                         )

ps.pair.bg.ps@sam_data$Compare <- gsub("_Cannb.", "", ps.pair.bg.ps@sam_data$Compare)
ps.pair.bg.ps@sam_data$Compare <- gsub("_Oat.", "", ps.pair.bg.ps@sam_data$Compare)
ps.pair.bg.ps@sam_data$Compare <- gsub("\\.", "", ps.pair.bg.ps@sam_data$Compare)
ps.pair.bg.ps@sam_data$Compare <- as.factor(ps.pair.bg.ps@sam_data$Compare)

sig.tabe.pair.bg.ps <- sig_table(ps.pair.bg.ps, ~Compare, c(30, 2.5))

draw_sig_table(sig.tabe.pair.bg.ps, "Genus")
```

![](EDA_breeds_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
tax.names <- c(rownames(sig.tabe.pair.bg.df), rownames(sig.tabe.pair.bg.ps)) %>% unique()
```

### Combined data


```r
plot_heatmap <- function(sig.ps, group = "SampleID", log.transform = TRUE){
  require(dplyr)
  require(phyloseq)
  require(ggplot2)
  
  sig.taxa.long <- psmelt(sig.ps) %>%
    arrange(Phylum) %>% 
    mutate(row = row_number()) %>% 
    mutate_at(c('Genus', 'Family', 'Order', 'Class', 'Phylum'), as.character)  # fuckin magic for ifelse tier
  
  sig.taxa.long$Abundance <- as.numeric(sig.taxa.long$Abundance)
  sig.taxa.long$Taxa <- ifelse(is.na(sig.taxa.long$Genus),
                               ifelse(is.na(sig.taxa.long$Family), 
                                      ifelse(is.na(sig.taxa.long$Order), 
                                             ifelse(is.na(sig.taxa.long$Class), 
                                                    paste(sig.taxa.long$Phylum, "// NA"), 
                                                    paste(sig.taxa.long$Class, "// NA")),
                                             paste(sig.taxa.long$Order, "// NA")),
                                      paste(sig.taxa.long$Family, "// NA")), 
                               sig.taxa.long$Genus)
  sig.taxa.long[sig.taxa.long == "Burkholderia-Caballeronia-Paraburkholderia"] <- "Burkholderia"
  sig.taxa.long[sig.taxa.long == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Pararhizobium"
  
  
  ggplot(sig.taxa.long, aes(x = sig.taxa.long[,group], y = reorder(Taxa, row))) + 
    {if(log.transform)geom_tile(aes(fill=log(Abundance)))} +
    {if(!log.transform)geom_tile(aes(fill=Abundance))} +
    scale_fill_distiller(palette = "OrRd", trans = "reverse") +
    facet_grid(rows = vars(Phylum), scales = "free_y", space = "free") +
    theme(strip.text.y = element_text(angle = 0),
          panel.spacing = unit(0.05,'lines')) +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90))
}


ps.plot <- prune_taxa(tax.names, ps.2 %>% subset_samples(Type %in% c("breed", "soil")))
ps.plot@sam_data$Compare <- gsub("\\.", "", ps.plot@sam_data$Compare)
ps.plot@sam_data$Facet <- substr(ps.plot@sam_data$Compare, 1, 2)

plot_heatmap(ps.plot, "Compare") + facet_grid(~Facet, scales = "free")
```

![](EDA_breeds_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

# Usage

## Beta-diversity


```r
beta_plot(ps.2 %>% subset_samples(Type %in% c("breed", "use")), "PCoA", "bray", color="Primer", shape = "Substrate")
```

![](EDA_breeds_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


```r
beta_plot(ps.use, "PCoA", "bray", color="Primer", shape = "Substrate")
```

![](EDA_breeds_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

```r
dist <- phyloseq::distance(ps.use, method = "bray")
sample_df <- data.frame(sample_data(ps.use))
 
permanova <- adonis2(dist ~ Primer + Substrate, data = sample_df)
permanova
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = dist ~ Primer + Substrate, data = sample_df)
##            Df SumOfSqs      R2      F Pr(>F)    
## Primer      6   10.670 0.25657  9.634  0.001 ***
## Substrate   2    9.689 0.23299 26.246  0.001 ***
## Residual  115   21.228 0.51044                  
## Total     123   41.588 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Nutritions and their correlations


```r
corrdata <- ps.use@sam_data %>% 
  data.frame() %>% 
  select(-Type, -Replica, -Nitrates) %>% 
  mutate(across(c(-Primer, -Substrate), as.double))

# corrdata

aov(Decomposition ~ Primer + Substrate, corrdata) %>% summary()
```

```
##              Df Sum Sq Mean Sq F value Pr(>F)    
## Primer        6   4009     668   29.36 <2e-16 ***
## Substrate     2  14514    7257  318.94 <2e-16 ***
## Residuals   115   2617      23                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
cor <- Hmisc::rcorr(as.matrix(corrdata %>% select(-Primer, -Substrate)))
cor$r
```

```
##               Decomposition   Humidity AshContent         рН          C
## Decomposition     1.0000000 -0.8517039 -0.6564786 -0.5338455 -0.7566614
## Humidity         -0.8517039  1.0000000  0.3916799  0.3513446  0.6481284
## AshContent       -0.6564786  0.3916799  1.0000000  0.2896756  0.3560688
## рН               -0.5338455  0.3513446  0.2896756  1.0000000  0.8282271
## C                -0.7566614  0.6481284  0.3560688  0.8282271  1.0000000
```


```r
veganifyOTU <- function(physeq){
  require(phyloseq)
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

ps.top1k <- prune_taxa(names(sort(taxa_sums(ps.use), TRUE)[1:1000]), ps.use)
X <- veganifyOTU(ps.top1k)
  
vare.cca <- vegan::cca(X ~ Decomposition + Humidity + AshContent + рН + C, data=corrdata)
anova(vare.cca)
```

```
## Permutation test for cca under reduced model
## Permutation: free
## Number of permutations: 999
## 
## Model: cca(formula = X ~ Decomposition + Humidity + AshContent + рН + C, data = corrdata)
##           Df ChiSquare      F Pr(>F)    
## Model      5    1.5019 6.2305  0.001 ***
## Residual 118    5.6890                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anova(vare.cca, by="terms")
```

```
## Permutation test for cca under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## Model: cca(formula = X ~ Decomposition + Humidity + AshContent + рН + C, data = corrdata)
##                Df ChiSquare       F Pr(>F)    
## Decomposition   1    0.4345  9.0120  0.001 ***
## Humidity        1    0.2028  4.2057  0.001 ***
## AshContent      1    0.2090  4.3357  0.001 ***
## рН              1    0.4953 10.2743  0.001 ***
## C               1    0.1603  3.3246  0.001 ***
## Residual      118    5.6890                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## CCA


```r
species.data <- vare.cca$CCA$v %>% 
               data.frame() %>% 
               mutate(ASV = rownames(.)) %>% 
               inner_join(data.frame(ASV = names(taxa_sums(ps.top1k)),
                                     Total.abund = taxa_sums(ps.top1k),
                                     ps.top1k@tax_table[,2], # Phylum
                                     ps.top1k@tax_table[,3], # Class
                                     ps.top1k@tax_table[,4],
                                     ps.top1k@tax_table[,5]),
                          by = "ASV")

samples.data <- vare.cca$CCA$u %>% 
  data.frame() %>% 
  mutate(Names = rownames(.)) %>% 
  inner_join(ps.top1k@sam_data %>% 
               data.frame() %>% 
               mutate(Samples = rownames(.)), by = c("Names" = "Samples"))
```


```r
major.phyla <- species.data %>% 
  group_by(Phylum) %>% 
  summarize(sum = sum(Total.abund)) %>% 
  arrange(desc(sum)) %>% 
  select(Phylum) %>% 
  head(10) %>% 
  as.vector()


ggplot() +
  geom_point(data=samples.data,
             aes(x=CCA1, y=CCA2, color=Primer, shape = Substrate), size=3, alpha=0.7) +
  geom_segment(data = vare.cca$CCA$biplot %>% data.frame(), 
               aes(x = 0, xend = CCA1, y = 0, yend = CCA2), 
               alpha=0.8, color = "black",arrow = arrow(angle = 3)) +
  geom_text(data = vare.cca$CCA$biplot %>% 
                    data.frame() %>% 
                    mutate(Label = rownames(.)), 
            aes(x=CCA1, y=CCA2, label= Label,
                hjust = -0.5), size=4) +
  theme_light() +
  ggtitle("B. Samples")
```

![](EDA_breeds_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
for (i in major.phyla$Phylum) {
  p <- ggplot() +
    geom_point(data=species.data,
               aes(x=CCA1, y=CCA2, size=Total.abund), alpha=0.2, color="grey80") +
    geom_point(data=species.data %>% filter(Phylum == i),
               aes(x=CCA1, y=CCA2, color=Order, size=Total.abund), alpha=0.9) +
    geom_segment(data = vare.cca$CCA$biplot %>% data.frame(), 
                 aes(x = 0, xend = CCA1, y = 0, yend = CCA2), 
                 alpha=0.8, color = "black",arrow = arrow(angle = 3)) +
    geom_text(data = vare.cca$CCA$biplot %>% 
                      data.frame() %>% 
                      mutate(Label = rownames(.)), 
              aes(x=CCA1, y=CCA2, label= Label,
                  hjust = -0.5), size=4) +
    theme_light() +
    ggtitle(i)
  print(p)
}
```

![](EDA_breeds_files/figure-html/unnamed-chunk-20-2.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-3.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-4.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-5.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-6.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-7.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-8.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-9.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-10.png)<!-- -->![](EDA_breeds_files/figure-html/unnamed-chunk-20-11.png)<!-- -->

