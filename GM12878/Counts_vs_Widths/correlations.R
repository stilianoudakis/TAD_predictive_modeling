ab <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/ab.rds")
dgv <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/dgv.rds")
gerp <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/gerp.rds")
vmr <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/vmr.rds")
se <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/se.rds")
broad <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/broad.rds")
combined <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/combined.rds")
dnase <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/dnase.rds")
histone <- readRDS("/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/histone.rds")

gm12878_5kb <- cbind.data.frame(ab,
								dgv[,-1],
								gerp[,-1],
								vmr[,-1],
								se[,-1],
								broad[,-1],
								combined[,-1],
								dnase[,-1],
								histone[,-1])
								
dim(gm12878_5kb)

saveRDS(gm12878_5kb, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/gm12878_5kb.rds")


toMatch <- c("dist", "count")

matches <- unique (grep(paste(toMatch,collapse="|"), 
                        names(gm12878_5kb), value=TRUE))
matches <- c("y",matches)



#Calculating correlations between width and count overlap variables
cordf <- data.frame(Feature = names(gm12878_5kb[,-which(names(gm12878_5kb) %in% matches)]),
					Correlation = c(cor(gm12878_5kb$A,gm12878_5kb$A_count),
									cor(gm12878_5kb$B,gm12878_5kb$B_count),
									cor(gm12878_5kb$complex,gm12878_5kb$complex_count),
									cor(gm12878_5kb$deletion,gm12878_5kb$deletion_count),
									cor(gm12878_5kb$duplication,gm12878_5kb$duplication_count),
									cor(gm12878_5kb$gain_loss,gm12878_5kb$gain_loss_count),
									cor(gm12878_5kb$insertion,gm12878_5kb$insertion_count),
									cor(gm12878_5kb$inversion,gm12878_5kb$inversion_count),
									cor(gm12878_5kb$mobile_element_insertion,gm12878_5kb$mobile_element_insertion_count),
									cor(gm12878_5kb$novel_sequence_insertion,gm12878_5kb$novel_sequence_insertion_count),
									cor(gm12878_5kb$sequence_alteration,gm12878_5kb$sequence_alteration_count),
									cor(gm12878_5kb$tandem_duplication,gm12878_5kb$tandem_duplication_count),
									cor(gm12878_5kb$GERP,gm12878_5kb$GERP_count),
									cor(gm12878_5kb$GM12878_se,gm12878_5kb$GM12878_se_count),
									cor(gm12878_5kb$VMR,gm12878_5kb$VMR_count),
									cor(gm12878_5kb$`Gm12878-TxnElongation`,gm12878_5kb$`Gm12878-TxnElongation_count`),
									cor(gm12878_5kb$`Gm12878-WeakTxn`,gm12878_5kb$`Gm12878-WeakTxn_count`),
									cor(gm12878_5kb$`Gm12878-Repressed`,gm12878_5kb$`Gm12878-Repressed_count`),
									cor(gm12878_5kb$`Gm12878-Heterochromlo`,gm12878_5kb$`Gm12878-Heterochromlo_count`),
									cor(gm12878_5kb$`Gm12878-RepetitiveCNV`,gm12878_5kb$`Gm12878-RepetitiveCNV_count`),
									cor(gm12878_5kb$`Gm12878-ActivePromoter`,gm12878_5kb$`Gm12878-ActivePromoter_count`),
									cor(gm12878_5kb$`Gm12878-WeakPromoter`,gm12878_5kb$`Gm12878-WeakPromoter_count`),
									cor(gm12878_5kb$`Gm12878-PoisedPromoter`,gm12878_5kb$`Gm12878-PoisedPromoter_count`),
									cor(gm12878_5kb$`Gm12878-StrongEnhancer`,gm12878_5kb$`Gm12878-StrongEnhancer_count`),
									cor(gm12878_5kb$`Gm12878-WeakEnhancer`,gm12878_5kb$`Gm12878-WeakEnhancer_count`),
									cor(gm12878_5kb$`Gm12878-Insulator`,gm12878_5kb$`Gm12878-Insulator_count`),
									cor(gm12878_5kb$`Gm12878-TxnTransition`,gm12878_5kb$`Gm12878-TxnTransition_count`),
									cor(gm12878_5kb$`Gm12878-CTCF`,gm12878_5kb$`Gm12878-CTCF_count`),
									cor(gm12878_5kb$`Gm12878-E`,gm12878_5kb$`Gm12878-E_count`),
									cor(gm12878_5kb$`Gm12878-PF`,gm12878_5kb$`Gm12878-PF_count`),
									cor(gm12878_5kb$`Gm12878-R`,gm12878_5kb$`Gm12878-R_count`),
									cor(gm12878_5kb$`Gm12878-T`,gm12878_5kb$`Gm12878-T_count`),
									cor(gm12878_5kb$`Gm12878-TSS`,gm12878_5kb$`Gm12878-TSS_count`),
									cor(gm12878_5kb$`Gm12878-WE`,gm12878_5kb$`Gm12878-WE_count`),
									cor(gm12878_5kb$`Gm12878-DNase`,gm12878_5kb$`Gm12878-DNase_count`),
									cor(gm12878_5kb$`Gm12878-H2az`,gm12878_5kb$`Gm12878-H2az_count`),
									cor(gm12878_5kb$`Gm12878-H3k27ac`,gm12878_5kb$`Gm12878-H3k27ac_count`),
									cor(gm12878_5kb$`Gm12878-H3k27me3`,gm12878_5kb$`Gm12878-H3k27me3_count`),
									cor(gm12878_5kb$`Gm12878-H3k36me3`,gm12878_5kb$`Gm12878-H3k36me3_count`),
									cor(gm12878_5kb$`Gm12878-H3k4me1`,gm12878_5kb$`Gm12878-H3k4me1_count`),
									cor(gm12878_5kb$`Gm12878-H3k4me2`,gm12878_5kb$`Gm12878-H3k4me2_count`),
									cor(gm12878_5kb$`Gm12878-H3k4me3`,gm12878_5kb$`Gm12878-H3k4me3_count`),
									cor(gm12878_5kb$`Gm12878-H3k79me2`,gm12878_5kb$`Gm12878-H3k79me2_count`),
									cor(gm12878_5kb$`Gm12878-H3k9ac`,gm12878_5kb$`Gm12878-H3k9ac_count`),
									cor(gm12878_5kb$`Gm12878-H3k9me3`,gm12878_5kb$`Gm12878-H3k9me3_count`),
									cor(gm12878_5kb$`Gm12878-H4k20me1`,gm12878_5kb$`Gm12878-H4k20me1_count`)))

saveRDS(cordf, "/home/stilianoudakisc/TAD_data_analysis/model_filtering/count_vs_width/cordf.rds")