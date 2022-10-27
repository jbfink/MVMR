# /genetics/PAREG/perrotn/automated_scripts/ADD_FEATURES_MVMR_RESULTS_FILE.r

rm(list = ls())

library("plyr")
library("optparse")

option_list <- list(
  make_option(c("--file"),
    type = "character", default = NULL,
    help = "MVMR results file path", metavar = "character"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file))
{
  message("--file should be define")
  exit()
}


########################## TEST #########################
# opt=NULL
# opt$file="/storage/genetics_work2/perrotn/PURE/MVMR_BM_consortia/NO_MHC_MISSENSE_SPLICING/CAD_cardiogram_c4d_ukb_2017/5_MR_RESULTS/LD_0.1/MR_RESULTS.LD_PRUNED_0.1_EUROPEAN_CMET_F11_F11_CIS_200000_PVALUE_0.01_F11_EUROPEAN_CMET_CCL14_CCL14_vs_CAD_cardiogram_c4d_ukb_2017.txt"
#########################################################


# ADD TO MR RESULTS:
# KEEP F-stat >=10
# PRIO MR METHOD (IVW OR EGGER)
# PRIO P-value

MR_RESULTS=data.table::fread(opt$file, data.table=F, stringsAsFactors=F)

TO_SKIP=NA
ALL_RESULTS=NULL

EXP=""
OUT=""
for (EXP in unique(MR_RESULTS$id.exposure))
{
    
    if (!is.na(TO_SKIP))
    {
        if (EXP == TO_SKIP)
        {
        ALL_RESULTS = rbind(ALL_RESULTS, MR_RESULTS[which(MR_RESULTS$id.exposure == EXP),])
        TO_SKIP=NA
        next()
        }
    }
        
    for (OUT in unique(MR_RESULTS$id.outcome))
    {
        
        if (length(which(MR_RESULTS$id.exposure == EXP & MR_RESULTS$id.outcome == OUT & MR_RESULTS$method %in% c("Inverse variance weighted", "MR Egger", "Multivariable inverse-variance weighted", "Multivariable MR-Egger")))==0)
        {
            next()
        }
        
        DATA=MR_RESULTS[which(MR_RESULTS$id.exposure == EXP & MR_RESULTS$id.outcome == OUT & MR_RESULTS$method %in% c("Inverse variance weighted", "MR Egger", "Multivariable inverse-variance weighted", "Multivariable MR-Egger")),]

        if (!is.na(unique(DATA$F_STAT)))
        {
            if((unique(na.omit(DATA$F_STAT)))<10)
            {
                TO_SKIP=EXP
                next()
            }
        }

        PRIO=DATA[1,]
        PRIO[,which(colnames(PRIO) %in% c("method","b", "se", "pval", "Q", "Q_df", "Q_pval", "egger_intercept", "egger_intercept_se", "egger_intercept_pval"))]=NA

        if (DATA$egger_intercept_pval[which(DATA$method == unique(DATA$method)[grep("Egger", unique(DATA$method))])]<0.05)
        {
            PRIO$method=paste0("MR METHOD PRIORITIZED: ", unique(DATA$method)[grep("Egger", unique(DATA$method))])
            PRIO$b=DATA$b[which(DATA$method == unique(DATA$method)[grep("Egger", unique(DATA$method))])]
            PRIO$se=DATA$se[which(DATA$method == unique(DATA$method)[grep("Egger", unique(DATA$method))])]
            PRIO$pval=DATA$pval[which(DATA$method == unique(DATA$method)[grep("Egger", unique(DATA$method))])]
        } else {
            PRIO$method=paste0("MR METHOD PRIORITIZED: ", unique(DATA$method)[grep("variance", unique(DATA$method))])
            PRIO$b=DATA$b[which(DATA$method == unique(DATA$method)[grep("variance", unique(DATA$method))])]
            PRIO$se=DATA$se[which(DATA$method == unique(DATA$method)[grep("variance", unique(DATA$method))])]
            PRIO$pval=DATA$pval[which(DATA$method == unique(DATA$method)[grep("variance", unique(DATA$method))])]
        }

        ALL_RESULTS = rbind(ALL_RESULTS, rbind(MR_RESULTS[which(MR_RESULTS$id.exposure == EXP & MR_RESULTS$id.outcome == OUT),],PRIO))

    } # END OUT LOOP

} # END EXP LOOP



# SAVE _PRIO_MR_METHOD.txt 
data.table::fwrite(x=ALL_RESULTS,file=gsub(opt$file, pattern = ".txt", replacement = "_PRIO_MR_METHOD.txt"), sep="\t", row.names=F, col.names=T, quote=F)



