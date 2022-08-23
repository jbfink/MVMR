# /genetics/PAREG/perrotn/automated_scripts/MVMR.r

rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)


EXPOSURES_INPUT <- as.character(args[1]) # INPUT_DIR FOR INPUT_FILE
OUTCOMES_INPUT <- as.character(args[2])
OUTPUT_FILE <- as.character(args[3])
POPULATION <- as.character(args[4])
PANEL <- as.character(args[5])
ASSAY <- as.character(args[6])
GENE <- as.character(args[7])


library("TwoSampleMR")
library("ggplot2")
library("tidyr")
library("psych")
library("mr.raps")
library("plyr")
library("MendelianRandomization")

options(warn = 1)

########################## TEST ##########################
# POPULATION <- "EUROPEAN"
# PANEL <- "CMET"
# ASSAY <- "F11"
# GENE <- "F11"
# EXPOSURES_INPUT <- "/storage/genetics_work2/perrotn/PURE/MVMR_consortia/NO_MHC_MISSENSE_SPLICING/PARENTAL_LIFESPAN_LIFEGEN_2017/5_MR_RESULTS/LD_0.1/EXPOSURES.LD_PRUNED_0.1_EUROPEAN_CMET_F11_F11_CIS_200000_PVALUE_0.01_6_chr_pos"
# OUTCOMES_INPUT <- "/storage/genetics_work2/perrotn/PURE/MVMR_consortia/NO_MHC_MISSENSE_SPLICING/PARENTAL_LIFESPAN_LIFEGEN_2017/5_MR_RESULTS/LD_0.1/OUTCOME.LD_PRUNED_0.1_EUROPEAN_CMET_F11_F11_CIS_200000_PVALUE_0.01_6_chr_pos"
# OUTPUT_FILE <- "/storage/genetics_work2/perrotn/PURE/MVMR_consortia/NO_MHC_MISSENSE_SPLICING/PARENTAL_LIFESPAN_LIFEGEN_2017/5_MR_RESULTS/LD_0.1/MR_RESULTS.LD_PRUNED_0.1_EUROPEAN_CMET_F11_F11_CIS_200000_PVALUE_0.01_6.txt"
#########################################################


message(paste0("PANEL: ", PANEL))
message(paste0("EXPOSURES_INPUT: ", EXPOSURES_INPUT))
message(paste0("OUTCOMES_INPUT: ", OUTCOMES_INPUT))


# source("/genetics/PAREG/perrotn/automated_scripts/custom_functions.R")

EXPOSURES <- TwoSampleMR::read_exposure_data(EXPOSURES_INPUT, sep = "\t")
# EXPOSURE$units.exposure_dat<-"mmol/L"
# Read in Outcomes
OUTCOME <- TwoSampleMR::read_outcome_data(OUTCOMES_INPUT, sep = "\t")
# OUTCOME$units.outcome_dat<-"1SD"
# Harmonize Data

EXPOSURES = EXPOSURES[which(EXPOSURES$SNP %in% OUTCOME$SNP),]
OUTCOME = OUTCOME[which(OUTCOME$SNP %in% EXPOSURES$SNP),]

MR_MODEL=paste0(c(paste0(POPULATION, "_", PANEL, "_", ASSAY, "_", GENE), paste0(unique(EXPOSURES$exposure[EXPOSURES$exposure!=paste0(POPULATION, "_", PANEL, "_", ASSAY, "_", GENE)]), collapse = " ")), collapse = " ")

EXP_TO_OUT <- EXPOSURES[EXPOSURES$exposure != paste0(POPULATION, "_", PANEL, "_", ASSAY, "_", GENE), ]
colnames(EXP_TO_OUT) <- gsub(colnames(EXP_TO_OUT), pattern = "exposure", replacement = "outcome")

OUTCOMES <- rbind(OUTCOME, EXP_TO_OUT[, colnames(OUTCOME)])

if (file.exists(OUTPUT_FILE)) {
    file.remove(OUTPUT_FILE)
}

if (file.exists(gsub(OUTPUT_FILE, pattern = "RESULTS\\.LD", replacement = "IVs\\.LD"))) {
    file.remove(gsub(OUTPUT_FILE, pattern = "RESULTS\\.LD", replacement = "IVs\\.LD"))
}

EXP = paste0(POPULATION, "_", PANEL, "_", ASSAY, "_", GENE)

    message(
        "EXPOSURE IN PROGRESS: ",
        EXP,
        " (",
        which(unique(EXPOSURES$exposure) == EXP),
        "/",
        length(unique(EXPOSURES$exposure)),
        ")"
    )
    message(paste0("PANEL: ", PANEL))
    message(paste0("EXPOSURES_INPUT: ", EXPOSURES_INPUT))
    message(paste0("OUTCOMES_INPUT: ", OUTCOMES_INPUT))

COUNT=0
    for (OUT in unique(OUTCOMES$outcome)[which(unique(OUTCOMES$outcome) != EXP)])
    {
        DATA <- TwoSampleMR::harmonise_data(
            exposure_dat = EXPOSURES[which(EXPOSURES$exposure == EXP), ],
            outcome_dat = OUTCOMES[which(OUTCOMES$outcome == OUT), ],
            action = 1
        ) # Keep ambiguous snps because we know SNPs are aligned to + strand
        # DATA<-DATA[which(DATA$beta.outcome!=0),]

        # MR
        message("MR")

        if (length(which(DATA$mr_keep == T)) == 0)
        {
        next()
        }

        all_MR_res <- character()


        if (exists("MR_res")) {
            rm(MR_res)
        }

        if (length(which(DATA$mr_keep == T)) < 3) {
            message("LESS THAN 3 SNPS")
            if (length(which(DATA$mr_keep == T)) == 1) {
                message("ONE SNP")

                MR_res <- TwoSampleMR::mr(
                    DATA,
                    method_list = "mr_wald_ratio",
                    parameters = list(
                        over.dispersion = T,
                        loss.function = "tukey",
                        nboot = 1000
                    )
                )

                MR_res <- cbind(
                    MR_res,
                    over.dispersion = NA,
                    Q = NA,
                    Q_df = NA,
                    Q_pval = NA,
                    egger_intercept = NA,
                    egger_intercept_se = NA,
                    egger_intercept_pval = NA,
                    N_outliers = NA,
                    Distortion_Pval = NA,
                    stringsAsFactors = F
                )
            } else {
                message("TWO SNPS")
                MR_res <- TwoSampleMR::mr(
                    DATA,
                    method_list = c("mr_ivw", "mr_weighted_median"),
                    parameters = list(
                        over.dispersion = T,
                        loss.function = "tukey",
                        nboot = 1000
                    )
                )

                MR_res$over.dispersion <- NA
                MR_hetero <- TwoSampleMR::mr_heterogeneity(DATA,
                    method_list = "mr_ivw"
                )

                MR_res <- merge(MR_res, MR_hetero, all.x = T)

                MR_res <- cbind(
                    MR_res,
                    egger_intercept = NA,
                    egger_intercept_se = NA,
                    egger_intercept_pval = NA,
                    N_outliers = NA,
                    Distortion_Pval = NA,
                    stringsAsFactors = F
                )
            }
        } else {
            message("MORE THAN 2 SNPS")
            # >2 SNPs
            # TwoSampleMR::mr_method_list()

            MR_res <- TwoSampleMR::mr(
                DATA,
                method_list = c("mr_ivw", "mr_weighted_median"),
                parameters = list(
                    over.dispersion = T,
                    loss.function = "tukey",
                    nboot = 1000
                )
            )

            MR_EGGER <- tryCatch(
                {
                    TwoSampleMR::mr(
                        DATA,
                        method_list = "mr_egger_regression",
                        parameters = list(
                            over.dispersion = T,
                            loss.function = "tukey",
                            nboot = 1000
                        )
                    )
                },
                warning = function(w) {
                    if (w$message == "Collinearities in MR Egger, try LD pruning the exposure variables.") {
                        return(w$message)
                    }
                }
            )

            MR_RAPS <- tryCatch(
                {
                    cbind(TwoSampleMR::mr(
                        DATA,
                        method_list = c("mr_raps"),
                        parameters = list(
                            over.dispersion = T,
                            loss.function = "tukey",
                            nboot = 100000,
                            shrinkage = T
                        )
                    ),
                    over.dispersion = T
                    )
                },
                warning = function(w) {
                    if (w$message == "Did not converge when solving the estimating equations. Consider to increase niter or decrease tol.") {
                        MR_RAPS <- MR_res[1, , drop = F]
                        MR_RAPS[, "method"] <- "Robust adjusted profile score (RAPS)"
                        MR_RAPS[, "nsnp"] <- NA
                        MR_RAPS[, "b"] <- "Did not converge when solving the estimating equations (nboot=100000)"
                        MR_RAPS[, "se"] <- "Did not converge when solving the estimating equations (nboot=100000)"
                        MR_RAPS[, "pval"] <- "Did not converge when solving the estimating equations (nboot=100000)"
                        MR_RAPS$over.dispersion <- NA
                        return(MR_RAPS)
                    } else {
                        if (w$message == "The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion." |
                            w$call == "sqrt(tau2 + se_out^2 + se_exp^2 * beta^2)") {
                            cbind(
                                TwoSampleMR::mr(
                                    DATA,
                                    method_list = c("mr_raps"),
                                    parameters = list(
                                        over.dispersion = F,
                                        loss.function = "tukey",
                                        nboot = 1000,
                                        shrinkage = T
                                    )
                                ),
                                over.dispersion = F
                            )
                        } else {
                            if (w$message == "Estimated overdispersion seems abnormaly large") {
                                MR_RAPS <- MR_res[1, , drop = F]
                                MR_RAPS[, "method"] <- "Robust adjusted profile score (RAPS)"
                                MR_RAPS[, "nsnp"] <- NA
                                MR_RAPS[, "b"] <- "Estimated overdispersion seems abnormaly large"
                                MR_RAPS[, "se"] <- "Estimated overdispersion seems abnormaly large"
                                MR_RAPS[, "pval"] <- "Estimated overdispersion seems abnormaly large"
                                MR_RAPS$over.dispersion <- NA
                                return(MR_RAPS)
                            }
                        }
                    }
                },
                error = function(e) {
                    MR_RAPS <- MR_res[1, , drop = F]
                    MR_RAPS[, "method"] <- "Robust adjusted profile score (RAPS)"
                    MR_RAPS[, "nsnp"] <- NA
                    MR_RAPS[, "b"] <- e$message
                    MR_RAPS[, "se"] <- NA
                    MR_RAPS[, "pval"] <- NA
                    MR_RAPS$over.dispersion <- NA
                    return(MR_RAPS)
                }
            )



            if (all(class(MR_EGGER) == "character")) {
                MR_res <- plyr::rbind.fill(MR_res, cbind(
                    MR_res[1, c("id.exposure", "id.outcome", "outcome", "exposure")],
                    data.frame(
                        method = "MR Egger",
                        b = MR_EGGER,
                        stringsAsFactors = F
                    )
                ))
            } else {
                MR_res <- rbind(MR_res, MR_EGGER)
            }

            MR_res$over.dispersion <- NA

            MR_res <- rbind(MR_res, MR_RAPS)

            MR_hetero <- TwoSampleMR::mr_heterogeneity(DATA,
                method_list = c(
                    "mr_egger_regression",
                    "mr_ivw"
                )
            )


            MR_res <- merge(MR_res, MR_hetero, all.x = T)


            MR_res$egger_intercept <- NA
            MR_res$egger_intercept_se <- NA
            MR_res$egger_intercept_pval <- NA
            MR_pleio <- TwoSampleMR::mr_pleiotropy_test(DATA)
            MR_res$egger_intercept[which(as.character(MR_res$method) == "MR Egger")] <- MR_pleio$egger_intercept
            MR_res$egger_intercept_se[which(as.character(MR_res$method) == "MR Egger")] <- MR_pleio$se
            MR_res$egger_intercept_pval[which(as.character(MR_res$method) == "MR Egger")] <- MR_pleio$pval


            MR_res$N_outliers <- NA
            MR_res$Distortion_Pval <- NA


            MR_PRESSO <- function(x) {
                MRPRESSO::mr_presso(
                    BetaOutcome = "beta.outcome",
                    BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome",
                    SdExposure = "se.exposure",
                    data = DATA,
                    OUTLIERtest = T,
                    DISTORTIONtest = T,
                    SignifThreshold = 0.05,
                    NbDistribution = x,
                    seed = NULL
                )
            }


            if (any(MR_res$Q_pval < 0.05, na.rm = T)) {
                message("Q P-VALUE LOWER THAN 0.05")

                for (DISTRIBUTION in c(100, 500, 1000))
                {
                    MR_PRESSO_res <- tryCatch(
                        {
                            MR_PRESSO(DISTRIBUTION)
                        },
                        error = function(e) {

                        }
                    )

                    if (!is.null(MR_PRESSO_res)) {
                        break()
                    }
                } # END DISTRIBUTION LOOP

                if (!is.null(MR_PRESSO_res)) {
                    message("MR_PRESSO_res != NULL")

                    MR_PRESSO_res <- cbind(
                        MR_res[1:2, c("id.exposure", "id.outcome", "outcome", "exposure")],
                        data.frame(
                            method = paste0("MR-PRESSO ", as.character(
                                levels(MR_PRESSO_res$`Main MR results`$"MR Analysis")
                            )),
                            nsnp = rep(NA, 2),
                            b = MR_PRESSO_res$`Main MR results`$"Causal Estimate",
                            se = MR_PRESSO_res$`Main MR results`$"Sd",
                            pval = MR_PRESSO_res$`Main MR results`$"P-value",
                            over.dispersion = rep(NA, 2),
                            Q = rep(NA, 2),
                            Q_df = rep(NA, 2),
                            Q_pval = rep(NA, 2),
                            egger_intercept = rep(NA, 2),
                            egger_intercept_se = rep(NA, 2),
                            egger_intercept_pval = rep(NA, 2),
                            N_outliers = rep(
                                length(
                                    MR_PRESSO_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`[-which(
                                        MR_PRESSO_res$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices` ==
                                            "No significant outliers"
                                    )]
                                ),
                                2
                            ),
                            Distortion_Pval = ifelse(
                                is.null(
                                    MR_PRESSO_res$`MR-PRESSO results`$`Distortion Test`$`Pvalue`
                                ),
                                rep(NA, 2),
                                rep(
                                    MR_PRESSO_res$`MR-PRESSO results`$`Distortion Test`$`Pvalue`,
                                    2
                                )
                            ),
                            stringsAsFactors = F
                        )
                    )
                    MR_res <- rbind(MR_res, MR_PRESSO_res)


                    if (unique(na.omit(MR_res$N_outliers)) > 0) {
                        message("OUTLIERS SNPs ")

                        message(unique(MR_res$exposure))
                        message(unique(MR_res$outcome))
                        data.table::fwrite(
                            list(c(
                                as.character(unique(MR_res$exposure)), as.character(unique(MR_res$outcome))
                            )),
                            paste0(OUTPUT_FILE, "_outliers"),
                            quote = F,
                            row.names = F,
                            col.names = T,
                            sep = "\t",
                            append = F
                        )
                    }
                } else {
                    message("MR_PRESSO_res == NULL")

                    MR_res <- rbind(MR_res, cbind(
                        MR_res[1, c("id.exposure", "id.outcome", "outcome", "exposure")],
                        data.frame(
                            method = "MR PRESSO",
                            nsnp = NA,
                            b = "Not enough elements to compute empirical P-values",
                            se = NA,
                            pval = NA,
                            over.dispersion = NA,
                            Q = NA,
                            Q_df = NA,
                            Q_pval = NA,
                            egger_intercept = NA,
                            egger_intercept_se = NA,
                            egger_intercept_pval = NA,
                            N_outliers = NA,
                            Distortion_Pval = NA,
                            stringsAsFactors = F
                        )
                    ))
                }
            } else {
                message("Q P-VALUE HIGHER THAN 0.05")



                MR_res <- plyr::rbind.fill(MR_res, cbind(
                    MR_res[1, c("id.exposure", "id.outcome", "outcome", "exposure")],
                    data.frame(
                        method = "MR PRESSO",
                        b = "NOT PERFORMED Q_pval > 0.05",
                        stringsAsFactors = F
                    )
                ))

                # MR_res = rbind(MR_res,  cbind(
                #   MR_res[1, c("id.exposure", "id.outcome", "outcome", "exposure")],
                #   data.frame(
                #     method = "MR PRESSO",
                #     nsnp = NA,
                #     b = "NOT PERFORMED Q_pval > 0.05",
                #     se = NA,
                #     pval = NA,
                #     over.dispersion = NA,
                #     Q = NA,
                #     Q_df = NA,
                #     Q_pval = NA,
                #     egger_intercept = NA,
                #     egger_intercept_se = NA,
                #     egger_intercept_pval = NA,
                #     N_outliers = NA,
                #     Distortion_Pval = NA,
                #     stringsAsFactors = F
                #   )
                # ))
            }
        }

        if (exists("MR_res")) {
            MR_res$mean.samplesize.outcome <- mean(DATA$samplesize.outcome)
            MR_res$min.samplesize.outcome <- min(DATA$samplesize.outcome)
            MR_res$max.samplesize.outcome <- max(DATA$samplesize.outcome)
            MR_res$ncase.outcome <- unique(DATA$ncase.outcome)
            MR_res$ncontrol.outcome <- unique(DATA$ncontrol.outcome)
            MR_res$samplesize.exposure <- unique(DATA$samplesize.exposure)
            all_MR_res <- rbind(all_MR_res, MR_res)

            if (is.factor(all_MR_res$exposure)) {
                all_MR_res$exposure <- as.character(levels(all_MR_res$exposure))[all_MR_res$exposure]
            }
        }


        #### Directionality test ####
        if (is.null(DATA$ncase.exposure)| all(is.na(DATA$ncase.exposure))) {
            DATA$r.exposure <- mapply(
                DATA$samplesize.exposure,
                DATA$pval.exposure,
                FUN = function(sample_size,
                               pval) {
                    TwoSampleMR::get_r_from_pn(p = pval, n = sample_size)
                }
            )
        } else {
            DATA$r.exposure <- mapply(
                DATA$beta.exposure,
                DATA$eaf.exposure,
                DATA$ncase.exposure,
                DATA$ncontrol.exposure,
                FUN = function(beta,
                               eaf,
                               ncase,
                               ncontrol) {
                    TwoSampleMR::get_r_from_lor(
                        lor = beta,
                        af = eaf,
                        ncase = ncase,
                        ncontrol = ncontrol,
                        prevalence = ncase / ncontrol
                    )
                }
            )
        }

        # if ncase/ncontrol/samplesize == NA skip directionality test & F-stat
        if ((is.na(unique(DATA$ncase.outcome)) &
            is.na(unique(DATA$ncontrol.outcome)) &
            is.na(unique(DATA$samplesize.outcome))) | length(which(is.na(DATA$eaf.outcome))) > 0 | length(which(DATA$eaf.outcome == "")) > 0)
        # no n case/ctrl/total
            {
                all_MR_res$snp_r2.exposure <- NA
                all_MR_res$snp_r2.outcome <- NA
                all_MR_res$correct_causal_direction <- NA
                all_MR_res$steiger_pval <- NA
                all_MR_res$F_STAT <- NA

                message("NO N AVAILABLE, SKIP DIRECTIONALITY TEST & F-STAT")
            } else {
            # N case/ctrl or total
            if (is.na(unique(DATA$ncase.outcome)) &
                is.na(unique(DATA$ncontrol.outcome)))
            # no N case/ctrl only total == quantitative trait
                {
                    DATA$r.outcome <- mapply(
                        DATA$samplesize.outcome,
                        DATA$pval.outcome,
                        FUN = function(sample_size,
                                       pval) {
                            TwoSampleMR::get_r_from_pn(p = pval, n = sample_size)
                        }
                    )
                } else {
                # N case/ctrl == qualitative trait
                DATA$r.outcome <- mapply(
                    DATA$beta.outcome,
                    DATA$eaf.outcome,
                    DATA$ncase.outcome,
                    DATA$ncontrol.outcome,
                    DATA$samplesize.outcome,
                    DATA$pval.outcome,
                    FUN = function(beta,
                                   eaf,
                                   ncase,
                                   ncontrol,
                                   sample_size,
                                   pval) {
                        ifelse(
                            is.na(ncase),
                            TwoSampleMR::get_r_from_pn(p = pval, n = sample_size),
                            TwoSampleMR::get_r_from_lor(
                                lor = beta,
                                af = eaf,
                                ncase = ncase,
                                ncontrol = ncontrol,
                                prevalence = ncase / ncontrol
                            )
                        )
                    }
                )
            }

            res_steiger <- TwoSampleMR::directionality_test(DATA)
            all_MR_res <- merge(x = all_MR_res, y = res_steiger, all.x = T)

            #### F-STAT ####
            all_MR_res$F_STAT <-
                apply(all_MR_res[, c("snp_r2.exposure", "samplesize.exposure", "nsnp")], 1, function(x) {
                    (x[1] * (x[2] - 1 - x[3])) / ((1 - x[1]) * x[3])
                })
        }

        if (EXP == paste0(POPULATION, "_", PANEL, "_", ASSAY, "_", GENE))
        {
            ## POPULATION ##
            all_MR_res$POPULATION <- POPULATION


            ## PANEL ##
            all_MR_res$PANEL <- PANEL

            ## BIOMARKER ##
            all_MR_res$BIOMARKER <- ASSAY

            ## GENE ##
            all_MR_res$GENE <- GENE

        } else {
              ## POPULATION ##
            all_MR_res$POPULATION <- NA

            ## PANEL ##
            all_MR_res$PANEL <- NA

            ## BIOMARKER ##
            all_MR_res$BIOMARKER <- NA

            ## GENE ##
            all_MR_res$GENE <- NA

        }
    
        all_MR_res$MODEL = MR_MODEL

        #### Save results ####
        COUNT <- COUNT + 1
        data.table::fwrite(
            all_MR_res,
            OUTPUT_FILE,
            quote = F,
            row.names = F,
            col.names = COUNT == 1,
            sep = "\t",
            append = COUNT != 1
        )

        #### Save IVs ####
        if (length(which(DATA$mr_keep == F)) > 0) {
            DATA <- DATA[which(DATA$mr_keep == T), ]
        }

        data.table::fwrite(
            DATA,
            gsub(OUTPUT_FILE, pattern = "RESULTS\\.LD", replacement = "IVs\\.LD"),
            quote = F,
            row.names = F,
            col.names = COUNT == 1,
            sep = "\t",
            append = COUNT != 1
        )
    } # END OUT LOOP


TSMR_MV_DATA=mv_harmonise_data(exposure_dat = EXPOSURES, outcome_dat= OUTCOME, harmonise_strictness=1)

MV_DATA = MendelianRandomization::mr_mvinput(
bx = TSMR_MV_DATA$exposure_beta,
bxse = TSMR_MV_DATA$exposure_se,
by = TSMR_MV_DATA$outcome_beta,
byse = TSMR_MV_DATA$outcome_se,
exposure = TSMR_MV_DATA$expname["exposure"],
outcome = TSMR_MV_DATA$outname["outcome"],
snps = rownames(TSMR_MV_DATA$exposure_beta),
effect_allele = sapply(rownames(TSMR_MV_DATA$exposure_beta), FUN=function(x){toupper(unlist(strsplit(x, "_"))[4])}),
other_allele = sapply(rownames(TSMR_MV_DATA$exposure_beta), FUN=function(x){tolower(unlist(strsplit(x, "_"))[3])})
)

MR_MVIVW = mr_mvivw(
MV_DATA,
model = "default",
robust = FALSE,
correl = FALSE,
)

MR_MVIVW_ROB = mr_mvivw(
MV_DATA,
model = "default",
robust = TRUE,
correl = FALSE,
)

MR_MVLASSO=mr_mvlasso(
MV_DATA,
orientate = 1,
distribution = "normal")

MR_MVMEDIAN=mr_mvmedian(
MV_DATA,
distribution = "normal",
iterations = 10000)

MR_MVEGGER=mr_mvegger(
MV_DATA,
orientate = 1,
correl = FALSE,
distribution = "normal",
)


MVMR_RES <- data.frame(
    id.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    id.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    outcome = rep(MR_MVIVW$Outcome, nrow(TSMR_MV_DATA$expname["exposure"])),
    exposure = TSMR_MV_DATA$expname["exposure"],
    method = rep("Multivariable inverse-variance weighted", nrow(TSMR_MV_DATA$expname["exposure"])),
    nsnp = rep(MR_MVIVW$SNPs, nrow(TSMR_MV_DATA$expname["exposure"])),
    b = MR_MVIVW$Estimate,
    se = MR_MVIVW$StdError,
    pval = MR_MVIVW$Pvalue,
    over.dispersion = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q = rep(MR_MVIVW$Heter.Stat[1], nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_df = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_pval = rep(MR_MVIVW$Heter.Stat[2], nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_se = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    N_outliers = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Distortion_Pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    mean.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    min.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    max.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncase.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncontrol.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    samplesize.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    correct_causal_direction = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    steiger_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    F_STAT = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    POPULATION = rep(POPULATION, nrow(TSMR_MV_DATA$expname["exposure"])),
    PANEL = rep(PANEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    BIOMARKER = rep(ASSAY, nrow(TSMR_MV_DATA$expname["exposure"])),
    GENE = rep(GENE, nrow(TSMR_MV_DATA$expname["exposure"])),
    MODEL = rep(MR_MODEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    stringsAsFactors = F
)


MVMR_RES <- rbind(MVMR_RES, data.frame(
    id.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    id.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    outcome = rep(MR_MVIVW_ROB$Outcome, nrow(TSMR_MV_DATA$expname["exposure"])),
    exposure = TSMR_MV_DATA$expname["exposure"],
    method = rep("Multivariable inverse-variance weighted (Robust)", nrow(TSMR_MV_DATA$expname["exposure"])),
    nsnp = rep(MR_MVIVW_ROB$SNPs, nrow(TSMR_MV_DATA$expname["exposure"])),
    b = MR_MVIVW_ROB$Estimate,
    se = MR_MVIVW_ROB$StdError,
    pval = MR_MVIVW_ROB$Pvalue,
    over.dispersion = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_df = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_pval = rep(MR_MVIVW_ROB$Heter.Stat[2], nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_se = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    N_outliers = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Distortion_Pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    mean.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    min.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    max.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncase.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncontrol.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    samplesize.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    correct_causal_direction = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    steiger_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    F_STAT = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    POPULATION = rep(POPULATION, nrow(TSMR_MV_DATA$expname["exposure"])),
    PANEL = rep(PANEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    BIOMARKER = rep(ASSAY, nrow(TSMR_MV_DATA$expname["exposure"])),
    GENE = rep(GENE, nrow(TSMR_MV_DATA$expname["exposure"])),
    MODEL = rep(MR_MODEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    stringsAsFactors = F
))

MVMR_RES <- rbind(MVMR_RES, data.frame(
    id.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    id.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    outcome = rep(MR_MVLASSO$Outcome, nrow(TSMR_MV_DATA$expname["exposure"])),
    exposure = TSMR_MV_DATA$expname["exposure"],
    method = rep("Multivariable MR-Lasso", nrow(TSMR_MV_DATA$expname["exposure"])),
    nsnp = rep(MR_MVLASSO$Valid, nrow(TSMR_MV_DATA$expname["exposure"])),
    b = MR_MVLASSO$Estimate,
    se = MR_MVLASSO$StdError,
    pval = MR_MVLASSO$Pvalue,
    over.dispersion = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_df = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_se = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    N_outliers = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Distortion_Pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    mean.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    min.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    max.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncase.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncontrol.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    samplesize.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    correct_causal_direction = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    steiger_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    F_STAT = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    POPULATION = rep(POPULATION, nrow(TSMR_MV_DATA$expname["exposure"])),
    PANEL = rep(PANEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    BIOMARKER = rep(ASSAY, nrow(TSMR_MV_DATA$expname["exposure"])),
    GENE = rep(GENE, nrow(TSMR_MV_DATA$expname["exposure"])),
    MODEL = rep(MR_MODEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    stringsAsFactors = F
))

MVMR_RES <- rbind(MVMR_RES, data.frame(
    id.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    id.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    outcome = rep(MR_MVMEDIAN$Outcome, nrow(TSMR_MV_DATA$expname["exposure"])),
    exposure = TSMR_MV_DATA$expname["exposure"],
    method = rep("Multivariable median-based", nrow(TSMR_MV_DATA$expname["exposure"])),
    nsnp = rep(MR_MVMEDIAN$SNPs, nrow(TSMR_MV_DATA$expname["exposure"])),
    b = MR_MVMEDIAN$Estimate,
    se = MR_MVMEDIAN$StdError,
    pval = MR_MVMEDIAN$Pvalue,
    over.dispersion = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_df = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_se = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    N_outliers = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Distortion_Pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    mean.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    min.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    max.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncase.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncontrol.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    samplesize.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    correct_causal_direction = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    steiger_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    F_STAT = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    POPULATION = rep(POPULATION, nrow(TSMR_MV_DATA$expname["exposure"])),
    PANEL = rep(PANEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    BIOMARKER = rep(ASSAY, nrow(TSMR_MV_DATA$expname["exposure"])),
    GENE = rep(GENE, nrow(TSMR_MV_DATA$expname["exposure"])),
    MODEL = rep(MR_MODEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    stringsAsFactors = F
))


MVMR_RES <- rbind(MVMR_RES, data.frame(
    id.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    id.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    outcome = rep(MR_MVEGGER$Outcome, nrow(TSMR_MV_DATA$expname["exposure"])),
    exposure = TSMR_MV_DATA$expname["exposure"],
    method = rep("Multivariable MR-Egger", nrow(TSMR_MV_DATA$expname["exposure"])),
    nsnp = rep(MR_MVEGGER$SNPs, nrow(TSMR_MV_DATA$expname["exposure"])),
    b = MR_MVEGGER$Estimate,
    se = MR_MVEGGER$StdError.Est,
    pval = MR_MVEGGER$Pvalue.Est,
    over.dispersion = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q = rep(MR_MVEGGER$Heter.Stat[1], nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_df = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Q_pval = rep(MR_MVEGGER$Heter.Stat[2], nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept = MR_MVEGGER$Estimate,
    egger_intercept_se = rep(MR_MVEGGER$StdError.Int, nrow(TSMR_MV_DATA$expname["exposure"])),
    egger_intercept_pval = rep(MR_MVEGGER$Pvalue.Int, nrow(TSMR_MV_DATA$expname["exposure"])),
    N_outliers = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    Distortion_Pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    mean.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    min.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    max.samplesize.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncase.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    ncontrol.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    samplesize.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.exposure = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    snp_r2.outcome = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    correct_causal_direction = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    steiger_pval = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    F_STAT = rep(NA, nrow(TSMR_MV_DATA$expname["exposure"])),
    POPULATION = rep(POPULATION, nrow(TSMR_MV_DATA$expname["exposure"])),
    PANEL = rep(PANEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    BIOMARKER = rep(ASSAY, nrow(TSMR_MV_DATA$expname["exposure"])),
    GENE = rep(GENE, nrow(TSMR_MV_DATA$expname["exposure"])),
    MODEL = rep(MR_MODEL, nrow(TSMR_MV_DATA$expname["exposure"])),
    stringsAsFactors = F
))


data.table::fwrite(
            MVMR_RES,
            OUTPUT_FILE,
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t",
            append = T
        )



  message("MR DONE")

