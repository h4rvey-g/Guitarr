library(tidyverse)
library(psych)
library(Rgb)
library(tximport)
cells <- c("GS689", "PC3E")
combinations <- expand.grid(cells = cells, suffix = c("_1", "_2", "_3"), stringsAsFactors = FALSE)
#####################
# 01. stringtie without imputation
#####################
# create a filelist for all t_data.ctab files in data/04.Tx_quantification/01.Stringtie
stringtie_tx_files <- list.files("data/04.Tx_quantification/01.Stringtie", pattern = "t_data.ctab", full.names = TRUE, recursive = TRUE)

# use tximport to read all t_data.ctab files in stringtie_tx_files
tx2gene <- read_tsv(stringtie_tx_files[1]) %>% select(t_name, gene_name)
stringtie_tx <- list(
    GS689 = tximport(stringtie_tx_files[grepl("GS689", stringtie_tx_files)], type = "stringtie", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, txOut = TRUE)$abundance,
    PC3E = tximport(stringtie_tx_files[grepl("PC3E", stringtie_tx_files)], type = "stringtie", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, txOut = TRUE)$abundance
) %>%
    map(~ as.data.frame(.x) %>%
        rownames_to_column("Transcript_ID") %>%
        as_tibble())

#####################
# 02. stringtie with imputation
#####################
stringtie_impute_files <- list.files("data/04.Tx_quantification/02.Stringtie_with_impute", pattern = "t_data.ctab", full.names = TRUE, recursive = TRUE)
tx2gene <- read_tsv(stringtie_impute_files[1]) %>% select(t_name, gene_name)
stringtie_impute <- list(
    GS689 = tximport(stringtie_impute_files[grepl("GS689", stringtie_impute_files)], type = "stringtie", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, txOut = TRUE)$abundance,
    PC3E = tximport(stringtie_impute_files[grepl("PC3E", stringtie_impute_files)], type = "stringtie", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, txOut = TRUE)$abundance
) %>%
    map(~ as.data.frame(.x) %>%
        rownames_to_column("Transcript_ID") %>%
        as_tibble())
#####################
# 03. Guitar
#####################
# read tsv from all cells like data/04.Tx_quantification/03.Guitar/GS689_1Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt
guitar <- pmap(combinations, function(cells, suffix) {
    read_tsv(paste("data/04.Tx_quantification/03.Guitar/", cells, suffix, "Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt", sep = ""))
})
names(guitar) <- paste(combinations$cells, combinations$suffix, sep = "")

# # keep Transcript_ID, and Transcript_Exp columns
# guitar <- map(guitar, ~ select(.x, Transcript_ID, Transcript_Exp, Transcript_length))

# # calculate TPM based on Transcript_Exp and Transcript_length
# # calculate RPK
# guitar <- map(guitar, ~ mutate(.x, RPK = Transcript_Exp / (Transcript_length / 1000)))

# # calculate total RPK
# total_RPK <- sum(sapply(guitar, function(x) sum(x$RPK)))

# # calculate TPM
# guitar <- map(guitar, ~ mutate(.x, TPM = (RPK / total_RPK) * 10^6))
# guitar <- map(guitar, ~ select(.x, Transcript_ID, TPM))

# bind columns of bam_impute respectively for GS689 and PC3E
guitar_combined <- list()
guitar_combined$GS689 <- inner_join(guitar$GS689_1, guitar$GS689_2, by = "Transcript_ID", suffix = c(".1", ".2")) %>%
    inner_join(., guitar$GS689_3, by = "Transcript_ID", suffix = c("", ".3")) %>%
    select(Transcript_ID, Transcript_Exp.1, Transcript_Exp.2)
guitar_combined$PC3E <- inner_join(guitar$PC3E_1, guitar$PC3E_2, by = "Transcript_ID", suffix = c(".1", ".2")) %>%
    inner_join(., guitar$PC3E_3, by = "Transcript_ID", suffix = c("", ".3")) %>%
    select(Transcript_ID, Transcript_Exp.1, Transcript_Exp.2)

#####################
# 04. three_gen
#####################
# read tsv test/3rd_gen/3rd_gen.txt, and drop columns begin with HEK293T
three_gen <- read_tsv("data/reference/3rd_gen/3rd_gen.txt") %>%
    select(-starts_with("HEK293T")) %>%
    as_tibble()
three_gen_gtf <- read.gtf("data/reference/3rd_gen/3rd_gen.gtf") %>%
    as_tibble() %>%
    filter(feature == "transcript")

# calculate length of each transcript
three_gen_gtf <- three_gen_gtf %>%
    mutate(transcript_length = three_gen_gtf$end - three_gen_gtf$start) %>%
    select(transcript_id, transcript_length)

# # remove _[digit] from Transcript_ID
# three_gen$transcript_ID <- three_gen$transcript_ID %>% str_replace("_\\d+$", "")
three_gen_list <- list()
three_gen_list$GS689 <- three_gen %>%
    select(transcript_ID, starts_with("GS689")) %>%
    dplyr::rename(long_read_1 = GS689_1, long_read_2 = GS689_2, long_read_3 = GS689_3) %>%
    inner_join(., three_gen_gtf, by = c("transcript_ID" = "transcript_id")) %>%
    mutate(transcript_ID = str_remove(transcript_ID, "_\\d+$"))
three_gen_list$PC3E <- three_gen %>%
    select(transcript_ID, starts_with("PC3E")) %>%
    dplyr::rename(long_read_1 = PC3E_1, long_read_2 = PC3E_2, long_read_3 = PC3E_3) %>%
    inner_join(., three_gen_gtf, by = c("transcript_ID" = "transcript_id")) %>%
    mutate(transcript_ID = str_remove(transcript_ID, "_\\d+$"))

# 计算 RPK
three_gen_list$GS689 <- three_gen_list$GS689 %>%
    mutate(across(starts_with("long_read"), ~ . / transcript_length * 1000, .names = "RPK_{.col}"))

three_gen_list$PC3E <- three_gen_list$PC3E %>%
    mutate(across(starts_with("long_read"), ~ . / transcript_length * 1000, .names = "RPK_{.col}"))

# 计算总的 RPK
total_RPK_GS689 <- select(three_gen_list$GS689, starts_with("RPK_")) %>%
    sapply(sum) %>%
    sum()
total_RPK_PC3E <- select(three_gen_list$PC3E, starts_with("RPK_")) %>%
    sapply(sum) %>%
    sum()

# 计算 TPM
three_gen_list$GS689 <- three_gen_list$GS689 %>%
    mutate(across(starts_with("RPK_"), ~ (. / total_RPK_GS689) * 10^6, .names = "TPM_{str_remove(.col, 'RPK_')}")) %>%
    select(transcript_ID, starts_with("TPM"))

three_gen_list$PC3E <- three_gen_list$PC3E %>%
    mutate(across(starts_with("RPK_"), ~ (. / total_RPK_PC3E) * 10^6, .names = "TPM_{str_remove(.col, 'RPK_')}")) %>%
    select(transcript_ID, starts_with("TPM"))
#####################
# 05. Calculate correlations
#####################
calculate_correlations <- function(short_list, name) {
    # rename columns [2:4] to short_read_1, short_read_2, short_read_3
    colnames(short_list$GS689)[2:4] <- paste("short_read", 1:3, sep = "_")
    colnames(short_list$PC3E)[2:4] <- paste("short_read", 1:3, sep = "_")

    # Combine three_gen and bam_impute$GS689, based on transcript_ID and Transcript_ID, add prefix
    GS689 <- inner_join(three_gen_list$GS689, short_list$GS689, by = c("transcript_ID" = "Transcript_ID"))
    PC3E <- inner_join(three_gen_list$PC3E, short_list$PC3E, by = c("transcript_ID" = "Transcript_ID"))

    remove_outliers <- function(df, col_name) {
        Q1 <- quantile(df[[col_name]], 0.25)
        Q3 <- quantile(df[[col_name]], 0.75)
        IQR <- Q3 - Q1
        filter(df, between(df[[col_name]], Q1 - 1.5 * IQR, Q3 + 1.5 * IQR))
    }

    # get average of long_read_1, long_read_2, long_read_3, short_read_1, short_read_2, short_read_3 for GS689 and PC3E
    GS689 <- mutate(GS689, long_read_avg = rowMeans(select(GS689, starts_with("TPM_long_read"))), short_read_avg = rowMeans(select(GS689, starts_with("short_read")))) %>%
        remove_outliers(., "short_read_avg") %>%
        remove_outliers(., "long_read_avg")
    PC3E <- mutate(PC3E, long_read_avg = rowMeans(select(PC3E, starts_with("TPM_long_read"))), short_read_avg = rowMeans(select(PC3E, starts_with("short_read")))) %>%
        remove_outliers(., "short_read_avg") %>%
        remove_outliers(., "long_read_avg")

    # write GS689 and PC3E to tsv files
    write_tsv(GS689, paste("data/05.correlations/", name, "_GS689.tsv", sep = ""))
    write_tsv(PC3E, paste("data/05.correlations/", name, "_PC3E.tsv", sep = ""))

    # calculate correlation between long_read_avg and short_read_avg for GS689 and PC3E
    GS689_cor <- cor.test(GS689$long_read_avg, GS689$short_read_avg, method = "spearman")
    PC3E_cor <- cor.test(PC3E$long_read_avg, PC3E$short_read_avg, method = "spearman")

    # draw scatter plot for long_read_avg and short_read_avg for GS689 and PC3E
    GS689_plot <- ggplot(GS689, aes(x = long_read_avg, y = short_read_avg)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        ggtitle(paste("GS689, correlation = ", round(GS689_cor$estimate, 2), ", p-value = ", format(GS689_cor$p.value, scientific = TRUE)))
    ggsave(paste("data/05.correlations/", name, "_GS689.png", sep = ""), GS689_plot)
    PC3E_plot <- ggplot(PC3E, aes(x = long_read_avg, y = short_read_avg)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        ggtitle(paste("PC3E, correlation = ", round(PC3E_cor$estimate, 2), ", p-value = ", format(PC3E_cor$p.value, scientific = TRUE)))

    ggsave(paste("data/05.correlations/", name, "_PC3E.png", sep = ""), PC3E_plot)

    # return the correlation results
    # list(GS689_correlation = GS689_cor, PC3E_correlation = PC3E_cor)
}

calculate_correlations(stringtie_tx, "stringtie_vs_three_gen")
calculate_correlations(stringtie_impute, "stringtie_impute_vs_three_gen")
calculate_correlations(guitar_combined, "guitar_vs_three_gen")

calculate_correlations <- function(short_list1, short_list2, name) {
    # rename columns [2:4] to short_read_1, short_read_2, short_read_3
    colnames(short_list1$GS689)[2:4] <- paste("short_read", 1:3, sep = "_")
    colnames(short_list1$PC3E)[2:4] <- paste("short_read", 1:3, sep = "_")
    colnames(short_list2$GS689)[2:4] <- paste("short_read", 1:3, sep = "_")
    colnames(short_list2$PC3E)[2:4] <- paste("short_read", 1:3, sep = "_")

    # Combine three_gen and bam_impute$GS689, based on transcript_ID and Transcript_ID, add prefix
    GS689 <- inner_join(short_list1$GS689, short_list2$GS689, by = "Transcript_ID", suffix = c(".1", ".2"))
    PC3E <- inner_join(short_list1$PC3E, short_list2$PC3E, by = "Transcript_ID", suffix = c(".1", ".2"))

    # get average of long_read_1, long_read_2, long_read_3, short_read_1, short_read_2, short_read_3 for GS689 and PC3E
    GS689 <- mutate(GS689, short_read_avg_1 = rowMeans(select(GS689, starts_with("short_read_1"))), short_read_avg_2 = rowMeans(select(GS689, starts_with("short_read_2"))))
    PC3E <- mutate(PC3E, short_read_avg_1 = rowMeans(select(PC3E, starts_with("short_read_1"))), short_read_avg_2 = rowMeans(select(PC3E, starts_with("short_read_2"))))

    # write GS689 and PC3E to tsv files
    write_tsv(GS689, paste("data/05.correlations/", name, "_GS689.tsv", sep = ""))
    write_tsv(PC3E, paste("data/05.correlations/", name, "_PC3E.tsv", sep = ""))

    # calculate correlation between long_read_avg and short_read_avg for GS689 and PC3E
    GS689_cor <- cor.test(GS689$short_read_avg_1, GS689$short_read_avg_2)
    PC3E_cor <- cor.test(PC3E$short_read_avg_1, PC3E$short_read_avg_2)

    # draw scatter plot for long_read_avg and short_read_avg for GS689 and PC3E
    GS689_plot <- ggplot(GS689, aes(x = short_read_avg_1, y = short_read_avg_2)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        ggtitle(paste("GS689, correlation = ", round(GS689_cor$estimate, 2), ", p-value = ", format(GS689_cor$p.value, scientific = TRUE)))
    ggsave(paste("data/05.correlations/", name, "_GS689.png", sep = ""), GS689_plot)
    PC3E_plot <- ggplot(PC3E, aes(x = short_read_avg_1, y = short_read_avg_2)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        ggtitle(paste("PC3E, correlation = ", round(PC3E_cor$estimate, 2), ", p-value = ", format(PC3E_cor$p.value, scientific = TRUE)))
    ggsave(paste("data/05.correlations/", name, "_PC3E.png", sep = ""), PC3E_plot)

    # return the correlation results
    # list(GS689_correlation = GS689_cor, PC3E_correlation = PC3E_cor)
}
calculate_correlations(stringtie_tx, stringtie_impute, "stringtie_vs_stringtie_impute")
calculate_correlations(stringtie_tx, guitar_combined, "stringtie_vs_guitar")
calculate_correlations(stringtie_impute, guitar_combined, "stringtie_impute_vs_guitar")
