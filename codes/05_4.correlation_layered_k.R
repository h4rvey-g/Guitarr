library(psych)
library(Rgb)
library(tximport)
library(outliers)
library(biomaRt)
library(tidyverse)
k_values <- read_tsv("data/reference/GRCh38_chr_K_value/kvalues.out")
cells <- c("GS689", "PC3E")
combinations <- expand.grid(cells = cells, suffix = c("_1", "_2", "_3"), stringsAsFactors = FALSE)
#####################
# 01. stringtie without imputation
#####################
# create a filelist for all *.tab files in data/04.Tx_quantification/01.Stringtie
stringtie_tx_files <- list.files("data/04.Tx_quantification/01.Stringtie", pattern = ".*\\.tab", full.names = TRUE, recursive = TRUE)

# use tximport to read all t_data.ctab files in stringtie_tx_files
tx2gene <- read_tsv(stringtie_tx_files[1]) %>% select("Gene ID", "TPM")

stringtie_tx <- list()
# read all files in stringtie_tx_files, and store them in stringtie_tx, use read_tsv
stringtie_tx <- map(stringtie_tx_files, read_tsv) %>%
    map(~ select(.x, "Gene ID", "TPM"))

# bind columns of bam_impute respectively for GS689 and PC3E
stringtie_tx_combined$GS689 <- inner_join(stringtie_tx$GS689_1, stringtie_tx$GS689_2, by = "Gene ID", suffix = c(".1", ".2")) %>%
    inner_join(., stringtie_tx$GS689_3, by = "Gene ID", suffix = c("", ".3")) %>%
    dplyr::rename(
        V1 = TPM, V2 = TPM.1, V3 = TPM.2
    ) %>%
    mutate("Gene ID" = str_remove(`Gene ID`, "\\..*"))
stringtie_tx_combined$PC3E <- inner_join(stringtie_tx$PC3E_1, stringtie_tx$PC3E_2, by = "Gene ID", suffix = c(".1", ".2")) %>%
    inner_join(., stringtie_tx$PC3E_3, by = "Gene ID", suffix = c("", ".3")) %>%
    dplyr::rename(
        V1 = TPM, V2 = TPM.1, V3 = TPM.2
    ) %>%
    mutate("Gene ID" = str_remove(`Gene ID`, "\\..*"))

#####################
# 02. stringtie with imputation
#####################
stringtie_impute_files <- list.files("data/04.Tx_quantification/02.Stringtie_with_impute", pattern = ".*\\.tab", full.names = TRUE, recursive = TRUE)
tx2gene_impute <- read_tsv(stringtie_impute_files[1]) %>% select("Gene ID", "TPM")
stringtie_impute <- list()
stringtie_impute <- map(stringtie_impute_files, read_tsv) %>%
    map(~ select(.x, "Gene ID", "TPM"))

stringtie_impute_combined$GS689 <- inner_join(stringtie_impute$GS689_1, stringtie_impute$GS689_2, by = "Gene ID", suffix = c(".1", ".2")) %>%
    inner_join(., stringtie_impute$GS689_3, by = "Gene ID", suffix = c("", ".3")) %>%
    dplyr::rename(
        V1 = TPM, V2 = TPM.1, V3 = TPM.2
    ) %>%
    mutate("Gene ID" = str_remove(`Gene ID`, "\\..*"))
stringtie_impute_combined$PC3E <- inner_join(stringtie_impute$PC3E_1, stringtie_impute$PC3E_2, by = "Gene ID", suffix = c(".1", ".2")) %>%
    inner_join(., stringtie_impute$PC3E_3, by = "Gene ID", suffix = c("", ".3")) %>%
    dplyr::rename(
        V1 = TPM, V2 = TPM.1, V3 = TPM.2
    ) %>%
    mutate("Gene ID" = str_remove(`Gene ID`, "\\..*"))
#####################
# 03. Guitar
#####################
# read tsv from all cells like data/04.Tx_quantification/03.Guitar/GS689_1Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt
guitar <- pmap(combinations, function(cells, suffix) {
    read_tsv(paste("data/04.Tx_quantification/03.Guitar/", cells, suffix, "Aligned.sortedByCoord.out.readsImpute.gene.expression.txt", sep = ""))
})
names(guitar) <- paste(combinations$cells, combinations$suffix, sep = "")

# keep Transcript_ID, and Transcript_Exp columns
guitar <- map(guitar, ~ dplyr::select(.x, Gene_ID, Reference, Start, End, Transcript_Exp, Transcript_ID))
# combine the Gene_ID and Reference columns
guitar <- map(guitar, ~ mutate(.x, Gene_ID_Ref = paste(Gene_ID, Reference, Start, End, sep = "_")))

# bind columns of bam_impute respectively for GS689 and PC3E
guitar_combined <- list()
guitar_combined$GS689 <- inner_join(guitar$GS689_1, guitar$GS689_2, by = "Gene_ID_Ref", suffix = c(".1", ".2")) %>%
    inner_join(., guitar$GS689_3, by = "Gene_ID_Ref", suffix = c("", ".3")) %>%
    dplyr::rename(V1 = Transcript_Exp, V2 = Transcript_Exp.1, V3 = Transcript_Exp.2)
guitar_combined$PC3E <- inner_join(guitar$PC3E_1, guitar$PC3E_2, by = "Gene_ID_Ref", suffix = c(".1", ".2")) %>%
    inner_join(., guitar$PC3E_3, by = "Gene_ID_Ref", suffix = c("", ".3")) %>%
    dplyr::rename(V1 = Transcript_Exp, V2 = Transcript_Exp.1, V3 = Transcript_Exp.2)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
BM_GS689 <- guitar_combined$GS689 %>%
    mutate(Transcript_ID_one = str_split(Transcript_ID, ",") %>% map_chr(~ .x[1]) %>% str_remove("\\.\\d+$")) %>%
    getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"), filters = "ensembl_transcript_id", values = .$Transcript_ID_one, mart = mart) %>%
    as_tibble()
guitar_combined$GS689 <- inner_join(guitar_combined$GS689 %>% mutate(Transcript_ID_one = str_split(Transcript_ID, ",") %>% map_chr(~ .x[1]) %>% str_remove("\\.\\d+$")),
    BM_GS689,
    by = c("Transcript_ID_one" = "ensembl_transcript_id")
) %>%
    dplyr::rename("Gene ID" = ensembl_gene_id)

#####################
# 04. three_gen
#####################
# read tsv test/3rd_gen/3rd_gen.txt, and drop columns begin with HEK293T
three_gen <- read_tsv("data/reference/3rd_gen/3rd_gen.txt") %>%
    select(-starts_with("HEK293T"))
# group_by(gene_ID) %>%
# summarise(across(GS689_1:PC3E_3, sum, na.rm = TRUE))
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
    select(transcript_ID, gene_ID, starts_with("GS689")) %>%
    dplyr::rename(long_read_1 = GS689_1, long_read_2 = GS689_2, long_read_3 = GS689_3) %>%
    inner_join(., three_gen_gtf, by = c("transcript_ID" = "transcript_id")) %>%
    mutate(transcript_ID = str_remove(transcript_ID, "_\\d+$"))
three_gen_list$PC3E <- three_gen %>%
    select(transcript_ID, gene_ID, starts_with("PC3E")) %>%
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
    group_by(gene_ID) %>%
    summarise(across(starts_with("TPM"), sum, na.rm = TRUE)) %>%
    select(gene_ID, starts_with("TPM"))

three_gen_list$PC3E <- three_gen_list$PC3E %>%
    mutate(across(starts_with("RPK_"), ~ (. / total_RPK_PC3E) * 10^6, .names = "TPM_{str_remove(.col, 'RPK_')}")) %>%
    group_by(gene_ID) %>%
    summarise(across(starts_with("TPM"), sum, na.rm = TRUE)) %>%
    select(gene_ID, starts_with("TPM"))

#####################
# 05. calculate correlation
#####################
k_values <- k_values_ori %>%
    mutate(Gene = str_remove(Gene, "\\..*")) %>%
    filter(!is.na(Kvalue)) %>%
    filter(Gene %in% all_gene_ID)
# remove outliers of k_values$Kvalue, which are outside of the range of (Q1 - 1.5 * IQR, Q3 + 1.5 * IQR)
Q1 <- quantile(k_values$Kvalue, 0.25)
Q3 <- quantile(k_values$Kvalue, 0.75)
IQR <- Q3 - Q1
k_values <- filter(k_values, between(Kvalue, Q1 - 1.5 * IQR, Q3 + 1.5 * IQR))

k_values_quantile <- k_values %>%
    filter(Kvalue > 1) %>%
    pull(Kvalue) %>%
    quantile(., probs = c(0.33, 0.66)) %>%
    round(2)
# separate k_values$Kvalue to 4 groups, one is =1, rest are divided by k_values_quantile
k_values <- k_values %>%
    mutate(k = case_when(
        Kvalue == 1 ~ "=1",
        Kvalue <= k_values_quantile[1] ~ paste0("1-", k_values_quantile[1]),
        Kvalue <= k_values_quantile[2] ~ paste0(k_values_quantile[1], "-", k_values_quantile[2]),
        TRUE ~ paste0(k_values_quantile[2], "-")
    )) %>%
    mutate(
        k_num = case_when(
            Kvalue == 1 ~ 1,
            Kvalue <= k_values_quantile[1] ~ k_values_quantile[1],
            Kvalue <= k_values_quantile[2] ~ k_values_quantile[2],
            TRUE ~ Inf
        ),
        k = fct_reorder(k, k_num)
    )

calculate_correlations <- function(short_list, name, k_values) {
    k_values_label <- levels(k_values$k)

    # rename columns V1, V2, V3 to short_read_1, short_read_2, short_read_3
    short_list <- map(short_list, ~ rename(.x, short_read_1 = V1, short_read_2 = V2, short_read_3 = V3))

    # Combine three_gen and bam_impute$GS689, based on transcript_ID and Transcript_ID, add prefix
    GS689 <- inner_join(three_gen_list$GS689, short_list$GS689, by = c("gene_ID" = "Gene ID"))
    PC3E <- inner_join(three_gen_list$PC3E, short_list$PC3E, by = c("gene_ID" = "Gene ID"))

    # get average of long_read_1, long_read_2, long_read_3, short_read_1, short_read_2, short_read_3 for GS689 and PC3E
    GS689 <- mutate(GS689, long_read_avg = rowMeans(select(GS689, starts_with("TPM_long_read"))), short_read_avg = rowMeans(select(GS689, starts_with("short_read"))))
    PC3E <- mutate(PC3E, long_read_avg = rowMeans(select(PC3E, starts_with("TPM_long_read"))), short_read_avg = rowMeans(select(PC3E, starts_with("short_read"))))

    remove_outliers <- function(df, col_name) {
        Q1 <- quantile(df[[col_name]], 0.25)
        Q3 <- quantile(df[[col_name]], 0.75)
        IQR <- Q3 - Q1

        filter(df, between(df[[col_name]], Q1 - 1.5 * IQR, Q3 + 1.5 * IQR))
    }

    GS689 <- remove_outliers(GS689, "long_read_avg")
    GS689 <- remove_outliers(GS689, "short_read_avg")

    PC3E <- remove_outliers(PC3E, "long_read_avg")
    PC3E <- remove_outliers(PC3E, "short_read_avg")

    # write GS689 and PC3E to tsv files
    write_tsv(GS689, paste("data/05.correlations/sep_by_k/", name, "_GS689_k.tsv", sep = ""))
    write_tsv(PC3E, paste("data/05.correlations/sep_by_k/", name, "_PC3E_k.tsv", sep = ""))

    # inner join k_values to GS689 and PC3E, based on gene_ID and Gene, and group by k
    GS689 <- inner_join(GS689, k_values, by = c("gene_ID" = "Gene")) %>%
        group_by(k) %>%
        nest()
    PC3E <- inner_join(PC3E, k_values, by = c("gene_ID" = "Gene")) %>%
        group_by(k) %>%
        nest()

    get_result <- function(GS689, PC3E, limits_label, sep_method) {
        GS689_cor <- cor.test(GS689$long_read_avg, GS689$short_read_avg, method = "spearman")
        PC3E_cor <- cor.test(PC3E$long_read_avg, PC3E$short_read_avg, method = "spearman")
        GS689_plot <- ggplot(GS689, aes(x = long_read_avg, y = short_read_avg)) +
            geom_point() +
            geom_smooth(method = "lm", se = FALSE) +
            ggtitle(
                paste(
                    "GS689, correlation = ", round(GS689_cor$estimate, 2), ", p-value = ", format(GS689_cor$p.value, scientific = TRUE),
                    ", transcript number = ", nrow(GS689), ", range = ", limits_label
                )
            ) +
            ylab(paste("short_read_avg (", name, ")", sep = ""))
        ggsave(paste("data/05.correlations/", sep_method, "/", name, "_vs_3rd_gen_GS689_", limits_label, ".png", sep = ""), GS689_plot, width = 9)

        PC3E_plot <- ggplot(PC3E, aes(x = long_read_avg, y = short_read_avg)) +
            geom_point() +
            geom_smooth(method = "lm", se = FALSE) +
            ggtitle(
                paste(
                    "PC3E, correlation = ", round(PC3E_cor$estimate, 2), ", p-value = ", format(PC3E_cor$p.value, scientific = TRUE),
                    ", transcript number = ", nrow(PC3E), ", range = ", limits_label
                )
            ) +
            ylab(paste("short_read_avg (", name, ")", sep = ""))
        ggsave(paste("data/05.correlations/", sep_method, "/", name, "_vs_3rd_gen_PC3E_", limits_label, ".png", sep = ""), PC3E_plot, width = 9)

        # return the correlation coefficient
        c(GS689_cor$estimate, PC3E_cor$estimate)
    }
    single_cor <- pmap(list(GS689$data, PC3E$data, k_values_label, "sep_by_k"), get_result)
    # names(single_cor) <- c("<1000", "1000-5000", "5000-10000", ">10000")
}

cor_sum <- list()
cor_sum[["Stringtie"]] <- calculate_correlations(stringtie_tx_combined, "Stringtie", k_values)
cor_sum[["Stringtie_impute"]] <- calculate_correlations(stringtie_impute_combined, "Stringtie_impute", k_values)
cor_sum[["Guitar"]] <- calculate_correlations(guitar_combined, "Guitar", k_values)

# use list_cbind to combine all cor_sum
cor_sum_df <- map(cor_sum, ~ map_dfr(.x, ~ as.data.frame(t(.x)))) %>%
    # use map to rename colnames to c("GS689", "PC3E")
    map(~ setNames(.x, c("GS689", "PC3E"))) %>%
    list_cbind() %>%
    unnest(cols = c(Stringtie, Stringtie_impute, Guitar), names_sep = "_") %>%
    # add a col of c("<1000", "1000-5000", "5000-10000", ">10000") to the result
    cbind(Kvalue = k_values$k %>% levels() %>% factor(., levels = .)) %>%
    # replace $ with _ in colnames
    rename_all(~ str_replace_all(., "\\$", "_"))

# use geom_dotplot for cor_sum_df, with x = Stringtie, Stringtie_impute, Guitar, y = value, color = length, facet = cell
p <- cor_sum_df %>%
    gather(key = "method", value = "value", -Kvalue) %>%
    # split method to two cols, method and cell, sep by the last _, use str_split_fixed
    mutate(cell = str_extract(method, "(?<=_)[[:alnum:]]+$"), method = str_extract(method, ".*(?=_)")) %>%
    ggplot(aes(x = method, y = value, fill = Kvalue)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
    # add line to connect the dots of the same fill
    geom_line(aes(group = interaction(Kvalue, cell), color = Kvalue)) +
    # make xlab follow the order of Stringtie, Stringtie_impute, Guitar
    scale_x_discrete(limits = c("Stringtie", "Stringtie_impute", "Guitar")) +
    facet_wrap(~cell) +
    ylab("Spearman correlation coefficient")
ggsave("data/05.correlations/sep_by_k/cor_sum.png", p, width = 8)
