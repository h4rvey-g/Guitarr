library(tidyverse)
library(psych)
library(Rgb)
library(tximport)
library(outliers)
cells <- c("GS689", "PC3E")
combinations <- expand.grid(cells = cells, suffix = c("_1", "_2", "_3"), stringsAsFactors = FALSE)
#####################
# 01. stringtie without imputation
#####################
# create a filelist for all t_data.ctab files in data/04.Tx_quantification/01.Stringtie
stringtie_tx_files <- list.files("data/04.Tx_quantification/01.Stringtie", pattern = "t_data.ctab", full.names = TRUE, recursive = TRUE)

# use tximport to read all t_data.ctab files in stringtie_tx_files
tx2gene <- read_tsv(stringtie_tx_files[1]) %>% select(t_name, gene_name, length, num_exons)
stringtie_tx <- list(
    GS689 = tximport(stringtie_tx_files[grepl("GS689", stringtie_tx_files)], type = "stringtie", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, txOut = TRUE)$abundance,
    PC3E = tximport(stringtie_tx_files[grepl("PC3E", stringtie_tx_files)], type = "stringtie", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, txOut = TRUE)$abundance
) %>%
    map(~ as.data.frame(.x) %>%
        rownames_to_column("Transcript_ID") %>%
        as_tibble())
# innerjoin stringtie_tx$GS689 and stringtie_tx$PC3E with tx2gene respectively, by Transcript_ID and t_name
stringtie_tx <- map(stringtie_tx, ~ inner_join(.x, tx2gene, by = c("Transcript_ID" = "t_name")))

# summary appearing times by gene_name, and add this column to stringtie_tx$GS689 and stringtie_tx$PC3E
stringtie_tx <- map(stringtie_tx, ~ mutate(.x, Transcript_Index = as.integer(factor(gene_name, levels = unique(gene_name)))))

#####################
# 02. stringtie with imputation
#####################
stringtie_impute_files <- list.files("data/04.Tx_quantification/02.Stringtie_with_impute", pattern = "t_data.ctab", full.names = TRUE, recursive = TRUE)
tx2gene <- read_tsv(stringtie_impute_files[1]) %>% select(t_name, gene_name, length, num_exons)
stringtie_impute <- list(
    GS689 = tximport(stringtie_impute_files[grepl("GS689", stringtie_impute_files)], type = "stringtie", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, txOut = TRUE)$abundance,
    PC3E = tximport(stringtie_impute_files[grepl("PC3E", stringtie_impute_files)], type = "stringtie", countsFromAbundance = "lengthScaledTPM", tx2gene = tx2gene, txOut = TRUE)$abundance
) %>%
    map(~ as.data.frame(.x) %>%
        rownames_to_column("Transcript_ID") %>%
        as_tibble())
stringtie_impute <- map(stringtie_impute, ~ inner_join(.x, tx2gene, by = c("Transcript_ID" = "t_name")))
stringtie_impute <- map(stringtie_impute, ~ mutate(.x, Transcript_Index = as.integer(factor(gene_name, levels = unique(gene_name)))))
#####################
# 03. Guitar
#####################
# read tsv from all cells like data/04.Tx_quantification/03.Guitar/GS689_1Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt
guitar <- pmap(combinations, function(cells, suffix) {
    read_tsv(paste("data/04.Tx_quantification/03.Guitar/", cells, suffix, "Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt", sep = ""))
})
names(guitar) <- paste(combinations$cells, combinations$suffix, sep = "")

# keep Transcript_ID, and Transcript_Exp columns
# guitar <- map(guitar, ~ select(.x, Transcript_ID, Transcript_Exp, Transcript_length, Transcript_Index, Exons)) %>%
#     map(~ rename(.x, length = Transcript_length, num_exons = Exons))

# # calculate TPM based on Transcript_Exp and Transcript_length
# # calculate RPK
# guitar <- map(guitar, ~ mutate(.x, RPK = Transcript_Exp / (length / 1000)))

# # calculate total RPK
# total_RPK <- sum(sapply(guitar, function(x) sum(x$RPK)))

# # calculate TPM
# guitar <- map(guitar, ~ mutate(.x, TPM = (RPK / total_RPK) * 10^6))
# guitar <- map(guitar, ~ select(.x, Transcript_ID, TPM, length, num_exons, Transcript_Index))

# bind columns of bam_impute respectively for GS689 and PC3E
guitar_combined <- list()
guitar_combined$GS689 <- inner_join(guitar$GS689_1, guitar$GS689_2, by = "Transcript_ID", suffix = c(".1", ".2")) %>%
    inner_join(., guitar$GS689_3, by = "Transcript_ID", suffix = c("", ".3")) %>%
    rename(V1 = Transcript_Exp, V2 = Transcript_Exp.1, V3 = Transcript_Exp.2, length = Transcript_length, num_exons = Exons)
guitar_combined$PC3E <- inner_join(guitar$PC3E_1, guitar$PC3E_2, by = "Transcript_ID", suffix = c(".1", ".2")) %>%
    inner_join(., guitar$PC3E_3, by = "Transcript_ID", suffix = c("", ".3")) %>%
    rename(V1 = Transcript_Exp, V2 = Transcript_Exp.1, V3 = Transcript_Exp.2, length = Transcript_length, num_exons = Exons)

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
calculate_correlations <- function(short_list, name, limits) {
    # generate labels based on limits
    limits_label <- paste0(limits[-length(limits)], "-", limits[-1])

    # rename columns V1, V2, V3 to short_read_1, short_read_2, short_read_3
    short_list <- map(short_list, ~ rename(.x, short_read_1 = V1, short_read_2 = V2, short_read_3 = V3))

    # Combine three_gen and bam_impute$GS689, based on transcript_ID and Transcript_ID, add prefix
    GS689 <- inner_join(three_gen_list$GS689, short_list$GS689, by = c("transcript_ID" = "Transcript_ID"))
    PC3E <- inner_join(three_gen_list$PC3E, short_list$PC3E, by = c("transcript_ID" = "Transcript_ID"))

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
    write_tsv(GS689, paste("data/05.correlations/", name, "_GS689.tsv", sep = ""))
    write_tsv(PC3E, paste("data/05.correlations/", name, "_PC3E.tsv", sep = ""))

    # group by length, to <500, 500-1000, 1000-2000, >2000, break to separate tibbles, and wrap to a list
    GS689_by_length <- GS689 %>%
        mutate(length_group = cut(length, breaks = limits)) %>%
        group_by(length_group) %>%
        nest()
    PC3E_by_length <- PC3E %>%
        mutate(length_group = cut(length, breaks = limits)) %>%
        group_by(length_group) %>%
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
                    ", transcript number = ", nrow(GS689)
                )
            ) +
            ylab(paste("short_read_avg (", name, ")", sep = ""))
        ggsave(paste("data/05.correlations/", sep_method, "/", name, "_vs_3rd_gen_GS689_", limits_label, ".png", sep = ""), GS689_plot, width = 8)
        PC3E_plot <- ggplot(PC3E, aes(x = long_read_avg, y = short_read_avg)) +
            geom_point() +
            geom_smooth(method = "lm", se = FALSE) +
            ggtitle(
                paste(
                    "PC3E, correlation = ", round(PC3E_cor$estimate, 2), ", p-value = ", format(PC3E_cor$p.value, scientific = TRUE),
                    ", transcript number = ", nrow(PC3E)
                )
            ) +
            ylab(paste("short_read_avg (", name, ")", sep = ""))
        ggsave(paste("data/05.correlations/", sep_method, "/", name, "_vs_3rd_gen_PC3E_", limits_label, ".png", sep = ""), PC3E_plot, width = 8)
        # return the correlation coefficient
        c(GS689_cor$estimate, PC3E_cor$estimate)
    }
    single_cor <- pmap(list(GS689_by_length$data, PC3E_by_length$data, limits_label, "sep_by_length"), get_result)
    # names(single_cor) <- c("<1000", "1000-5000", "5000-10000", ">10000")
}
# test
# short_list <- stringtie_tx
# name <- "Stringtie"
limits <- c(0, 500, 1000, 2000, Inf)
limits <- c(0, 1000, 5000, 10000, Inf)
limits_label <- paste0(limits[-length(limits)], "-", limits[-1]) %>% as_factor()
cor_sum <- list()
cor_sum[["Stringtie"]] <- calculate_correlations(stringtie_tx, "Stringtie", limits)
cor_sum[["Stringtie_impute"]] <- calculate_correlations(stringtie_impute, "Stringtie_impute", limits)
cor_sum[["Guitar"]] <- calculate_correlations(guitar_combined, "Guitar", limits)

# use list_cbind to combine all cor_sum
cor_sum_df <- map(cor_sum, ~ map_dfr(.x, ~ as.data.frame(t(.x)))) %>%
    # use map to rename colnames to c("GS689", "PC3E")
    map(~ setNames(.x, c("GS689", "PC3E"))) %>%
    list_cbind() %>%
    unnest(cols = c(Stringtie, Stringtie_impute, Guitar), names_sep = "_") %>%
    # add a col of c("<1000", "1000-5000", "5000-10000", ">10000") to the result
    cbind(length = limits_label) %>%
    # replace $ with _ in colnames
    rename_all(~ str_replace_all(., "\\$", "_"))
write_tsv(cor_sum_df, "data/05.correlations/sep_by_length/cor_sum.tsv")

# use geom_dotplot for cor_sum_df, with x = Stringtie, Stringtie_impute, Guitar, y = value, color = length, facet = cell
p <- cor_sum_df %>%
    gather(key = "method", value = "value", -length) %>%
    # split method to two cols, method and cell, sep by the last _, use str_split_fixed
    mutate(cell = str_extract(method, "(?<=_)[[:alnum:]]+$"), method = str_extract(method, ".*(?=_)")) %>%
    ggplot(aes(x = method, y = value, fill = length)) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
    # add line to connect the dots of the same fill
    geom_line(aes(group = interaction(length, cell), color = length)) +
    # make xlab follow the order of Stringtie, Stringtie_impute, Guitar
    scale_x_discrete(limits = c("Stringtie", "Stringtie_impute", "Guitar")) +
    facet_wrap(~cell) +
    ylab("Spearman correlation coefficient")
ggsave("data/05.correlations/sep_by_length/cor_sum.png", p, width = 8)
