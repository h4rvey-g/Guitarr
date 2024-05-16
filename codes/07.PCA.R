library(tidyverse)
library(psych)
library(Rgb)
library(tximport)
library(scales)
library(ggVennDiagram)
library(UpSetR)
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

# remove rows when V1, V2, V3 are all NA or 0
stringtie_tx <- map(stringtie_tx, ~ filter(.x, rowSums(select(.x, starts_with("V")) == 0) < 3))
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

stringtie_impute <- map(stringtie_impute, ~ filter(.x, rowSums(select(.x, starts_with("V")) == 0) < 3))
#####################
# 03. Guitar
#####################
# read tsv from all cells like data/04.Tx_quantification/03.Guitar/GS689_1Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt
guitar <- pmap(combinations, function(cells, suffix) {
    read_tsv(paste("data/04.Tx_quantification/03.Guitar/", cells, suffix, "Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt", sep = ""))
})
names(guitar) <- paste(combinations$cells, combinations$suffix, sep = "")

# bind columns of bam_impute respectively for GS689 and PC3E
guitar_combined <- list()
guitar_combined$GS689 <- inner_join(guitar$GS689_1, guitar$GS689_2, by = "Transcript_ID", suffix = c(".1", ".2")) %>%
    inner_join(., guitar$GS689_3, by = "Transcript_ID", suffix = c("", ".3")) %>%
    dplyr::rename(V1 = Transcript_Exp.1, V2 = Transcript_Exp.1, V3 = Transcript_Exp.2) %>%
    # remove rows when V1, V2, V3 are all NA or 0
    filter(rowSums(select(., starts_with("V")) == 0) < 3)

guitar_combined$PC3E <- inner_join(guitar$PC3E_1, guitar$PC3E_2, by = "Transcript_ID", suffix = c(".1", ".2")) %>%
    inner_join(., guitar$PC3E_3, by = "Transcript_ID", suffix = c("", ".3")) %>%
    dplyr::rename(V1 = Transcript_Exp.1, V2 = Transcript_Exp.1, V3 = Transcript_Exp.2) %>%
    filter(rowSums(select(., starts_with("V")) == 0) < 3)

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
    select(transcript_ID, starts_with("TPM")) %>%
    filter(rowSums(select(., starts_with("TPM")) == 0) < 3) %>%
    dplyr::rename(Transcript_ID = transcript_ID)

three_gen_list$PC3E <- three_gen_list$PC3E %>%
    mutate(across(starts_with("RPK_"), ~ (. / total_RPK_PC3E) * 10^6, .names = "TPM_{str_remove(.col, 'RPK_')}")) %>%
    select(transcript_ID, starts_with("TPM")) %>%
    filter(rowSums(select(., starts_with("TPM")) == 0) < 3) %>%
    dplyr::rename(Transcript_ID = transcript_ID)

#####################
# 05. PCA
#####################
library(ggfortify)
stringtie_tx_pca_data <- stringtie_tx %>%
    map(~ select(.x, starts_with("V"))) %>%
    imap(~ mutate(.x, cell = .y)) %>%
    reduce(bind_rows)
stringtie_tx_pca <- prcomp(stringtie_tx_pca_data %>% select(-cell), scale = TRUE)
p <- autoplot(stringtie_tx_pca, data = stringtie_tx_pca_data, colour = "cell")
ggsave("data/07.PCA/stringtie_tx_pca.png", p, width = 10, height = 10, dpi = 300)

stringtie_impute_pca_data <- stringtie_impute %>%
    map(~ select(.x, starts_with("V"))) %>%
    imap(~ mutate(.x, cell = .y)) %>%
    reduce(bind_rows)
stringtie_impute_pca <- prcomp(stringtie_impute_pca_data %>% select(-cell), scale = TRUE)
p <- autoplot(stringtie_impute_pca, data = stringtie_impute_pca_data, colour = "cell")
ggsave("data/07.PCA/stringtie_impute_pca.png", p, width = 10, height = 10, dpi = 300)

guitar_combined_pca_data <- guitar_combined %>%
    map(~ select(.x, starts_with("V"))) %>%
    imap(~ mutate(.x, cell = .y)) %>%
    reduce(bind_rows)
guitar_combined_pca <- prcomp(guitar_combined_pca_data %>% select(-cell), scale = TRUE)
p <- autoplot(guitar_combined_pca, data = guitar_combined_pca_data, colour = "cell")
ggsave("data/07.PCA/guitar_combined_pca.png", p, width = 10, height = 10, dpi = 300)

three_gen_list_pca_data <- three_gen_list %>%
    map(~ select(.x, starts_with("TPM"))) %>%
    imap(~ mutate(.x, cell = .y)) %>%
    reduce(bind_rows)
three_gen_list_pca <- prcomp(three_gen_list_pca_data %>% select(-cell), scale = TRUE)
p <- autoplot(three_gen_list_pca, data = three_gen_list_pca_data, colour = "cell")
ggsave("data/07.PCA/three_gen_list_pca.png", p, width = 10, height = 10, dpi = 300)