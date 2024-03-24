library(tidyverse)
library(psych)
# read tsv test/3rd_gen/3rd_gen.txt, and drop columns begin with HEK293T
three_gen <- read_tsv("test/3rd_gen/3rd_gen.txt") %>% select(-starts_with("HEK293T"))
cells <- c("GS689", "PC3E")
# read tsv from all cells like test/bam_impute/GS689_1Aligned.sortedByCoord.out/GS689_1Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt
combinations <- expand.grid(cells = cells, suffix = c("_1", "_2", "_3"), stringsAsFactors = FALSE)
bam_impute <- pmap(combinations, function(cells, suffix) {
  read_tsv(paste("test/bam_impute/", cells, suffix, "Aligned.sortedByCoord.out/", cells, suffix, "Aligned.sortedByCoord.out.readsImpute.transcript.expression.txt", sep = ""))
})
names(bam_impute) <- paste(combinations$cells, combinations$suffix, sep = "")
# keep Gene_ID, Transcript_ID, and Transcript_Exp columns
bam_impute <- map(bam_impute, ~ select(.x, Gene_ID, Transcript_ID, Transcript_Exp))
# bind columns of bam_impute respectively for GS689 and PC3E
bam_impute_combined <- list()
bam_impute_combined$GS689 <- inner_join(bam_impute$GS689_1, bam_impute$GS689_2, by = "Transcript_ID", suffix = c(".1", ".2")) %>%
    inner_join(.,bam_impute$GS689_3, by = "Transcript_ID", suffix = c("", ".3")) %>%
    select(-starts_with("Gene_ID"))
bam_impute_combined$PC3E <- inner_join(bam_impute$PC3E_1, bam_impute$PC3E_2, by = "Transcript_ID", suffix = c(".1", ".2")) %>%
    inner_join(.,bam_impute$PC3E_3, by = "Transcript_ID", suffix = c("", ".3")) %>%
    select(-starts_with("Gene_ID"))
# rename columns of bam_impute_combined from Transcript_Exp to short_read_1, short_read_2, short_read_3
bam_impute_combined$GS689 <- bam_impute_combined$GS689 %>% rename(short_read_1 = Transcript_Exp.1, short_read_2 = Transcript_Exp.2, short_read_3 = Transcript_Exp)
bam_impute_combined$PC3E <- bam_impute_combined$PC3E %>% rename(short_read_1 = Transcript_Exp.1, short_read_2 = Transcript_Exp.2, short_read_3 = Transcript_Exp)
# remove _[digit] from Transcript_ID
three_gen$transcript_ID <- three_gen$transcript_ID %>% str_replace("_\\d+$", "")
# split three_gen to a list of two data frames,based on column names containing GS689 or PC3E
three_gen_list <- list()
three_gen_list$GS689 <- three_gen %>% select(transcript_ID, starts_with("GS689")) %>% rename(long_read_1 = GS689_1, long_read_2 = GS689_2, long_read_3 = GS689_3)
three_gen_list$PC3E <- three_gen %>% select(transcript_ID, starts_with("PC3E")) %>% rename(long_read_1 = PC3E_1, long_read_2 = PC3E_2, long_read_3 = PC3E_3)
# Combine three_gen and bam_impute$GS689, based on transcript_ID and Transcript_ID, add prefix
gs689 <- inner_join(three_gen_list$GS689, bam_impute_combined$GS689, by = c("transcript_ID" = "Transcript_ID"))
PC3E <- inner_join(three_gen_list$PC3E, bam_impute_combined$PC3E, by = c("transcript_ID" = "Transcript_ID"))
calc_cor <- function(row) {
  row2 <- print(as.numeric(row[2:7]))
  cor_result <- corr.test(row2[1:3], row2[4:6], method = "spearman")
}
# use apply to calculate correlation for each row of gs689
gs689_cor <- apply(gs689, 1, calc_cor)
# take the correlation result from the list, and draw a density plot
gs689_cor <- map_dbl(gs689_cor, ~ .x$r)
p <- ggplot(data.frame(gs689_cor = gs689_cor), aes(x = gs689_cor)) + geom_density(fill = "blue", alpha = 0.5) + labs(title = "GS689 correlation density plot", x = "correlation", y = "density")
ggsave("gs689_cor.png", p)
PC3E_cor <- apply(PC3E, 1, calc_cor)
