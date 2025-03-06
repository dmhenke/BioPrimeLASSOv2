# GeneRIF
# doc2vec exploring

# LIBRARIES ####
lapply(c("data.table","biomaRt","doc2vec"),require, character.only = TRUE)





# DATA IN ####
## Gene RIF ####
generif <- fread("./Data/generifs_basic.gz",sep='\t')


# Create doc2vec model ####
# texttextvec <- txt_clean_word2vec(
#   generif$terms, 
#   ascii = TRUE, alpha = TRUE, tolower = TRUE, trim = TRUE)
# db <- data.frame(doc_id = generif$gene, text = text) data.frame(doc_id = generif$gene, text = text)
# clean input text
df_d2v <- generif[,c("Gene ID","GeneRIF text")];df_d2v$`GeneRIF text` <- gsub('\n','\t',df_d2v$`GeneRIF text`);colnames(df_d2v) <- c("doc_id","text")
df_d2v <- df_d2v[!duplicated(df_d2v),]
df_d2v <- df_d2v[grep("observational study of gene disease",df_d2v$text,invert = T),]

df_d2vsplt <- do.call(rbind,lapply(split(df_d2v,df_d2v$doc_id),function(x){
  
  
}))


df_d2v$text   <- tolower(df_d2v$text)
df_d2v$text   <- gsub("[^[:alpha:]]", " ", df_d2v$text)
df_d2v$text   <- gsub("[[:space:]]+", " ", df_d2v$text)
df_d2v$text   <- trimws(df_d2v$text)
df_d2v$nwords <- txt_count_words (df_d2v$text)
df_d2v <- subset(df_d2v, nwords < 1000 & nchar(text) > 0)

## Build the model
# system.time(model <- paragraph2vec(x = df_d2v[1:30000,], type = "PV-DM",   dim = 15,  iter = 5))
system.time(model <- paragraph2vec(x = df_d2v, type = "PV-DM",   dim = 15,  iter = 5)) # elapsed tiem 56 seconds
str(model)
embedding <- as.matrix(model, which = "words")
embedding <- as.matrix(model, which = "docs")
head(embedding)
## Get vocabulary
vocab <- summary(model, type = "vocabulary",  which = "docs")
vocab <- summary(model, type = "vocabulary",  which = "words")
