load("C:/Users/po320237/Downloads/expression.cb1.rda")
load("C:/Users/po320237/Downloads/expression.cb2.rda")
load("C:/Users/po320237/Downloads/clinical.cb.rda")

expression <- rbind(expression.cb1, expression.cb2)
colnames(expression)
clinical.cb[,1]
clinical.cb[,1] <- gsub("\\-", "\\.", clinical.cb[,1])

install.packages("dplyr")
library(dplyr)

class(expression)
class(clinical.cb)

clinical.cb <- as.matrix(clinical.cb)
expression <- as.matrix(expression)
expression_tmp <- t(rbind(colnames(expression), expression))

dim(expression_tmp)
dim(expression)
colnames(expression_tmp)[1] <- "PatientID"
colnames(clinical.cb)[1] <- "PatientID"

colnames(expression_tmp) <- expression_tmp[1,]
expression_tmp <- expression_tmp[-1,]

merged_table <- merge(expression_tmp, clinical.cb[,c("PatientID", "X_cohort")], by="PatientID")

colnames(merged_table) <- gsub("\\?\\|","G", colnames(merged_table))
colnames(merged_table) <- paste("G", colnames(merged_table), sep="")
dim(merged_table)
gen_count <- dim(merged_table)[2]-2
merged_table <- as.factor(merged_table)
form <- paste(colnames(merged_table)[2],"X_cohort", sep="~")
aov(as.formula(form), data=merged_table)
aov(as.factor(G100133144)~as.factor(GX_cohort), data=merged_table)
p_values <- unlist(lapply(1:gen_count, function(i) {aov(colnames(merged_table)[i+1]~X_cohort, data=merged_table}))

