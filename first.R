################################### creating table ###################################################

load("expression.cb1.rda")
load("expression.cb2.rda")
load("clinical.cb.rda")

expression <- rbind(expression.cb1, expression.cb2)
expression[1:5,1:5]
clinical.cb[,1]
clinical.cb[,1] <- gsub("\\-", "\\.", clinical.cb[,1])
clinical.cb[,1]

expression_tmp <- t(expression[,-1])
expression_tmp <- as.data.frame(expression_tmp)
expression[1:5,1:5]
expression_tmp[1:5,1:5]
colnames(expression_tmp) <- expression[,1]
expression[1:5,1:5]
expression_tmp[1:5,1:5]

expression_tmp <- cbind(rownames(expression_tmp), expression_tmp)
expression_tmp[1:5,1:5]

dim(expression_tmp)
dim(expression)

colnames(expression_tmp)[1] <- "PatientID"
colnames(clinical.cb)[1] <- "PatientID"

expression_tmp[1:5,1:5]
clinical.cb[1:5,1:5]


merged_table <- merge(expression_tmp, clinical.cb[,c("PatientID", "X_cohort")], by="PatientID")
merged_table[1:5, 1:5]
merged_table[1:5, 16110:16117]

colnames(merged_table) <- gsub("\\?\\|","G", colnames(merged_table))
colnames(merged_table) <- paste("G", colnames(merged_table), sep="")
colnames(merged_table) <- gsub("\\-","", colnames(merged_table))
colnames(merged_table) <- gsub("\\,","", colnames(merged_table))
colnames(merged_table) <- gsub(" ","", colnames(merged_table))


merged_table[1:5, 1:5]
merged_table[1:5, 16110:16117]
dim(merged_table)
merged_table[,"GX_cohort"] <- as.factor(merged_table[,"GX_cohort"])


form <- paste(colnames(merged_table)[16110],"GX_cohort", sep="~")
r <- aov(as.formula(form), data=merged_table)

######################################################################################################################



########################################## testing ##################################################################

pval_test <- summary(aov(as.formula(form), data=merged_table))[[1]][["Pr(>F)"]]
summary(r)
merged_table_2 <- merged_table[,-4442]

########################################## statistics ##############################################################

range(merged_table[,2:(ncol(merged_table)-1)], na.rm = TRUE)
plot(by(merged_table[,100], merged_table$GX_cohort, mean))

boxplot(change.2191~Mutation, data=AML)
abline(h=0, col='grey', lwd=2)

####################################################################################################################
small_gens_vector <- 2:(ncol(small_merged_table)-1)
scaled_merged_table = cbind(merged_table$GPatientID, merged_table[,small_gens_vector]*100, merged_table$GX_cohort)

n <- 500
small_merged_table = merged_table[,c(1:n, ncol(merged_table))]

p_values <- unlist(lapply(c(2:(ncol(small_merged_table)-1)), function(i) {
  summary(aov(as.formula(paste(colnames(small_merged_table)[i],"GX_cohort", sep="~")), data=small_merged_table))[[1]][["Pr(>F)"]]
}
))

p_values <- unlist(lapply(c(2:(ncol(scaled_merged_table)-1)), function(i) {
  summary(aov(as.formula(paste(colnames(scaled_merged_table)[i],"merged_table$GX_cohort", sep="~")), data=scaled_merged_table))[[1]][["Pr(>F)"]][1]
}
))

hist(p_values*10e50)

length(which(p_values==0))
sort(p_values, decreasing = FALSE)

p_values <- unlist(lapply(c(2:(ncol(merged_table)-1)), function(i) {
  aov(as.formula(paste(colnames(merged_table)[i],"GX_cohort", sep="~")), data=merged_table)
}
))


p_values <- lapply(2:(ncol(merged_table)-1), function(i) {
  aov(as.formula(paste(colnames(merged_table)[i],"GX_cohort", sep="~")), data=merged_table)
}
)

###############################liczenie p-warto########################################
formulas <- lapply(2:(ncol(merged_table)-1), function(i) {
  as.formula(paste(colnames(merged_table)[i],"GGX_cohort", sep="~"))
}
)

p_values <- unlist(lapply(formulas, function(i) {
  summary(aov(i, data=merged_table))[[1]][["Pr(>F)"]][1]
}
))

#########################################################################################
pairwise.t.test("GGG10431", "GGX_cohort", data=merged_table, na.omit=T)





#########################################################################################

summary(aov(formulas[[2]], data=merged_table))[[1]][["Pr(>F)"]][1]

p_values <- unlist(lapply(2:(ncol(merged_table)-1), function(i) {
  aov(as.formula(paste(colnames(merged_table)[i],"GX_cohort", sep="~")), data=merged_table)
    }
  ))


