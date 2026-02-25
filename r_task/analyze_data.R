data <- read.csv("sample_data.csv")
mean_score <- mean(data$Score)
cat("Среднее значение показателя Score", mean_score, "\n")
treatment_data <- data[data$Group == "Treatment", ]
max_tr_score <- max(treatment_data$Score)
cat("Максимальное значение Score в группе Treatment", max_tr_score, "\n")
png("score_boxplot.png", width=800, height=600)
boxplot(Score ~ Group, data=data, main="Score Distribution by Group", xlab="Group", ylab="Score")
dev.off()

