library(glmnet)   #这个包用于二元分类任务的数据准备
load("step2.Rdata")


train_data <- t(cpg1[g11,]) #
train_y <- ifelse(group_cpg1=="Normal",0,1) #
test_data <- t(cpg2[g11,])  #
test_data1 <- t(dia1[g11,]) 
test_y1 <- ifelse(group_dia1=="normal",0,1)  #
test_y <- ifelse(group_cpg2=="normal",0,1)

set.seed(1)  #这一行设置随机数生成器的种子，使得每次运行代码得到的结果都是一样的
library(glmnet)

model_lasso <- glmnet(train_data,train_y,family = "binomial", maxit = 10000)
print(model_lasso)   
cv_fit <- cv.glmnet(train_data,train_y, nlambda = 1000,alpha = 1)
plot(cv_fit)  
fit <- glmnet(train_data,train_y,family = "binomial", maxit = 10000) #
plot(fit, xvar = "lambda", label = TRUE)  #
model_lasso_min <- glmnet(train_data,train_y, alpha = 1, lambda=cv_fit$lambda.min)   #
model_lasso_1se <- glmnet(train_data,train_y, alpha = 1, lambda=cv_fit$lambda.1se)   #
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]  #
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]  #
length(choose_gene_min)     #
length(choose_gene_1se)




lasso.prob <- predict(cv_fit, newx=train_data , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(train_y ,lasso.prob)  #
re=as.data.frame(re)   #
colnames(re)=c('event','prob_min','prob_1se')  #
re$event=as.factor(re$event)   #
library(ggpubr) 
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
library(patchwork)
p1+p2
library(ROCR)
library(caret)
#min
pred_min <- prediction(re[,2], re[,1])  #
auc_min = performance(pred_min,"auc")@y.values[[1]] #
perf_min <- performance(pred_min,"tpr","fpr") #
plot(perf_min,colorize=FALSE, col="blue")  #
lines(c(0,1),c(0,1),col = "gray", lty = 4 ) # 
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))  #
#1se
pred_1se <- prediction(re[,3], re[,1]) #
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_1se,3)))



lasso.prob <- predict(cv_fit, newx=test_data , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re1=cbind(test_y ,lasso.prob)
re1=as.data.frame(re)
colnames(re1)=c('event','prob_min','prob_1se')
re1$event=as.factor(re1$event)
library(ggpubr) 
library(stringr)
re1$event <- ifelse(str_detect(re1$event,"0"),"CON","ACP")   

p1 = ggboxplot(re1, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re1, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
library(patchwork)
p1
p2
library(ROCR)
library(caret)
#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_1se,3)))


lasso.prob <- predict(cv_fit, newx=test_data1 , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(test_y1 ,lasso.prob)
re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
re$event=as.factor(re$event)
library(ggpubr) 
library(stringr)
re1 <- re
re1$event <- ifelse(str_detect(re1$event,"0"),"CON","DIA")
p1 = ggboxplot(re1, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re1, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
library(patchwork)
p1
p2
library(ROCR)
library(caret)

#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_1se,3)))


save(choose_gene_1se,file = "choose_gene_1se.Rdata")
