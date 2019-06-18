#this file store codes for feature selection cross-validation  performance with glm 
library (pROC)

glmpred <- function (msubset, X_tr, y_tr, X_ts)
{
	if (length (msubset) == 0) byl_pred <- mean (y_tr - 1) 
	else{
		x_tr_sel <- X_tr[,msubset, drop = F]
		x_ts_sel <- X_ts[,msubset,drop = F]

		byl <- bayesglm(y_tr - 1~., data = data.frame(x=x_tr_sel), family = binomial(link="logit"))

		byl_pred <- predict(byl, newdata= data.frame(x=x_ts_sel), type="response") 
	}
	byl_pred

}

gen_subsets <- function(f_order, no_subsets)
{

  fsubsets <- rep(list(""), no_subsets)

  for(i in 1:no_subsets)
  {
    fsubsets[[i]] <- f_order[1:i]
  }


  fsubsets
}



############################################################


cv_glm <- function(fsub, x_total, y_total)
{
#x_total: nxp, y_total: 0 or 1
pred_prob <- NULL

no_folds <- length(y_total)

#print(fsub)
for(i in 1:no_folds){

x_tr_sel <- x_total[-i, fsub]
x_ts_sel <- x_total[i, fsub]

y_tr <- y_total[-i]
y_ts <- y_total[i]




glm <- glm(  y_tr~., data=data.frame(x=x_tr_sel), family=binomial(link="logit") )

glm_pred <- predict(glm, newdata= data.frame(x=x_total[,fsub]), type="response")


temp_prob <- glm_pred[i]


pred_prob <- c(pred_prob,  temp_prob)



}


pred_prob <- cbind (1-pred_prob, pred_prob)

pred_prob

}




cv_bys_glm <- function(fsub, x_total, y_total)
{

if (length (fsub) == 0) pred_prob <- rep (mean (y_total), length (y_total))
else {
library(arm)

#x_total: nxp, y_total: 0 or 1
pred_prob <- NULL
no_folds <- length(y_total)

for(i in 1:no_folds){

#print (fsub)
x_tr_sel <- x_total[-i, fsub]
x_ts_sel <- x_total[i, fsub]

y_tr <- y_total[-i]
y_ts <- y_total[i]


byl <- bayesglm(y_tr~., data = data.frame(x=x_tr_sel), family = binomial(link="logit"))

byl_pred <- predict(byl, newdata= data.frame(x=x_total[,fsub]), type="response")


temp_prob <- byl_pred[i]

pred_prob <- c(pred_prob,  temp_prob)



}

}
pred_prob <- cbind (1-pred_prob, pred_prob)

pred_prob

}





#############################################################


cv_table <- function( fts_array , x_total, y_total , method = "bys_glm")
{

##########################################
#x_total: nxp, y_total: 0 or 1
##########################################
#source("CV.r")

amlp_list <- NULL
er_list   <- NULL
auc_list  <- NULL

top_fsub <- length(fts_array)

table <- data.frame(matrix(0,top_fsub,4))

colnames(table) <- c("fsubsets", "amlp","er", "auc")

# after find interesting feature subsets,
#generate feature subsets prediction table for logistic: 
finish <- 0
pb <- txtProgressBar(min = 0,  max = top_fsub, style = 3)

for(i_fsub in 1:top_fsub){

 fsubset_index <-  fts_array[[i_fsub]]
 #for each feature subset
  
  if(method == "glm")
  {
    pred_prob <- cv_glm(c(fsubset_index), x_total, y_total )
    #x_total: nxp, y_total: 0 or 1
  }

  else{
    pred_prob <- cv_bys_glm(c(fsubset_index), x_total, y_total )
  }


  amlp <- comp_amlp(pred_prob,y_total+1)
  
  er   <- comp_loss(pred_prob, y_total+1)*length(y_total)

  auc <- as.numeric( roc( y_total ~pred_prob[,2])$auc )


  #record all values

  amlp_list <- c(amlp_list, amlp)
  er_list   <- c(er_list,     er)
  auc_list  <- c(auc_list,   auc)

  table[i_fsub,1] <- paste(fsubset_index, collapse =",")
  table[i_fsub,2] <- round(amlp,2)
  table[i_fsub,3] <- er
  table[i_fsub,4] <- round(auc,2)
  
  finish <- finish + 1
  setTxtProgressBar(pb, finish)
}

cat ("\n")

table

}







