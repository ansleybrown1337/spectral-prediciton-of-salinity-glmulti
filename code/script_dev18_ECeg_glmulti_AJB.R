# GLmulti Learning, created by Brian Craig and A.J. Brown in Nov. 2022

# Script developed for use with data frames containing spectral indices to 
# develop a generalized multiple linear regression (GLM) model to predict 
# EC. The model is then evaluated
# and cross-validated via "leave one field out" root mean squared prediction
# Error.

#For glmulti, we're looking at 2018 Dev ECe NOT interpolated

library(MuMIn)
library(stats)
library(hydroGOF)
library(ggplot2)
library(caret)
library(tdr)
library(glmulti)
#library(performance)


# Step 0: Importing and cleaning data
# import data file here (in this case "all2018.csv")
original_df = read.csv(file.choose(), na.strings = c(-9999, 'NaN'))
# remove rows with N/As and -9999s
df = original_df[complete.cases(original_df),]
# change na. action; dredge will now fail if n/a
options(na.action = 'na.fail')
# clean up Field column so that fields are correctly labeled
levels(df$Field)[levels(df$Field) == 2] <- 'Muth2'
levels(df$Field)[levels(df$Field) == 3] <- 'Muth3'
levels(df$Field)[levels(df$Field) == 8] <- 'Muth8'

# Create log transformations of values
df$ln_ECeg18 = log(df$ECeg18)
df$ln_ND = log(df$ND)
df$ln_CD = log(df$CD)
df$ln_ED = log(df$ED)
df$ln_CRD = log(df$CRD)
df$ln_CCD = log(df$CCD)
#df$ln_CRDF = log(df$CRDF)
#df$ln_ETDF = log(df$ETDF)


# Step 1: Use glmulti to create best linear model

  glmulti("ln_ECeg18",c("ND","CD","ED","ln_ND","ln_CD","ln_ED"),df,level=2,method="d",crit="aicc",maxsize = 8)
  GLresults <- glmulti("ln_ECeg18",c("ND","CD","ED","ln_ND","ln_CD","ln_ED"),df,level=2,method="g",crit="aicc",maxsize=8)

  #save(GLresults,file = "Q:/rcode/Model/2018/ResultsTemp/glmulti_dev18ECeg_maxsize8.r")

  final_model <- GLresults@objects[[1]]
  second_model <- GLresults@objects[[2]]
  print(final_model)
  print(final_model$formula)

  summary(final_model)

  # the below code is not needed when hydroGOF works, which it should
  #compare_performance(final_model,second_model)

  # Check model assumptions
  # Coefficients and significance
  # Summary plots
  par(mfrow=c(2,2))
  plot(final_model)

  #This is the end of the glmulti section
  ##############################################################################

# Step 2: Take best model and create output predictions
  df$ln_ECeg18_pred = predict(final_model, newdata=df)
  df$ECeg18_pred = exp(df$ln_ECeg18_pred)  #back-transform from log
  # df with only muth2 rows
  m2df = subset(df, Field == 'Muth2')

# Step 3: Run GOF stats on predictions
  GOF_stats = function(df, ECeg18=T) {
    # GOF stats generated via functions from hydroGOF package:
    # https://www.rforge.net/doc/packages/hydroGOF/d.html
    # ifelse statement to select the prediction of ECeg or not
    if(ECeg18==T) {
      obs = df$ECeg18
      pred = df$ECeg18_pred
    } else {
      obs = df$ECe
      pred = df$ECe_pred
    }
    # Root Mean Squared Error (RMSE)
    RMSE_model = rmse(pred, obs)
    # Normalized RMSE (NRMSE); normalized by std. dev. by default
    NRMSE_model = nrmse(pred, obs)
    # Pearson's Correlation Statistic
    Cor_model = cor(pred, obs)
    # Willmott's index of agreement (IOA)
    IOA_model = d(pred, obs)
    
    print("full dataset statistics:###########################################")
    Summary = data.frame(
      Statistic = c('RMSE', 'NRMSE', 'Corr', 'IOA'),
      Value = c(RMSE_model, NRMSE_model, Cor_model, IOA_model)
    )
    print(Summary)
  }
  GOF_stats(df=df, ECeg18 = T) # get GOF stats on all fields
  GOF_stats(df=m2df,ECeg18 = T) # get GOF stats on just M2

  #Calculate the Mean Bias Error
  #MBE = tdStats(df$ECe18, df$ECe18_pred, functions = c("mbe"))
  MBE = tdStats(df$ECeg18, df$ECeg18_pred, functions = c("mbe"))
  MBE
  MBEm2 = tdStats(m2df$ECeg18, m2df$ECeg18_pred, functions = c("mbe"))
  MBEm2 # MBE for just Muth2 field

# Step 4: Cross-validate output predictions 
  #Step 4a: using Leave-one-field-out (LOFO) Root Mean Squared Prediction Error (RMSPE)
    lofo = function(df, fxn){
      rmspe_list = c()
      field_list = c()
      for(i in unique(df$Field)){
        test = subset(df, Field != i)
        lofo = subset(df, Field == i)
        # extract formula and apply to new data
        mdl = fxn
        y_pred = exp(predict(mdl, newdata = lofo)) #back-transform from log
        rmspe = sqrt(mean((lofo$ECeg18 - y_pred)^2))
        rmspe_list = append(rmspe_list, rmspe)
        field_list = append(field_list, i)
      }
      final_df = data.frame(Left_out_Field = field_list,
                            ECeg18_RMSPE_dSm = rmspe_list)
      print(final_df)
      print('Average RMSPE:')
      print(mean(final_df$ECeg18_RMSPE_dSm))
    }
    lofo(df = df, fxn = final_model)

  # Step 4b: Using repeated K-folds
    set.seed(123)
    train.control = trainControl(method = "repeatedcv", number = 10, repeats = 100)
    # Train the model
    model = train(final_model$formula, data = df, method = "lm",
                  trControl = train.control)
    # Summarize the results
    print(model)
    exp(model$results[2])

# Step 5 export predictions and original dataframe to CSV
  # Get working directory (i.e. where the CSV will be saved)
  getwd()
  # Save CSV: predicted ECe values will be a column named "ECe_pred"
  write.csv(df, paste('Q:/rcode/Model/2018/ResultsTemp','/result.csv', sep = ''))

# Step 6 (optional) 1:1 plot of results
  par(pty="s")
  qplot(ECeg18, ECeg18_pred, data = df, geom = 'point', color = Field) +
    geom_abline(slope = 1) +
    coord_fixed(xlim = c(0, max(df$ECeg18)), ylim = c(0, max(df$ECeg18))) +
    #coord_fixed(xlim = c(0, max(df$ECeg18)), ylim = c(0, 40)) +
    ggtitle("1:1 Plot of Observed and Model Predicted Vales \n2018 Dev, ECeg = f(NDVI, CRSI, ETD), log-log, No Interpo") +
    xlab(expression(Observed~EC[eg]~dS~m^{-1})) +
    ylab(expression(Predicted~EC[eg]~dS~m^{-1})) +
    theme(plot.title = element_text(hjust = 0.5))

  ggplot(df, aes(x=ECeg18, y=ECeg18_pred) ) +
    geom_bin2d() +
    geom_abline(slope = 1) +
    coord_fixed(xlim = c(0, max(df$ECeg18)), ylim = c(0, max(df$ECeg18))) +
    ggtitle("Dev ECeg18 = f(NDVI, CRSI, ED), log-log, no interpo") +
    xlab(expression(Observed~EC[eg]~dS~m^{-1})) +
    ylab(expression(Predicted~EC[eg]~dS~m^{-1})) +
    theme(plot.title = element_text(hjust = 0.5))

