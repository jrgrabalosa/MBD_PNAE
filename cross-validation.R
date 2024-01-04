# Cross-validation function ----------------------------------------------------

# Función para calcular la cross-validation, lo que hace es correr el modelo muchas veces para poder calcular estimas de ajuste y poder de predicción de los modelos

cros_val <- function(culex_df, func_f, nam = "variables", n = 10){
  
  # Empty dataframe
  cv_all <- data.frame()
  
  for (i in 1:n){
    
    print(i)
    
    
    # Aquí creas la tabla test y training para el modelo, lo haces para cada iteración
    
    sample <- sample(c(TRUE, FALSE), nrow(culex_df), replace=TRUE, prob=c(0.8, 0.2))
    train  <- culex_df[sample, ]
    test   <- culex_df[!sample, ]
    
    tmb_full <- glmmTMB(func_f,
                        data = train, family = nbinom2) # CUIDADO! Mira si este es el tipo de modelo que utilizas, o si necesitas utilizar otro diferente!
    
    
    pp <- round(
      predict(tmb_full, newdata = test, type = c("response"), allow.new.levels = TRUE),
      0)
    
    
    # Cada estadístico tiene su explicación. OJO! Aquí a lo mejor tienes que modificar lo de females, tienes que poner el nombre de tu variable respuesta
    
    cv <- data.frame(
      model = nam,
      AIC = AIC(tmb_full),
      R2 = R2(pp, test$females), #  Representing the squared correlation between the observed outcome values and the predicted values by the model. The higher the adjusted R2, the better the model.
      RMSE = RMSE(pp, test$females), # which measures the average prediction error made by the model in predicting the outcome for an observation. That is, the average difference between the observed known outcome values and the values predicted by the model. The lower the RMSE, the better the model.
      MAE = MAE(pp, test$females), # an alternative to the RMSE that is less sensitive to outliers. It corresponds to the average absolute difference between observed and predicted outcomes. The lower the MAE, the better the model
      pred.error = RMSE(pp, test$females)/mean(test$females)
    )
    
    cv_all <- rbind(cv_all, cv)
    
  }
  return(cv_all)
}

cv_all <- data.frame()

# Ejemplo --------------------------------------------------------------------

func_f <- females ~ ma_expected + (1 | trap_name)
nam <- c("females ~ ma_expected + (1 | trap_name)")

cv_one <- bind_rows(mclapply(1:10, function(i){ # Aquí estás paralelizando en realidad estas haciendo el modelo 100 iterations = 10*10
  
  print(i)
  
  cros_val(culex_df_ma_agg, func_f, nam = nam, n = 10)
}, mc.cores = 10))  %>%
  group_by(model) %>%
  summarise(
    AIC = mean(AIC),
    R2 = mean(mean(R2)),
    RMSE = mean(RMSE),
    MAE = mean(MAE),
    pred.error = mean(pred.error)
  )

cv_all <- rbind(cv_all, cv_one)