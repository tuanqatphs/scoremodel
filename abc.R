# Predictions

epl_predict = function(data_df, played_matches, upcoming_matches){
  
  implied_prob = implied_probabilities(data_df[, c('PSCH', 'PSCD', 'PSCA')], 
                                       method = 'basic')
  implied_expg = expg_from_probabilities(implied_prob$probabilities)
  data_df$xhg = implied_expg$expg[, 1]
  data_df$xag = implied_expg$expg[, 2]
  
  data_df$expg1 = rep(NA, nrow(data_df))
  data_df$expg2 = rep(NA, nrow(data_df))
  
  new_train_set = data_df[1:played_matches, ]
  
  for (i in 1:upcoming_matches){
    
    my_weights = weights_dc(as.Date(new_train_set$Date, format = '%d/%m/%Y'), xi=0.05)
    
    model = goalmodel(goals1 = round(new_train_set$xhg, 0), 
                      goals2 = round(new_train_set$xag, 0),
                      team1 = new_train_set$HomeTeam, 
                      team2 = new_train_set$AwayTeam,
                      dc = F, model = 'gaussian',
                      weights = my_weights)
    
    new_match = data_df[played_matches+i, ]
    
    new_pred = predict_expg(model, return_df = TRUE, 
                            team1 = new_match$HomeTeam, 
                            team2 = new_match$AwayTeam)
    
    data_df$expg1[played_matches+i] = new_pred$expg1
    data_df$expg2[played_matches+i] = new_pred$expg2
    
    new_train_set = rbind(new_train_set, new_match)
    
  }
  
  a = runif(1, min = 0, max = 0.5)
  b = 1 - a
  
  weight_rounding = round(round(data_df[, c('expg1', 'expg2')], 10) * a + 
                            round(data_df[, c('expg1', 'expg2')], 10) * b, 2)
  
  data_df$expg1 = weight_rounding$expg1
  data_df$expg2 = weight_rounding$expg2
  
  return(data_df)
  
}

epl_matchday_predict = function(data_df, played_matches, upcoming_match_days){
  
  implied_prob = implied_probabilities(data_df[, c('AvgCH', 'AvgCD', 'AvgCA')], 
                                       method = 'basic')
  implied_expg = expg_from_probabilities(implied_prob$probabilities)
  data_df$xhg = implied_expg$expg[, 1]
  data_df$xag = implied_expg$expg[, 2]
  
  data_df$expg1 = rep(NA, nrow(data_df))
  data_df$expg2 = rep(NA, nrow(data_df))
  
  new_train_set = data_df[1:played_matches, ]
  
  for (i in 1:upcoming_match_days){
    
    my_weights = weights_dc(as.Date(new_train_set$Date, format = '%d/%m/%Y'), xi=0.009*i)
    
    model = goalmodel(goals1 = round(new_train_set$xhg), 
                      goals2 = round(new_train_set$xag),
                      team1 = new_train_set$HomeTeam, 
                      team2 = new_train_set$AwayTeam,
                      dc = F, model = 'gaussian',
                      weights = NULL)
    
    lower_bound = played_matches + 1 + (i-1)*10
    upper_bound = played_matches + i*10
    
    new_match_day = data_df[lower_bound:upper_bound, ]
    
    new_pred = predict_expg(model, return_df = TRUE, 
                            team1 = new_match_day$HomeTeam, 
                            team2 = new_match_day$AwayTeam)
    
    data_df$expg1[lower_bound:upper_bound] = new_pred$expg1
    data_df$expg2[lower_bound:upper_bound] = new_pred$expg2
    
    new_train_set = rbind(new_train_set, new_match_day)
    
  }
  
  return(data_df)
  
}

epl_matchday_remodelling = function(data_df, played_matches, upcoming_match_days){
  
  implied_prob = implied_probabilities(data_df[, c('PSCH', 'PSCD', 'PSCA')], 
                                       method = 'basic')
  implied_expg = expg_from_probabilities(implied_prob$probabilities)
  data_df$xhg = implied_expg$expg[, 1]
  data_df$xag = implied_expg$expg[, 2]
  
  data_df$expg1 = rep(NA, nrow(data_df))
  data_df$expg2 = rep(NA, nrow(data_df))
  
  new_train_set = data_df[1:played_matches, ]
  
  for (i in 1:upcoming_match_days){
    
    my_weights = weights_dc(as.Date(new_train_set$Date, format = '%d/%m/%Y'), xi=0.05)
    
    model = goalmodel(goals1 = round(new_train_set$xhg, 0), 
                      goals2 = round(new_train_set$xag, 0),
                      team1 = new_train_set$HomeTeam, 
                      team2 = new_train_set$AwayTeam,
                      dc = F, model = 'gaussian',
                      weights = my_weights)
    
    model_params = data.frame(model$parameters)
    model_params = cbind(rownames(model_params), rownames(model_params), 
                         model_params)
    rownames(model_params) = NULL
    colnames(model_params) = c('HomeTeam', 'AwayTeam', 'attack', 'defense', 'hfa')
    
    lower_bound = played_matches + 1 + (i-1)*10
    upper_bound = played_matches + i*10
    
    new_match_day = data_df[lower_bound:upper_bound, ]
    new_train_set = rbind(new_train_set, new_match_day)
    
    bbb = new_train_set
    
    bbb$h_att = rep(NA, nrow(bbb))
    bbb$a_def = rep(NA, nrow(bbb))
    bbb$a_att = rep(NA, nrow(bbb))
    bbb$h_def = rep(NA, nrow(bbb))
    
    for (j in 1:nrow(bbb)){
      
      h_att_index = which(bbb$HomeTeam[j] == model_params$HomeTeam)
      bbb$h_att[j] = model_params$attack[h_att_index]
      
      a_def_index = which(bbb$AwayTeam[j] == model_params$AwayTeam)
      bbb$a_def[j] = model_params$defense[a_def_index]
      
      a_att_index = which(bbb$AwayTeam[j] == model_params$AwayTeam)
      bbb$a_att[j] = model_params$attack[a_att_index]
      
      h_def_index = which(bbb$HomeTeam[j] == model_params$HomeTeam)
      bbb$h_def[j] = model_params$defense[h_def_index]
    }
    
    # train a new model
    
    bbb_home_model = glm(xhg ~ h_att + a_def, 
                         family = gaussian,
                         data = bbb[1:(nrow(bbb)-10), ])
    pred_fthg = predict(bbb_home_model, tail(bbb, 10))
    
    bbb_away_model = glm(xag ~ a_att + h_def,
                         family = gaussian,
                         data = bbb[1:(nrow(bbb)-10), ])
    pred_ftag = predict(bbb_away_model, tail(bbb, 10))
    
    # append the predictions
    
    data_df$expg1[lower_bound:upper_bound] = pred_fthg
    data_df$expg2[lower_bound:upper_bound] = pred_ftag
    
  }
  
  #a = runif(1, min = 0, max = 0.5)
  #b = 1 - a
  #
  #weight_rounding = round(round(data_df[, c('expg1', 'expg2')], 10) * a + 
  #                          round(data_df[, c('expg1', 'expg2')], 10) * b, 2)
  #
  #data_df$expg1 = weight_rounding$expg1
  #data_df$expg2 = weight_rounding$expg2
  
  return(data_df)
  
}

epl_nn = function(data_df, played_matches, upcoming_match_days){
  
  #data_df$FTHG[data_df$FTHG > 5] = 5
  #data_df$FTAG[data_df$FTAG > 5] = 5
  
  implied_prob = implied_probabilities(data_df[, c('AvgCH', 'AvgCD', 'AvgCA')], 
                                       method = 'basic')
  implied_expg = expg_from_probabilities(implied_prob$probabilities)
  data_df$xhg = implied_expg$expg[, 1]
  data_df$xag = implied_expg$expg[, 2]
  
  data_df$expg1 = rep(NA, nrow(data_df))
  data_df$expg2 = rep(NA, nrow(data_df))
  
  new_train_set = data_df[1:played_matches, ]
  
  for (i in 1:upcoming_match_days){
    
    lower_bound = played_matches + 1 + (i-1)*10
    upper_bound = played_matches + i*10
    
    new_match_day = data_df[lower_bound:upper_bound, ]
    
    bbb_home_model = neuralnet(FTHG ~ xhg + xag, hidden = c(5, 5),
                               data = new_train_set,
                               threshold = 1,
                               stepmax = 1e7)
    pred_fthg = predict(bbb_home_model, new_match_day)
    
    bbb_away_model = neuralnet(FTAG ~ xhg + xag, hidden = c(5, 5),
                               data = new_train_set,
                               threshold = 1,
                               stepmax = 1e7)
    pred_ftag = predict(bbb_away_model, new_match_day)
    
    data_df$expg1[lower_bound:upper_bound] = pred_fthg
    data_df$expg2[lower_bound:upper_bound] = pred_ftag
    
    new_train_set = rbind(new_train_set, new_match_day)
    
  }
  
  return(data_df)
  
}

# OU

calculate_ou_00 = function(expg1, expg2, over_under){
  
  stopifnot(over_under%%1 == 0)
  
  all_possibilities = as.matrix(dpois(0:49, expg1) %o% dpois(0:49, expg2))
  total_goals = matrix(NA, nrow = 50, ncol = 50)
  
  for (i in 1:50){
    for (j in 1:50){
      total_goals[i, j] = (i-1) + (j-1)
    }
  }
  
  over_prob = sum(all_possibilities[total_goals > over_under])
  equal_prob = sum(all_possibilities[total_goals == over_under])
  under_prob = sum(all_possibilities[total_goals < over_under])
  
  print(over_prob + equal_prob + under_prob)
  
  new_over_prob = over_prob / (over_prob + under_prob)
  new_under_prob = under_prob / (over_prob + under_prob)
  
  over_odds = 1/new_over_prob
  under_odds = 1/new_under_prob
  
  return(c(over_odds, under_odds))
  
}

calculate_ou_25 = function(expg1, expg2, over_under){
  
  stopifnot(over_under%%1 == 0.25,
            over_under > 0)
  
  all_possibilities = as.matrix(dpois(0:49, expg1) %o% dpois(0:49, expg2))
  total_goals = matrix(NA, nrow = 50, ncol = 50)
  
  for (i in 1:50){
    for (j in 1:50){
      total_goals[i, j] = (i-1) + (j-1)
    }
  }
  
  over_odds_00 = calculate_ou_00(expg1, expg2, over_under - 0.25)[1]
  over_odds_05 = calculate_ou_05(expg1, expg2, over_under + 0.25)[1]
  
  over_odds = (over_odds_00 + over_odds_05)/2
  over_prob = 1/over_odds
  under_prob = 1 - over_prob
  under_odds = 1/under_prob
  
  return(c(over_odds, under_odds))
  
}

calculate_ou_05 = function(expg1, expg2, over_under){
  
  stopifnot(over_under%%1 == 0.5)
  
  all_possibilities = as.matrix(dpois(0:49, expg1) %o% dpois(0:49, expg2))
  total_goals = matrix(NA, nrow = 50, ncol = 50)
  
  for (i in 1:50){
    for (j in 1:50){
      total_goals[i, j] = (i-1) + (j-1)
    }
  }
  
  over_prob = sum(all_possibilities[total_goals > over_under])
  under_prob = sum(all_possibilities[total_goals < over_under])
  
  print(over_prob + under_prob)
  
  over_odds = 1/over_prob
  under_odds = 1/under_prob
  
  return(c(over_odds, under_odds))
  
}

calculate_ou_75 = function(expg1, expg2, over_under){
  
  stopifnot(over_under%%1 == 0.75,
            over_under > 0)
  
  all_possibilities = as.matrix(dpois(0:49, expg1) %o% dpois(0:49, expg2))
  total_goals = matrix(NA, nrow = 50, ncol = 50)
  
  for (i in 1:50){
    for (j in 1:50){
      total_goals[i, j] = (i-1) + (j-1)
    }
  }
  
  under_odds_05 = calculate_ou_05(expg1, expg2, over_under - 0.25)[2]
  under_odds_00 = calculate_ou_00(expg1, expg2, over_under + 0.25)[2]
  
  under_odds = (under_odds_05 + under_odds_00)/2
  under_prob = 1/under_odds
  over_prob = 1 - under_prob
  over_odds = 1/over_prob
  
  return(c(over_odds, under_odds))
  
}

## Handicap

calculate_ah_00 = function(expg1, expg2, home_handicap){
  
  stopifnot(home_handicap%%1 == 0)
  
  all_possibilities = as.matrix(dpois(0:49, expg1) %o% dpois(0:49, expg2))
  goal_diff = matrix(NA, nrow = 50, ncol = 50)
  
  for (i in 1:50){
    for (j in 1:50){
      goal_diff[i, j] = (i-1) - (j-1)
    }
  }
  
  min_goal_diff = 1 - home_handicap # minimum goal difference if the home team wants to win
  draw_goal_diff = - home_handicap # goal difference for draws
  
  home_win_prob = sum(all_possibilities[goal_diff >= min_goal_diff])
  draw_prob = sum(all_possibilities[goal_diff == draw_goal_diff])
  away_win_prob = sum(all_possibilities[goal_diff < draw_goal_diff])
  
  print(sum(home_win_prob + draw_prob + away_win_prob))
  
  new_home_prob = home_win_prob / (home_win_prob + away_win_prob)
  new_away_prob = away_win_prob / (home_win_prob + away_win_prob)
  
  h_odds = 1/new_home_prob
  a_odds = 1/new_away_prob
  
  return(c(h_odds, a_odds))
  
}

calculate_ah_05 = function(expg1, expg2, home_handicap){
  
  stopifnot(home_handicap%%1 == 0.5)
  
  all_possibilities = as.matrix(dpois(0:49, expg1) %o% dpois(0:49, expg2))
  goal_diff = matrix(NA, nrow = 50, ncol = 50)
  
  for (i in 1:50){
    for (j in 1:50){
      goal_diff[i, j] = (i-1) - (j-1)
    }
  }
  
  min_goal_diff = 0.5 - home_handicap # minimum goal difference if the home team wants to win
  home_win_prob = sum(all_possibilities[goal_diff >= min_goal_diff])
  
  h_odds = 1/home_win_prob
  a_odds = 1/(1 - home_win_prob)
  
  return(c(h_odds, a_odds))
  
}

calculate_ah_25 = function(expg1, expg2, home_handicap){
  
  stopifnot(abs(home_handicap)%%1 == 0.25)
  all_possibilities = as.matrix(dpois(0:49, expg1) %o% dpois(0:49, expg2))
  goal_diff = matrix(NA, nrow = 50, ncol = 50)
  
  for (i in 1:50){
    for (j in 1:50){
      goal_diff[i, j] = (i-1) - (j-1)
    }
  }
  
  if (home_handicap > 0){
    
    home_handicap_00 = home_handicap - 0.25
    h_00_odds = calculate_ah_00(expg1, expg2, home_handicap_00)[1]
    
    home_handicap_05 = home_handicap + 0.25
    h_05_odds = calculate_ah_05(expg1, expg2, home_handicap_05)[1]
    
  }
  
  if (home_handicap < 0){
    
    home_handicap_00 = home_handicap + 0.25
    h_00_odds = calculate_ah_00(expg1, expg2, home_handicap_00)[1]
    
    home_handicap_05 = home_handicap - 0.25
    h_05_odds = calculate_ah_05(expg1, expg2, home_handicap_05)[1]
    
  }
  
  h_25_odds = (h_00_odds + h_05_odds) / 2
  h_25_prob = 1/h_25_odds
  a_25_prob = 1 - h_25_prob
  a_25_odds = 1/a_25_prob
  
  return(c(h_25_odds, a_25_odds))
  
}

calculate_ah_75 = function(expg1, expg2, home_handicap){
  
  stopifnot(abs(home_handicap)%%1 == 0.75)
  all_possibilities = as.matrix(dpois(0:49, expg1) %o% dpois(0:49, expg2))
  goal_diff = matrix(NA, nrow = 50, ncol = 50)
  
  for (i in 1:50){
    for (j in 1:50){
      goal_diff[i, j] = (i-1) - (j-1)
    }
  }
  
  if (home_handicap > 0){
    
    home_handicap_00 = home_handicap + 0.25
    a_00_odds = calculate_ah_00(expg1, expg2, home_handicap_00)[2]
    
    home_handicap_05 = home_handicap - 0.25
    a_05_odds = calculate_ah_05(expg1, expg2, home_handicap_05)[2]
    
  }
  
  if (home_handicap < 0){
    
    home_handicap_00 = home_handicap - 0.25
    a_00_odds = calculate_ah_00(expg1, expg2, home_handicap_00)[2]
    
    home_handicap_05 = home_handicap + 0.25
    a_05_odds = calculate_ah_05(expg1, expg2, home_handicap_05)[2]
    
  }
  
  a_25_odds = (a_00_odds + a_05_odds) / 2
  a_25_prob = 1/a_25_odds
  h_25_prob = 1 - a_25_prob
  h_25_odds = 1/h_25_prob
  
  return(c(h_25_odds, a_25_odds))
  
}
