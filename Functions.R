# Functions

odds_to_prob = function(odds, method='basic', normalize=TRUE,
                        target_probability = 1,
                        grossmargin = 0, shin_method = 'js', 
                        shin_maxiter = 1000,
                        uniroot_options = NULL){
  
  stopifnot(length(method) == 1,
            tolower(method) %in% c('basic', 'shin', 'bb', 'wpo', 'or', 'power', 'additive', 'jsd'),
            all(odds >= 1, na.rm=TRUE),
            length(target_probability) == 1,
            target_probability > 0,
            grossmargin >= 0,
            shin_method %in% c('js', 'uniroot'),
            length(shin_method) == 1,
            length(shin_maxiter) == 1,
            shin_maxiter > 1,
            is.null(uniroot_options) | is.list(uniroot_options))
  
  
  if (!is.matrix(odds)){
    
    if ('data.frame' %in% class(odds)){
      odds = as.matrix(odds)
    } else {
      odds = matrix(odds, nrow=1,
                    dimnames = list(NULL, names(odds)))
    }
  }
  
  
  # Prepare the list that will be returned.
  out = vector(mode='list', length=2)
  names(out) = c('probabilities', 'margin')
  
  # Some useful quantities
  n_odds = nrow(odds)
  n_outcomes = ncol(odds)
  
  # Inverted odds and margins
  inverted_odds = 1 / odds
  inverted_odds_sum = rowSums(inverted_odds)
  out$margin = inverted_odds_sum - target_probability
  
  # Missing values
  missing_idx = apply(odds, MARGIN = 1,
                      FUN = function(x) any(is.na(x)))
  
  
  if (any(inverted_odds_sum[!missing_idx] < 1)){
    stop('Some inverse odds sum to less than 1.')
  }
  
  if (method == 'basic'){
    out$probabilities = (target_probability * inverted_odds) / inverted_odds_sum
    
  } 
  
  ## do a final normalization to make sure the probabilities
  ## sum to 1 without rounding errors.
  if (normalize){
    out$probabilities = (target_probability * out$probabilities) / rowSums(out$probabilities)
  }
  
  # Make sure the matrix of implied probabilities has column names.
  if (!is.null(colnames(odds))){
    colnames(out$probabilities) = colnames(odds)
  }
  
  # check if there are any probabilities outside the 0-1 range.
  problematic = apply(out$probabilities, MARGIN = 1, FUN=function(x){any(x > 1 | x < 0)})
  problematic[is.na(problematic)] = TRUE
  problematic[missing_idx] = NA
  
  if (any(problematic, na.rm=TRUE)){
    warning(sprintf('Probabilities outside the 0-1 range produced at %d instances.\n',
                    sum(problematic)))
  }
  
  out$problematic = problematic
  
  return(out)
}

expg_prob_sq_error <- function(pars, trgt_probs, rho, uprx){
  
  pars <- exp(pars) # trick to avoid negaive lambda parameters.
  hda_probs <- numeric(3)
  probmat <- stats::dpois(0:uprx, lambda=pars[1]) %o% stats::dpois(0:uprx, lambda=pars[2])
  
  ## DC adjustment
  if (rho != 0){
    correctionmat <- matrix(tau(c(0,1,0,1), c(0,0,1,1),
                                rep(pars[1], 4),
                                rep(pars[2], 4), rho), nrow=2)
    probmat[1:2, 1:2] <- probmat[1:2, 1:2] * correctionmat
  }
  
  hda_probs[2] <- sum(diag(probmat))
  hda_probs[1] <- sum(probmat[lower.tri(probmat)])
  hda_probs[3] <- 1 - sum(hda_probs[1:2])
  
  sum((hda_probs - trgt_probs)^2)
  
}

prob_to_expg = function(probabilities, rho=0, uprx=75){
  
  # Convert to matrix
  if (!is.matrix(probabilities)){
    if (is.numeric(probabilities)){
      probabilities <- matrix(probabilities, nrow=1,
                              dimnames = list(NULL, names(probabilities)))
    } else {
      probabilities <- as.matrix(probabilities)
    }
    
    
  }
  
  stopifnot(ncol(probabilities) == 3,
            all(abs(rowSums(probabilities) - 1) < 0.0001),
            uprx >= 1,
            length(uprx) == 1,
            length(rho) == 1)
  
  expg <- matrix(ncol=2, nrow=nrow(probabilities))
  sq_errors <- numeric(nrow(probabilities))
  
  for (ii in 1:nrow(probabilities)){
    
    optim_res <- stats::optim(c(0,0), fn=expg_prob_sq_error,
                              trgt_prob=probabilities[ii,],
                              rho = rho, uprx = uprx)
    
    expg[ii,] <- exp(optim_res$par)
    sq_errors[ii] <- optim_res$value
    
  }
  
  out <- list(expg = expg, sq_errors=sq_errors)
  
  return(out)
}

# The negative log-likelihood function 
negloglik = function(params, goals1, goals2, team1, team2,
                     x1, x2, hfa, model, param_skeleton, all_teams,
                     weights = NULL, fixed_params = NULL){
  
  # relist, to make things easier.
  plist = utils::relist(params, param_skeleton)
  
  # TODO: There is a bug here, i think, if one of the teams with
  # a fixed attack or defense parameter is the first team
  # in the all_teams vector.
  
  # Add sum to zero constraint on defense parameters.
  # The defense parameter for the first team is computed from the rest.
  if (length(plist$defense) != 0){
    # if 0, then all defense parameters are presumably fixed.
    plist$defense = c(sum(plist$defense)*-1, plist$defense)
    names(plist$defense)[1] = all_teams[1] # add name to first element.
  }
  
  if (length(plist$attack) != 0){
    # if 0, then all attack parameters are presumably fixed.
    plist$attack = c(sum(plist$attack)*-1, plist$attack)
    names(plist$attack)[1] = all_teams[1] # add name to first element.
  }
  
  # Add the fixed parameters to the parameter list.
  plist = fill_fixed_params(plist, fixed_params = fixed_params)
  
  ## Expected goals (poisson & nbin parameters)
  expg = lambda_pred(plist, team1, team2, x1, x2)
  
  # The log-likelihood
  if (model == 'poisson'){
    log_lik_1 = stats::dpois(goals1, lambda = expg$expg1, log=TRUE)
    log_lik_2 = stats::dpois(goals2, lambda = expg$expg2, log=TRUE)
  } else if (model == 'negbin'){
    dispersion_tmp = 1 / exp(plist$dispersion)
    log_lik_1 = stats::dnbinom(goals1, mu = expg$expg1, size = dispersion_tmp, log=TRUE)
    log_lik_2 = stats::dnbinom(goals2, mu = expg$expg2, size = dispersion_tmp, log=TRUE)
  } else if (model == 'gaussian'){
    log_lik_1 = stats::dnorm(goals1, mean = expg$expg1, sd = exp(plist$sigma), log=TRUE)
    log_lik_2 = stats::dnorm(goals2, mean = expg$expg2, sd = exp(plist$sigma), log=TRUE)
  } else if (model == 'cmp'){
    exp_log_upsilon = exp(plist$dispersion)
    
    # This is a dirty hack essentially setting hard upper and lower
    # bounds for the the dispersion parameter.
    if (exp_log_upsilon < 0.7){return(Inf)}
    if (exp_log_upsilon > 1.7){return(Inf)}
    
    log_lik_1 = dCMP(goals1, lambda = lambdaCMP(mu=expg$expg1,
                                                upsilon = exp_log_upsilon,
                                                method = 'fast'),
                     upsilon = exp_log_upsilon, log=TRUE)
    log_lik_2 = dCMP(goals2, lambda = lambdaCMP(mu=expg$expg2,
                                                upsilon = exp_log_upsilon,
                                                method = 'fast'),
                     upsilon = exp_log_upsilon, log=TRUE)
  }
  
  log_lik_terms = log_lik_1 + log_lik_2
  
  # Dixon-Coles adjustment
  if ('rho' %in% names(plist)){
    
    dc_adj = tau(goals1, goals2, expg$expg1, expg$expg2, rho = plist$rho)
    
    # Trick to avoid warnings.
    if (any(dc_adj <= 0.0)){
      return(Inf)
    }
    
    log_lik_terms = log_lik_terms + log(dc_adj)
  }
  
  
  if ('hurdle' %in% names(plist)){
    
    zz_idx = goals1 == 0 & goals2 == 0
    
    hurdle_adj = exp(plist$hurdle) * exp(-(expg$expg1 + expg$expg2))
    
    if (any(hurdle_adj >= 1)){
      return(Inf)
    }
    
    log_lik_terms[zz_idx] = log(hurdle_adj[zz_idx])
    log_lik_terms[!zz_idx] = log(1 - hurdle_adj[!zz_idx]) + log_lik_terms[!zz_idx]
    
  }
  
  # sum the log likelihood.
  if (!is.null(weights)){
    log_lik = sum(log_lik_terms*weights)
  } else {
    log_lik = sum(log_lik_terms)
  }
  
  return(log_lik*-1)
  
}

# Check that a list of parameters has correct layout. Returns TRUE if everything is OK.
check_plist = function(plist, all_teams = NULL){
  
  t1 = is.list(plist)
  
  plist_names = names(plist)
  
  t2 = all(plist_names %in% c('attack', 'defense', 'beta', 'intercept',
                              'sigma', 'rho', 'dispersion', 'gamma', 'hfa'))
  
  t3 = all(sapply(plist, is.numeric))
  
  # check the length of the vectors in the list.
  t4 = TRUE
  for (ii in 1:length(plist_names)){
    if (plist_names[ii] %in% c('intercept', 'sigma', 'dispersion', 'rho', 'gamma', 'hfa')){
      t4 = t4 & length(plist[[plist_names[ii]]] == 1)
    }
  }
  
  # (Optional) Check that the attack and defence parameters
  # are available for all teams.
  t5 = TRUE
  if (!is.null(all_teams)){
    if ('attack' %in% plist_names){
      t5 = t5 & all(names(plist$attack) %in% all_teams)
    }
    
    if ('defense' %in% plist_names){
      t5 = t5 & all(names(plist$defense) %in% all_teams)
    }
    
  }
  
  
  out = all(t1, t2, t3, t4)
  return(out)
}

is_connected = function(edgelist){
  
  all_nodes = sort(unique(c(edgelist[,1], edgelist[,2]))) # all nodes
  n_nodes = length(all_nodes) # number of nodes
  
  start_node = all_nodes[1]
  
  # initiate the lgocial vector for nodes that has been reachef from the
  # starting node.
  reachable_from_start = all_nodes == start_node
  
  # logical of length n_nodes, indicating nodes that has been visited.
  # These should not be added to queue, hence the name dequed.
  dequed_nodes = all_nodes == start_node
  
  # logical of length n_teams, indicating queue.
  node_queue = all_nodes == start_node
  
  while (sum(node_queue) != 0){
    
    # First node in the queue.
    cur_node = all_nodes[node_queue][1]
    
    reachable_from_current = union(edgelist[edgelist[,1] == cur_node,2], edgelist[edgelist[,2] == cur_node,1])
    reachable_from_current_idx = all_nodes %in% reachable_from_current & (!dequed_nodes)
    reachable_from_start[reachable_from_current_idx] = reachable_from_start[reachable_from_current_idx] | TRUE
    
    # list of nodes that has been visited, and had their neighbors visited.
    dequed_nodes[all_nodes == cur_node] = TRUE
    
    # Update node queue
    node_queue = (node_queue | (all_nodes %in% reachable_from_current)) & (!dequed_nodes)
    
    if (sum(reachable_from_start) == n_nodes){break}
    
  }
  
  sum(reachable_from_start) == n_nodes
  
}

gm_fit_glm = function(goals1, goals2, team1, team2, all_teams, model,
                      x1 = NULL, x2 = NULL,
                      hfa = TRUE,  weights = NULL){
  
  # prepare some
  n_teams = length(all_teams)
  team1_stacked = c(team1, team2)
  team2_stacked = c(team2, team1) # Opponent
  
  # Make response (y) vector
  yy = c(goals1, goals2)
  
  # Attack dummy matrix
  xmata = matrix(0, nrow=length(goals2)*2,
                 ncol = (n_teams-1))
  colnames(xmata) = all_teams[-1]
  
  # Defense dummy matrix
  xmatd = matrix(0, nrow=length(goals2)*2,
                 ncol = (n_teams-1))
  colnames(xmatd) = all_teams[-1]
  
  for (ii in 2:n_teams){
    
    t1_idx = team1_stacked == all_teams[ii]
    t2_idx = team2_stacked == all_teams[ii]
    xmata[t1_idx, all_teams[ii]] = 1
    xmatd[t2_idx, all_teams[ii]] = -1
  }
  
  # Sum-to-zero constraint for the first team.
  xmata[team1_stacked == all_teams[1]] = -1
  xmatd[team2_stacked == all_teams[1]] = 1
  
  # Combine the attack and defence matrices.
  colnames(xmata) = paste(colnames(xmata), 'Attack', sep='_')
  colnames(xmatd) = paste(colnames(xmatd), 'Defense', sep='_')
  
  # Add home field advantage.
  if (hfa){
    xmat = cbind(intercept = 1, hfa = rep(1:0, each=length(goals2)),
                 xmata, xmatd)
  } else {
    xmat = cbind(intercept = 1, xmata, xmatd)
  }
  
  
  if (model %in% c('poisson', 'negbin', 'cmp')){
    # Note that a poisson model is fitted if model = 'negbin' or 'cmp'.
    glm_family = stats::poisson(link='log')
  } else if (model == 'gaussian'){
    glm_family = stats::gaussian(link='log')
  }
  
  if (is.null(weights)){
    glm_res = stats::glm.fit(x=xmat, y=yy,
                             start=c(mean(log(yy+0.5)-0.2), rep(0, ncol(xmat)-1)),
                             family=glm_family,
                             control = list(maxit=100, epsilon = 1e-9),
                             intercept = FALSE)
  } else {
    glm_res = stats::glm.fit(x=xmat, y=yy,
                             weights=rep(weights, 2),
                             start=c(mean(log(yy+0.5)-0.2), rep(0, ncol(xmat)-1)),
                             family=glm_family,
                             control = list(maxit=100, epsilon = 1e-9),
                             intercept = FALSE)
  }
  
  attack_params = glm_res$coefficients[grepl('_Attack$', names(glm_res$coefficients))]
  names(attack_params) = sub('_Attack$', '', names(attack_params))
  attack_params = c(sum(attack_params)*-1, attack_params)
  names(attack_params)[1] = all_teams[1]
  
  defense_params = glm_res$coefficients[grepl('_Defense$', names(glm_res$coefficients))]
  names(defense_params) = sub('_Defense$', '', names(defense_params))
  defense_params = c(sum(defense_params)*-1, defense_params)
  names(defense_params)[1] = all_teams[1]
  
  param_list = list(attack = attack_params,
                    defense = defense_params,
                    intercept = glm_res$coefficients['intercept'])
  
  if (model == 'gaussian'){
    param_list$sigma = stats::sd(glm_res$residuals)
  }
  
  names(param_list$intercept) = NULL # The intercept vector should not be named.
  
  if (hfa){
    param_list$hfa = glm_res$coefficients['hfa']
    names(param_list$hfa) = NULL # The hfa vector should not be named.
  }
  
  stopifnot(check_plist(param_list))
  
  # Compute the log-likelihood.
  if (model %in% c('poisson', 'negbin', 'cmp')){
    if (is.null(weights)){
      loglik = sum(stats::dpois(x = c(goals1, goals2),
                                lambda=glm_res$fitted.values, log=TRUE))
    } else {
      loglik = sum(stats::dpois(x = c(goals1, goals2),
                                lambda=glm_res$fitted.values, log=TRUE)*rep(weights, 2))
    }
    
  } else if (model == 'gaussian'){
    if (is.null(weights)){
      loglik = sum(stats::dnorm(x = c(goals1, goals2),
                                mean=glm_res$fitted.values, log=TRUE))
    } else {
      loglik = sum(stats::dnorm(x = c(goals1, goals2),
                                mean=glm_res$fitted.values, log=TRUE)*rep(weights,2))
    }
  }
  
  
  return(list(parameters = param_list, loglikelihood=loglik,
              npar_fixed = 0, npar_est = ncol(xmat), aic=glm_res$aic,
              converged = glm_res$converged,
              boundary = glm_res$boundary))
  
}

expg_model = function(goals1, goals2, team1, team2,
                      x1 = NULL, x2 = NULL,
                      hfa = TRUE, dc = FALSE, rs = FALSE, hurdle = FALSE,
                      fixed_params = NULL, weights = NULL,
                      model, optim_method='BFGS'){
  
  warning_messages = c()
  
  stopifnot(length(goals1) == length(goals2),
            length(goals2) == length(team1),
            length(team1) == length(team2),
            length(goals1) >= 1,
            is.numeric(goals1), is.numeric(goals2),
            model %in% c('poisson', 'negbin', 'gaussian', 'cmp'))
  
  
  if (!is_connected(cbind(team1, team2))){
    wm_tmp = 'The data contains two or more separate clusters of teams that are not comparable. The results may be unreliable.'
    warning_messages = append(warning_messages, wm_tmp)
    warning(wm_tmp)
  }
  
  
  if (!is.null(weights)){
    stopifnot(is.numeric(weights),
              length(goals1)==length(weights),
              all(weights >= 0),
              all(!is.na(weights)),
              !all(weights == 0))
  }
  
  
  # Make sure the team vectors are of the charcter type.
  team1 = as.character(team1)
  team2 = as.character(team2)
  
  # Some useful quantities.
  all_teams = sort(unique(c(unique(team1), unique(team2))), decreasing = FALSE)
  n_teams = length(all_teams)
  ngames = length(goals1)
  
  # If it is sufficient to fit the model with glm.fit().
  mdefault = model %in% c('poisson', 'gaussian') & dc == FALSE & rs == FALSE & hurdle == FALSE & is.null(fixed_params)
  
  fitter = ifelse(mdefault, 'glm.fit', 'gm')
  
  # Fit a model with glm.fit. Some model classes are compatible with this function.
  # If not, the results from fitting this simples model is used as starting values.
  
  # TODO: If all attack and defence and hfa and intercept are fixed, this is not
  # needed for fitting the rest of the models. This step should be skiped
  # to save computing time.
  
  # TODO: Maybe even negbin models can be fitted with this approach, with
  # the glm.nb() function from the MASS pacakge.
  
  start_time = Sys.time()
  
  gm_fit_glm_res = gm_fit_glm(goals1=goals1, goals2=goals2,
                              team1=team1, team2=team2,
                              all_teams = all_teams,
                              model = model,
                              weights = weights,
                              x1 = x1, x2=x2, hfa=hfa)
  
  
  if (mdefault){
    
    if (!gm_fit_glm_res$converged){
      wm_tmp = 'Did not converge (glm.fit). Parameter estimates are unreliable.'
      warning_messages = append(warning_messages, wm_tmp)
      warning(wm_tmp)
    }
    
    if (gm_fit_glm_res$boundary){
      wm_tmp = 'glm.fit(): Fitted values on the boundary of the attainable values. Parameter estimates are unreliable.'
      warning_messages = append(warning_messages, wm_tmp)
      warning(wm_tmp)
    }
    
    parameter_list = gm_fit_glm_res$parameters
    loglikelihood = gm_fit_glm_res$loglikelihood
    npar_est = gm_fit_glm_res$npar_est
    npar_fixed = gm_fit_glm_res$npar_fixed
    aic = gm_fit_glm_res$aic
    converged = gm_fit_glm_res$converged
    
  } else {
    
    # initial values from gm_fit_glm_res.
    parameter_list_init = gm_fit_glm_res$parameters
    parameter_list_init$attack = parameter_list_init$attack[-1]
    parameter_list_init$defense = parameter_list_init$defense[-1]
    
    # Add additional initial values.
    if (dc){
      parameter_list_init$rho = 0.01
    }
    
    if (rs){
      parameter_list_init$gamma = 0.0
    }
    
    if (hurdle){
      parameter_list_init$hurdle = 0.0
    }
    
    if (model == 'negbin'){
      # on log scale during estimation.
      # start with a value close to 0, which is almost Poisson.
      parameter_list_init$dispersion = -10
    }
    
    if (model == 'gaussian'){
      # on log scale during estimation.
      parameter_list_init$sigma = log(parameter_list_init$sigma)
    }
    
    if (model == 'cmp'){
      # on log scale during estimation.
      # start with a value exp(0)=1, which is the same as Poisson.
      parameter_list_init$dispersion = 0
    }
    
    # Deal with fixed parameters.
    if (!is.null(fixed_params)){
      
      stopifnot(is.list(fixed_params))
      
      if (any(!names(fixed_params) %in% c('attack', 'defense', 'beta', 'intercept',
                                          'sigma', 'rho', 'dispersion', 'gamma', 'hfa'))){
        stop('In fixed_params: Invalid parameter name.')
      }
      
      # remove fixed parameters from the parameter_list_init, since they are
      # not optimized over.
      if ('attack' %in% names(fixed_params)){
        fixed_attack_params = names(parameter_list_init$attack) %in% names(fixed_params$attack)
        parameter_list_init$attack = parameter_list_init$attack[!fixed_attack_params]
      }
      
      if ('defense' %in% names(fixed_params)){
        fixed_defence_params = names(parameter_list_init$defense) %in% names(fixed_params$defense)
        parameter_list_init$defense = parameter_list_init$defense[!fixed_defence_params]
      }
      
      # remove the parameters that are not attack aand defence parameters.
      if (any(!names(fixed_params) %in% c('attack', 'defense'))){
        
        parameter_list_init[names(fixed_params)] = NULL
        
        if ('dispersion' %in% names(fixed_params)){
          # During estimation, the dispersion parameter is on the log scale
          # to avoid negative values.
          fixed_params$dispersion = log(fixed_params$dispersion)
          
          if (model == 'poisson'){
            wm_tmp = 'Dispersion parameter is fixed, but model is Poisson. The dispersion parameter will not have an effect.'
            warning_messages = append(warning_messages, wm_tmp)
            warning(wm_tmp)
          } else if (model == 'gaussian'){
            wm_tmp = 'Dispersion parameter is fixed, but model is Gaussian The dispersion parameter will not have an effect. The related parameter for the Gaussian model is sigma.'
            warning_messages = append(warning_messages, wm_tmp)
            warning(wm_tmp)
          }
        }
      }
      
    } # end fixed params
    
    # Commence estimation.
    parameter_vector = unlist(parameter_list_init)
    
    optim_res = stats::optim(par = parameter_vector, fn=negloglik,
                             goals1 = goals1, goals2 = goals2,
                             team1=team1, team2=team2,
                             x1 = x1, x2 = x2,
                             fixed_params=fixed_params,
                             model = model,
                             all_teams = all_teams,
                             param_skeleton=parameter_list_init,
                             weights = weights,
                             method = optim_method,
                             control = list(maxit = 250))
    
    converged = optim_res$convergence == 0
    
    if (!converged){
      wm_tmp = 'Did not converge (optim). Parameter estimates are unreliable.'
      warning_messages = append(warning_messages, wm_tmp)
      warning(wm_tmp)
    }
    
    # relist the parameter vector, calculate the missing attack and defense parameter.
    parameter_list_est = utils::relist(optim_res$par, parameter_list_init)
    
    if (length(parameter_list_est$defense) != 0){
      # TODO: There is probably a bug here, if there is only one parameter that
      # is not fixed. Also in the negloglik() function.
      
      # if 0, then all defense parameters are presumably fixed.
      parameter_list_est$defense = c(sum(parameter_list_est$defense)*-1, parameter_list_est$defense)
      names(parameter_list_est$defense)[1] = all_teams[1]
    }
    
    if (length(parameter_list_est$attack) != 0){
      # if 0, then all attack parameters are presumably fixed.
      parameter_list_est$attack = c(sum(parameter_list_est$attack)*-1, parameter_list_est$attack)
      names(parameter_list_est$attack)[1] = all_teams[1]
    }
    
    loglikelihood = optim_res$value*-1
    npar_est = length(optim_res$par)
    npar_fixed = length(unlist(fixed_params))
    aic = npar_est*2 - 2*loglikelihood
    
    # rescale dispersion
    if ('dispersion' %in% names(parameter_list_est)){
      parameter_list_est$dispersion = exp(parameter_list_est$dispersion)
    }
    
    # rescale sigma
    if ('sigma' %in% names(parameter_list_est)){
      parameter_list_est$sigma = exp(parameter_list_est$sigma)
    }
    
    # rescale hurlde parameter
    if ('hurdle' %in% names(parameter_list_est)){
      parameter_list_est$hudrle = exp(parameter_list_est$hurdle)
    }
    
    # Add the fixed parameters to the parameter list.
    parameter_list = fill_fixed_params(parameter_list_est, fixed_params = fixed_params)
    
  }
  
  end_time = Sys.time()
  est_time = difftime(end_time, start_time, units='secs')
  
  # Compute R squared.
  all_goals = c(goals1, goals2)
  if (is.null(weights)){
    mean_goals = mean(all_goals)
  } else {
    mean_goals = stats::weighted.mean(all_goals, w = rep(weights, 2))
  }
  
  
  ## Deviances needed for R squared.
  if (model == 'poisson'){
    if (is.null(weights)){
      loglikelihood_saturated = sum(stats::dpois(all_goals, lambda = all_goals, log=TRUE))
      loglikelihood_null = sum(stats::dpois(all_goals, lambda = mean_goals, log=TRUE))
    } else {
      loglikelihood_saturated = sum(stats::dpois(all_goals, lambda = all_goals, log=TRUE)*rep(weights,2))
      loglikelihood_null = sum(stats::dpois(all_goals, lambda = mean_goals, log=TRUE)*rep(weights,2))
    }
    
  }  else if (model == 'gaussian'){
    if (is.null(weights)){
      sigma0_tmp = stats::sd(all_goals)
      loglikelihood_saturated = sum(stats::dnorm(all_goals, mean = all_goals, sd=sigma0_tmp, log=TRUE))
      loglikelihood_null = sum(stats::dnorm(all_goals, mean = mean_goals, sd=sigma0_tmp, log=TRUE))
    } else {
      sigma0_tmp = sqrt(sum(rep(weights, 2) * (all_goals - mean_goals)^2))
      loglikelihood_saturated = sum(stats::dnorm(all_goals, mean = all_goals, sd=sigma0_tmp, log=TRUE)*rep(weights,2))
      loglikelihood_null = sum(stats::dnorm(all_goals, mean = mean_goals, sd=sigma0_tmp, log=TRUE)*rep(weights,2))
    }
  } 
  
  # TODO: How should deviances and R^2 be computed in the Dixon-Coles model?
  
  deviance = 2 * (loglikelihood_saturated - loglikelihood)
  deviance_null = 2 * (loglikelihood_saturated - loglikelihood_null)
  
  r_squared = 1 - (deviance / deviance_null)
  
  
  # Update the all_teams variable to include teams not in data, but with
  # fixed parameters.
  all_param_teams = unique(c(names(parameter_list$defense), names(parameter_list$attack)))
  if (any(!all_param_teams %in% all_teams)){
    all_teams = sort(all_param_teams)
  }
  
  # sort the attack and defence parameters alphabetically
  parameter_list$defense = parameter_list$defense[order(names(parameter_list$defense))]
  parameter_list$attack = parameter_list$attack[order(names(parameter_list$attack))]
  
  
  if (any(is.na(unlist(parameter_list)))){
    wm_tmp = 'Some parameters are NA.'
    warning_messages = append(warning_messages, wm_tmp)
    warning(wm_tmp)
  }
  
  # maxgoal. Useful for later predictions.
  maxgoal = max(all_goals)
  
  out = list(parameters = parameter_list,
             loglikelihood = loglikelihood, npar_est=npar_est,
             npar_fixed = npar_fixed, aic=aic, r_squared=r_squared,
             all_teams = all_teams, ngames = ngames,
             est_time = est_time, model = model,
             fixed_params = fixed_params,
             converged = converged,
             maxgoal = maxgoal,
             fitter = fitter,
             warnings = warning_messages)
  
  class(out) = 'goalmodel'
  
  return(out)
  
}

lambda_pred = function(plist, team1, team2, x1, x2){
  
  # Team 1 log-expected goals
  eta1 = plist$intercept + plist$attack[team1] - plist$defense[team2]
  
  if ('hfa' %in% names(plist)){
    eta1 = eta1 + plist$hfa
  }
  
  # Team 2 log- expected goals
  eta2 = plist$intercept + plist$attack[team2] - plist$defense[team1]
  
  # Psychological underestimation factor (Rue & Salvesen 1997).
  if ('gamma' %in% names(plist)){
    deltas = delta_ab(plist$attack[team1], plist$defense[team1],
                      plist$attack[team2], plist$defense[team2])
    gamma_delta = plist$gamma * deltas
    eta1 = eta1 - gamma_delta
    eta2 = eta2 + gamma_delta
  }
  
  
  # link function.
  lambda1 = exp(eta1)
  lambda2 = exp(eta2)
  
  out = list(expg1 = lambda1, expg2 = lambda2)
  return(out)
  
}

pred_expg = function(model_fit, team1, team2, 
                     x1 = NULL, x2 = NULL, return_df = FALSE){
  
  stopifnot(length(team1) == length(team2))
  
  team1 = as.character(team1)
  team2 = as.character(team2)
  
  ee = lambda_pred(model_fit$parameters, team1, team2, x1, x2)
  
  if (any(is.na(ee$expg1) | is.na(ee$expg2))){
    warning('Could not make predictions in some instances.')
  }
  
  if (return_df){
    out = data.frame(team1 = team1, team2 = team2,
                     expg1 = ee$expg1, expg2 = ee$expg2,
                     stringsAsFactors = FALSE, row.names = NULL)
  } else {
    out = ee
  }
  
  
  return(out)
}

weights_dc <- function(dates, xi=0, currentDate=NULL){
  stopifnot(xi >= 0,
            length(xi) == 1)
  dates <- as.Date(dates)
  
  if (is.null(currentDate)){
    currentDate <- max(dates)
  } else {
    currentDate <- as.Date(currentDate)
  }
  
  datediffs <- dates - as.Date(currentDate)
  datediffs <- as.numeric(datediffs *-1)
  w <- exp(-1*xi*datediffs)
  w[datediffs < 0] <- 0 #Future dates should have zero weights
  return(w)
}

matchday_predict = function(data_df, played_matches, upcoming_match_days){
  
  implied_prob = odds_to_prob(data_df[, c('AvgCH', 'AvgCD', 'AvgCA')], 
                                       method = 'basic')
  implied_expg = prob_to_expg(implied_prob$probabilities)
  
  data_df$xhg = implied_expg$expg[, 1]
  data_df$xag = implied_expg$expg[, 2]
  
  data_df$expg1 = rep(NA, nrow(data_df))
  data_df$expg2 = rep(NA, nrow(data_df))
  
  new_train_set = data_df[1:played_matches, ]
  
  for (i in 1:upcoming_match_days){
    
    my_weights = weights_dc(as.Date(new_train_set$Date, format = '%d/%m/%Y'), xi=0.009*i)
    
    model = expg_model(goals1 = round(new_train_set$xhg), 
                       goals2 = round(new_train_set$xag),
                       team1 = new_train_set$HomeTeam, 
                       team2 = new_train_set$AwayTeam,
                       dc = F, model = 'gaussian',
                       weights = NULL)
    
    lower_bound = played_matches + 1 + (i-1)*10
    upper_bound = played_matches + i*10
    
    new_match_day = data_df[lower_bound:upper_bound, ]
    
    new_pred = pred_expg(model, return_df = TRUE, 
                         team1 = new_match_day$HomeTeam, 
                         team2 = new_match_day$AwayTeam)
    
    data_df$expg1[lower_bound:upper_bound] = new_pred$expg1
    data_df$expg2[lower_bound:upper_bound] = new_pred$expg2
    
    new_train_set = rbind(new_train_set, new_match_day)
    
  }
  
  return(data_df)
  
}