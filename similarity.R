pmi <- function(q, t) {
  
  q <- q[!is.na(q)]
  t <- t[!is.na(t)]
  
  q_l <- length(q)
  t_l <- length(t)
  
  alphabet <- c(letters, LETTERS, 0:9, strsplit("!@#$%^&*()[]{}", "")[[1]])
  
  u <- as.integer(factor(c(q, t)))
  
  if(length(unique(u)) > length(alphabet)){
    warning("Too many levels for pairwiseAlignment")
    return(NA)
  }
  
  q_enc <- alphabet[u[1:q_l]] %>% paste(collapse = "")
  t_enc <- alphabet[u[(q_l+1):length(u)]] %>% paste(collapse = "")
  
  aligned <- pwalign::pairwiseAlignment(q_enc,
                                        t_enc,
                                        type = "global", # i.e., Needleman-Wunsch
                                        gapOpening = 12,
                                        gapExtension = 6)
  
  q_aligned <- aligned@pattern %>% as.character() %>% str_split("") %>% unlist()
  t_aligned <- aligned@subject %>% as.character() %>% str_split("") %>% unlist()
  browser()
  sum(q_aligned == t_aligned) / ((q_l + t_l)/2)
}
edit_sim_utf8 <- function(s, t) {
  if(!is.numeric(s) || !is.numeric(t)){
    logging::logerror("Called edit_sim_utf8 with non-numeric vectors")
    stop()
  }
  s_l <- length(s)
  t_l <- length(t)
  u <- as.integer(factor(c(s, t)))
  alphabet <- c(letters, LETTERS, 0:9, strsplit("!@#$%^&*()[]{}", "")[[1]])
  
  if(length(unique(u)) > length(alphabet)){
    warning("Too many levels for pairwiseAlignment")
    return(NA)
  }
  
  s_enc <- alphabet[u[1:s_l]] %>% paste(collapse = "")
  t_enc <- alphabet[u[(t_l+1):length(u)]] %>% paste(collapse = "")
  browser()
  s <- s[!is.na(s)]
  t <- t[!is.na(t)]
  
  offset <- min(c(s, t))
  s <- s -  offset + 256
  t <- t -  offset  + 256
  1 - utils::adist(intToUtf8(s),intToUtf8(t))[1,1]/max(length(s), length(t))
}
