library(tidyverse)

messagef <- function(...) message(sprintf(...))

example <- c('nm(pul)', 'dhq', 'sq', 'dhq', 'sq', 'dhq', 'sq', 'nm(pul)')

zipfify <- function(data, target, with_plot = T){
  ret  <- data %>% 
    count(!!sym(target), sort = T) %>% 
    mutate(r = 1:nrow(.), 
           f = n/sum(n), 
           log_r = log(r), 
           log_n = log(n), 
           log_f = log(f))
  if(with_plot){
    q <- ret %>% 
      ggplot(aes(x = log_r, y = log_n)) + geom_point() + geom_smooth(method = "lm", color = "blue")+ geom_smooth(color = "red")
    print(q)
  }
  ret
}

get_bigrams <- function(data, target, sep = ""){
  data %>% mutate(lead1 = lead(!!sym(target)), bigram = sprintf("%s%s%s", !!sym(target), sep, lead1)) 
}

get_bernoulli_sequence <- function(n = 100, size = 26){
  if(size <= 26){
    tibble(roll = letters[rbinom(size = size, n = n, p = .5)])
  }
  else{
    tibble(roll = sprintf("%05d", rbinom(size = size, n = n, p = .5)))
  }
}

make_spc <- function(data, target){
  tmp <- data %>% count(!!sym(target)) %>% count(n)
  browser()
  print(table(diff(tmp$n)))
  spc(m = tmp$n, Vm = tmp$nn, V = nrow(data %>% count(!!sym(target))), N = nrow(data))
}

get_packages <- function(x, N = 2, sep = ";"){
  if(N == 1) return(x)
   x <- as.character(x)
   l <- length(x)
   if(N > l){
     stop("N too large")
   }
   if(N == l){
     return(paste(x, collapse = sep))
   }
   n_complete <- floor(l/N)
   n_rest <- l %% N
   start_pos <- c((0:(n_complete-1)) * N + 1, l)
   segments <- character(length(start_pos))
   for(i in 1:(length(start_pos) - 1)){
     segments[i] <- x[start_pos[i]:(start_pos[i+1] - 1)] %>% paste(collapse = sep)
   }
   #browser()
   if(n_rest > 0){
     segments[length(start_pos)] <- x[(l - n_rest):l] %>% paste(collapse = sep)
   }
   segments
}

remove_repetitions <- function(x){
  r <- rle(x)
  #browser()
  start_pos <- cumsum(c(1, r$lengths[1:length(r$lengths)-1]))
  x[start_pos]  
}

check_for_power_law <- function(data, target = "segment", n_sims = 100, type = c("displ", "disexp"), with_plot = F){
  m_m <- data  %>%  count(!!sym(target), sort = T) %>% pull(n) 
  r <- tibble(n = m_m, r = 1:length(m_m)) %>% 
    mutate(log_r = log(r), log_n = log(n)) %>% 
    select(log_r, log_n) %>% 
    correlation::correlation()
  
  #browser()
  map_dfr(type, function(tp){
    if(tp == "displ"){
      m_m <- m_m %>% displ$new(.)
      m_m$setXmin(estimate_xmin(m_m))
    }
    else if(tp == "dislnorm"){
      stop("'dislnorm' will throw an error")
      m_m <- m_m %>% dislnorm$new(.)
      
    }
    else if(tp == "disexp"){
      m_m <- m_m %>% disexp$new(.)
    }
    #browser()
    bs_p <- bootstrap_p(m_m, no_of_sims = n_sims, threads = parallel::detectCores())
    if(with_plot){
      plot(m_m)
      invisible(readline(prompt="Press [enter] to continue"))
      q <- GGally::ggpairs(bs_p$bootstraps %>% as.data.frame() %>% as_tibble() %>% select(2:3))
      print(q)
    }
    tibble(p = bs_p$p, type = tp, par = estimate_pars(m_m)$pars, xmin = estimate_xmin(m_m)$xmin)
  }) %>% 
    bind_rows(r %>% select(r, p) %>%  mutate(par = r*r, type = "log_log_reg") 
              %>% select(-r)) %>% 
    mutate(p = round(p, 4))
   
}

get_segments <- function(x, thresh = .5, sep =",", target = "x"){
  #browser()
  if(is.data.frame(x)){
    bigrams <- get_bigrams(x, target = target, sep = sep) %>% rename(x  =bigram)
  }
  else if(is.character(x)){
    bigrams <- get_bigrams(tibble(x = x), target = "x", sep = sep)
  }
  else{
    stop("x must be data frame or character vector")
  }
  ret <- bigrams %>% 
    group_by(x) %>% 
    mutate(n = n()) %>% 
    ungroup() %>% 
    mutate(f = n/sum(n), log_n = log(n), log_f = log(f)) %>% 
    mutate(lag_f = lead(f, 1), 
           f_ratio = lag_f/f, 
           boundary_tmp = as.integer(f_ratio <= thresh)%>% replace_na(0)) %>% 
    mutate(boundary = lag(boundary_tmp, 2)  %>% replace_na(0)) 
  ret$boundary[1] <- 1
  #browser()
  start_pos <- c(which(ret$boundary == 1), nrow(ret))
  segments <- rep("", length(start_pos)-1)
  
  for(i in 1:(length(start_pos) - 1)){
    segments[i] <- ret[[target]][start_pos[i]:(start_pos[i+1] - 1)] %>% paste(collapse = sep)
  }
  #browser()
  
  #ret
  tibble(segment = segments, 
         len = diff(start_pos),
         start_pos = start_pos[1:(length(start_pos)-1)]
         ) %>% select(start_pos, len, segment)
}

parse_song_string <- function(song_string){
 map(1:length(song_string), function(i){
    song_string[i] %>% str_remove("\\[")%>% str_remove("\\]") %>% str_remove_all("\'") %>% str_split(",")  %>% sapply(trimws)
   #tibble(sound = tmp, sound_id = i)
 }) 
}

setup_workspace <- function(){
  whale_s1 <- read.table("whale_data/science.adq7055_data_s1.txt", header = T, sep = "\t") %>% 
    as_tibble() %>% 
    janitor::clean_names()
  assign("whales1", whale_s1, globalenv())
  whale_s2 <- readr::read_csv("whale_data/science.adq7055_data_s2.csv") %>% 
    rename(id = ...1) %>% 
    mutate(string = parse_song_string(string)) %>% 
    unnest(string) %>% 
    rename(sound = string) %>% 
    group_by(id) %>% 
    mutate(sound_id = 1:n()) %>% 
    ungroup() %>% 
    mutate(sound = as.character(sound))
  assign("whales2", whale_s2, globalenv())
  
}

symbolic_ac <- function(x, max_lag = 10){
  max_lag <- min(length(x), max_lag)
  map_dfr(0:max_lag, function(i){
    n <- length(lag(x, i))
    #browser()
    r <- mean(na.omit(x[1:n] == lead(x, i)))
    tibble(lag = i, r = r)
  })
}

parse_sounds <- function(x){
  properties <- c("ascending" = "a", 
                  "descending"  ="d",
                  "modulated" = "m",
                  "high" = "h",
                  "low" = "l",
                  "n-shaped" = "n", 
                  "u-shaped" = "u")
  tmp <- x %>% 
    str_replace("a-l", "a=l") %>% 
    str_replace("ch-l", "ch=l") %>% 
     str_replace("pul-l", "pul=l") %>% 
    str_replace("l-pul",  "pul=l") %>% 
    str_replace("pul-s", "pul=s") %>% 
    str_replace("arg", "agr") %>% 
    str_replace("mods(pe)", "modws(pe)")
    
  #browser()
  
  tmp <- tmp %>% tolower() %>% str_remove_all("\\)") %>% str_split("-") 
  
  get_property <- function(x){
    map_dfr(1:length(x), function(k){
    prop <- ""
    y <- x[k]
    if(substr(y, 1, 3) == "mod"){
      prop <- "modulated"
      y <- substr(y, 4, nchar(y))
    }  
    else if(y == "m"){
      return(tibble( main = y, property = ""))
    }
    
    for(i in 1:length(properties)){
      if(substr(y, 1, 1) == properties[i]){
        #browser()
        prop <- names(properties)[i]
        y <- substr(y, 2, nchar(y))
        break
      }
    }
    tibble(full = x[k], main = y, property = prop)
  })
  }
  
  ret <- map_dfr(1:length(tmp), function(i){
    map_dfr(1:length(tmp[[i]]), function(j){
      #browser()
      el <- tmp[[i]][j]  %>% str_split_fixed( "\\(", 2) 
      if(el[1] == "ti"){
        if(el[2] == "a=l"){
          #browser()
          el[2] <- "l"
          el[1] <- "ati"
        }
        else{
          el[1] <- sprintf("%s%s", el[2], el[1])
          el[2] <- ""
        }
      }
      tibble(main_id = i, sub_id = j,  original = x[i], y = el[1], spec = el[2]) 
      
    })
  })
  ret <- ret %>% bind_cols(get_property(ret$y)) %>% select(main_id, sub_id, original, property, main, spec)
  browser()
  ret[ret$main == "cd",]$main <- "dc"
  ret[ret$main == "gu",]$main <- "gr"
  ret[ret$main == "sti",]$spec <- "s"
  ret[ret$main == "sti",]$main <- "ti"
  ret[ret$spec == "ch=l",]$property <- sprintf("%s-%s", ret[ret$spec == "ch=l",]$property, "l")
  ret[ret$spec == "ch=l",]$spec <- "ch"

  ret[ret$spec == "pul=l",]$property <- sprintf("%s-%s", ret[ret$spec == "pul=l",]$property, "l")
  ret[ret$spec == "pul=l",]$spec <- "pul"

  ret[ret$spec == "pul=s",]$property <- sprintf("%s-%s", ret[ret$spec == "pul=s",]$property, "s")
  ret[ret$spec == "pul=s",]$spec <- "pul"
  ret[ret$main == "g",]$main <- "gr"
  ret[ret$main == "q",]$main <- "sq"
  ret
}

check_all_songs <- function(with_segments = T){
  map_dfr(unique(whales2$id), function(sid){
    tmp <- whales2 %>% filter(id == sid) 
    messagef("Checking song %s (year: %s, cycle: %s)", sid, tmp$song_type %>% unique(), tmp$cycle %>% unique())
    #browser()
    if(with_segments){
      segments <- tmp %>% 
      get_segments(target = "sound") 
    }
    else{
      segments <- tmp %>% rename(segment = sound)
    }
    if(nrow(segments) < 10){
      return(NULL)
    }
    segments %>% check_for_power_law(target = "segment", with_plot = F) %>% 
      mutate(song_id = sid)
  }) %>% 
    mutate(segment_based = with_segments)
}