### These functions are used in data_prep files

## Function to get most discrete value
most_frequent_discrete_value <- function(vec) {
  tbl <- table(vec)
  if (dim(tbl) > 0) {
    value <- names(tbl)[which.max(tbl)]
    return(value)
  } else {
    return(NA)
  }
}

## Function to create species name
get_gen_sp <- function(x) {
  if (is.na(x)) {
    return(NA)
  } else if (grepl("_cf_", x)) {
    namesplit <- unlist(strsplit(x, split = "_"))
    newname <- paste0(namesplit[1], "_", namesplit[3])
    return(newname)
  } else {
    namesplit <- unlist(strsplit(x, split = "_"))
    newname <- paste0(namesplit[1], "_", namesplit[2])
    return(newname)
  }
}


## Function to quantify branching and bracteoles
typeset <- function(df) {
  if (any(is.na(df))) {
    return(NA)
  } else {
    if (df[2]==TRUE & df[3]==FALSE) { # umbellike, no bracteoles
      return(0)
    }
    else if (df[2]==TRUE & df[3]==TRUE) { # umbellike w/ bracteoles
      return(1)
    }
else if (df[2]==FALSE & df[3]==TRUE) { # non umbel (branching) w/ bracteoles
      return(2)
    } else {
      return(NA)
    }
  }
}
