my_reset_outliers <- function(x, na.rm = TRUE ) { ## x is a numerical vector or matrix.
  qnt <- quantile(x, probs=c(0.05, 0.95) , type=1,  na.rm = na.rm )  
  y <- x
  y[x < qnt[1] ] <- qnt[1]
  y[x > qnt[2] ] <- qnt[2]    
  return(y)
}



## scale values to [lower, upper]
my_scale_matrix <- function( x, upper = 1, lower = -1 ) { ## x is a numerical vector or matrix. 
    y = my_reset_outliers(x, na.rm = TRUE)  
    z = lower + (upper - lower) * ( y - min(y) )/( max(y)- min(y) )
    return(z)
}

