# Index of Skewness
library(touchard)
library(latex2exp)

r= 3 
lambdas = seq(1, 20, by=0.1)
deltas =   c(1,1.5,2)   #c(-0.5,0,0.5)  c(-1,-1.5,-2)

lines <- matrix(NA, nrow=length(lambdas), ncol=3)

for (del in deltas) {
  e_r_vec = c() #vec of values for each delta
  for (lambda in lambdas) {
    
    t0 <- tau(lambda,del)
    t1 <- tau(lambda,del+1)
    t2 <- tau(lambda,del+2)
    ex <- (t1/t0) - 1
    va <- (t2/t0) -(t1/t0)*(t1/t0)
    
    e_r = 0
    for (j in seq(0,r,by=1)) {
      e_r = e_r + 
            choose(r, j)*(
              ((-1)^(r-j))*tau(lambda=lambda, delta=(del+j)) / tau(lambda=lambda, delta=del)
                          )
    }
    skew = ( e_r -3*ex*va - ex^3)/ va^(3/2)
    e_r_vec = c(e_r_vec, skew)
  }
  lines[,which(del == deltas)] = e_r_vec
}

t1= paste0("$\\delta =",deltas[1],"$")
t2= paste0("$\\delta =",deltas[2],"$")
t3= paste0("$\\delta =",deltas[3],"$")
y_lab = paste0("$ E(X^",r,")$")

plot(lambdas, lines[,1], type = "l", lwd=2, font.axis = 2, 
     ylab = 'Skewness', xlab = TeX("$\\lambda$"), col = rgb(0,0,0), ylim = c(0.2,0.9))  
lines(lambdas, lines[,2], type = "l", col = rgb(1,0,0), lwd=2)        # Add second line
lines(lambdas, lines[,3], type = "l", col = rgb(0,1,0), lwd=2)   
legend("topleft",  bty="n",                                     # Add legend to plot
       legend = c(TeX(t1), TeX(t2), TeX(t3)),
       col = c(rgb(0,0,0),rgb(1,0,0),rgb(0,1,0)),
       pch = c(16, 16, 16))
