library(touchard)
library(latex2exp)


lambdas = seq(1, 20, by=0.1) #values of lambda
deltas =  c(-1,-1.5,-2)   # c(1,1.5,2) c(-0.5,0,0.5)  values of delta

lines <- matrix(NA, nrow=length(lambdas), ncol=3) #store values in a matrix

for (del in deltas) {
  
  e_r_vec = c() #vec of values for each delta
  for (lambda in lambdas) {
    t0 <- tau(lambda,del) #Tau
    t1 <- tau(lambda,del+1)
    t2 <- tau(lambda,del+2)
    ex <- (t1/t0) - 1 #mu
    va <- (t2/t0) -(t1/t0)*(t1/t0) #variance
    #Fourth moment
    r = 4
    e_r_fth = 0
    for (j in seq(0,r, by=1)) {
      e_r_fth = e_r_fth + 
        choose(r, j)*(
          ((-1)^(r-j))*tau(lambda=lambda, delta=(del+j)) / tau(lambda=lambda, delta=del)
        )
    }
    #Third moment
    r = 3
    e_rthird = 0
    for (j in seq(0, r ,by=1)) {
      e_rthird = e_rthird + 
        choose(r, j)*(
          ((-1)^(r-j))*tau(lambda=lambda, delta=(del+j)) / tau(lambda=lambda, delta=del)
        )
    }
    
    kurt = (e_r_fth -4*e_rthird*ex +6*( va +ex^2)*ex*ex - 3* ex^4) / va^2 #Index of kurtosis
    e_r_vec = c(e_r_vec, kurt)
  }
  lines[,which(del == deltas)] = e_r_vec
}

t1= paste0("$\\delta =",deltas[1],"$") #Plot legends
t2= paste0("$\\delta =",deltas[2],"$")
t3= paste0("$\\delta =",deltas[3],"$")
y_lab = paste0("$ E(X^",r,")$")
#Plot data
plot(lambdas, lines[,1], type = "l", lwd=2, font.axis = 2, 
     ylab = 'Kurtosis', xlab = TeX("$\\lambda$"), col = rgb(0,0,0))  #, ylim = c(0,200)     
lines(lambdas, lines[,2], type = "l", col = rgb(1,0,0), lwd=2)      
lines(lambdas, lines[,3], type = "l", col = rgb(0,1,0), lwd=2)   
legend("topright",  bty="n",                                     # Add legend to plot
       legend = c(TeX(t1), TeX(t2), TeX(t3)),
       col = c(rgb(0,0,0),rgb(1,0,0),rgb(0,1,0)),
       pch = c(16, 16, 16))