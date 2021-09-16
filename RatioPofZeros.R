#Ratio of Probabilities of zeros 
library(touchard)
library(latex2exp)

lambda = seq(1, 20, by=0.01)
delta =c(1,1.5,2)# c(-0.5,0,0.5) c(-1,-1.5,-2)  

p_zero1 = dtouch(0, lambda=lambda, delta=delta[1]) / dpois(0, lambda=lambda)
p_zero2 = dtouch(0, lambda=lambda, delta=delta[2]) / dpois(0, lambda=lambda)
p_zero3 = dtouch(0, lambda=lambda, delta=delta[3]) / dpois(0, lambda=lambda)

t1= paste0("$\\delta =",delta[1],"$")
t2= paste0("$\\delta =",delta[2],"$")
t3= paste0("$\\delta =",delta[3],"$")

plot(lambda, p_zero1, type = "l", lwd=2, font.axis = 2, 
     ylab = "P(0) Ratio", xlab = "Lambda", col = rgb(0,0,0), ylim = c(0,0.8))    #, ylim = c(0,200) # Draw first line
lines(lambda, p_zero2, type = "l", col = rgb(0,0,1), lwd=2)        # Add second line
lines(lambda, p_zero3, type = "l", col = rgb(0,1,0), lwd=2)   
legend("topleft",  bty="n",                                     # Add legend to plot
       legend = c(TeX(t1), TeX(t2), TeX(t3)),
       col = c(rgb(0,0,0),rgb(0,0,1),rgb(0,1,0)),
       pch = c(16, 16, 16))
