library(touchard)
library(rootSolve)
library(ggplot2)
library(latex2exp)

csq = qchisq(0.95,1) #chisqaure quantile
n_sim = 1000 #Number of simulations
h_vec = seq(1,2.6,by=0.2) #values for h

#Values of N and Delta
N_vec <- c(5, 8, 10, 14, 16)  #c(20,35,50,75, 100)
delta_vec = c(-1.5, -1, 1, 1.5)

#solves for lambda for given mu
solve_lambda <- function(mu_in, delta) {
  f <- function (lbd) (tau(lbd,delta+1)/tau(lbd, delta)) - 1 - mu_in
  ld <- uniroot.all(f, c(0,100))
  return(ld)
}

#solves for variance
t_var <- function(l_in,d_in) {
  v_out <-   (tau(l_in,d_in+2)/tau(l_in,d_in)) -(tau(l_in,d_in+1)/tau(l_in,d_in))*(tau(l_in,d_in+1)/tau(l_in,d_in))
  return(v_out)
}


for ( N in N_vec) {
  for ( delta in delta_vec) {
    #############
    r_rej_vec_dev = c()
    r_rej_vec_sco = c()

    for ( h in h_vec) { #For loop for each h
      
      lambda_1 <- 6
      lambda_2 <- lambda_1 * h
      nrej_dev = 0
      nrej_sco = 0
      
      for ( i in seq(1,n_sim,by=1)) { #The 1000 Simulations for loop
        
        #The two random samples
        Y_1 <- rtouch(n = N, lambda=lambda_1, delta=delta)
        Y_2 <- rtouch(n = N, lambda=lambda_2, delta=delta)
        Y <-data.frame(Y_1,Y_2)
        
        #######Hypothesis############
        #Null Hypothesis
        mu_hat = (1/(2*N))*sum(Y)
        ld <- solve_lambda(mu_hat, delta)
        
        #Alternative Hypothesis
        mu_hat2 = (1/N)*sum(Y$Y_2)
        mu_hat1 = ((1/N)*sum(Y)) - mu_hat2
        
        ld_1 <- solve_lambda(mu_hat1, delta)
        ld_2 <- solve_lambda(mu_hat2, delta)
        
        #############Score######################################
        ######score######
        U = c(0, sum(Y$Y_2)-(N*mu_hat))
        J_inv = solve( t_var(ld, delta) * rbind(c(2*N,N),c(N,N)) )
        score = t(U)%*% J_inv %*% U
        
        #############Deviance######################################
        ##DElta D###
        y = c(as.matrix(Y))
        del_D = 2 *(sum(log( dtouch(y,  lambda = c(rep(ld_1,N),rep(ld_2,N)), delta = delta))) 
                    - 
                      sum(log( dtouch(y,  lambda =ld,     delta = delta)))
        )
        
        
        # Verify if delta_D is greater than the chiSquare 
        if(del_D  > csq ){
          nrej_dev = nrej_dev + 1
        } else {
          #do nothing
        }
        
        # Verify if score is greater than the chiSquare 
        if(score  > csq ){
          nrej_sco = nrej_sco + 1
        } else {
          #do nothing
        }
      }
      
      r_rej_dev = nrej_dev / n_sim
      r_rej_vec_dev = c(r_rej_vec_dev, r_rej_dev)
      
      r_rej_sco = nrej_sco / n_sim
      r_rej_vec_sco = c(r_rej_vec_sco, r_rej_sco)
      
    }
    
    #Plot data
    df <- data.frame(Statistics = c(rep("Score", each=length((h_vec))),rep("Deviance", each=length((h_vec)))),
                     h=c(h_vec,h_vec),
                     Power=c( r_rej_vec_sco , r_rej_vec_dev) )
    #Plot title and legend
    tt= paste0("$\\delta =",delta," , N =",N," , \\lambda_{2} = \\lambda_{1} * h "," , \\lambda_{1} =",lambda_1,"$")
    img_plot <- ggplot(df, aes(x=h, y=Power, group=Statistics)) + 
      geom_line(aes(linetype=Statistics, color=Statistics), size=1.2)+
      coord_cartesian(ylim=c(0, 1)) + 
      geom_hline(yintercept=0.05, linetype="dashed") +
      theme_bw() +
      ggtitle(TeX(tt))
    #Saves plots to file for each plot
    ggsave(img_plot, file=paste0("plot_del", delta,"N",N,".png"), width = 15.61, height = 11.22, units = "cm")
    

  }
}
