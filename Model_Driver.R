rm(list=ls(all=TRUE))
library(deSolve)
library(ggplot2)
library(gridExtra)

##Function parameters
## "constant" describes arbritray link between dying worms & optical density
parameters<- c(R0=2,a1=0.9,a2=15, a3=0.6, a4=0.1, mu.w=1/5, mu.h=0.04,c=88.95, beta=0.05, 
               d=15 , nw = 5, na=50, amax=50, gain.imm=.1, loss.imm = .2, eggrate=5.26)

## immunity function (depends on shape parameters & constant linking number of dying worms)
## therefore dynamically updated during numerical integration
imm.a<-function(parameters, immuno.stimulus) {
  with(as.list(c(parameters)), {
    j = seq(1,na)
    a = (j-1)*amax/na + amax/(2*na)
    da <- amax/na
    p.estab <- a1/(1+exp(a2*(immuno.stimulus-a3)))+a4  
    #p.estab <- exp(-a1*immuno.stimulus) + 
    p.estab# this returns the probability of establishment 
    # which varies between 1 and 0.1 (this will need to be revised)
  })
}

## demographic age structure
pi.a <- function(parameters) {
  with(as.list(c(parameters)), {
    j = seq(1,na)
    a = (j-1)*amax/na + amax/(2*na)
    da <- amax/na
    da*mu.h*exp(-mu.h*a)/(1-(exp(-mu.h*amax))) 
  })
}

## cercarial exposure by age
rho.a<- function(parameters) {
  with(as.list(c(parameters)), {
    j = seq(1,na)
    a = (j-1)*amax/na + amax/(2*na)
    a*c*exp(-(beta*a)) +d
    ## set to constant expore
    d
  })
}

## Transmission model
model <- function(times, y, parameters){
  with(as.list(c(y,parameters)),{
    ymat <- matrix(y, ncol=na, nrow=(nw+2))
    dymat <- matrix(0, ncol=na, nrow=(nw+2))
    ## indicator variables for elements of i corresponding to live, dead and optical density IgE
    live.worms.id <- seq(1:nw)
    dead.worms.id <- nw+1
    IgE.id <- nw + 2 # this compartment tracks the IgE optical density in different age groups 
    tot.comp <- IgE.id # total number of compartments
    #age step
    da <- 1/(amax/na) 
    # IgE OD, either a function of dying worms or cumulative experience of dead worms
    IgE.a <- ymat[IgE.id, ]
    #establishment probability by each age group
    probest.a <- imm.a(parameters, IgE.a) 
    pia <- pi.a(parameters) # age structure
    rhoa <- rho.a(parameters) # exposure function
    Wtot <- colSums(ymat[live.worms.id,]) #Total number of worms
    lambda.a <- (mu.w*R0*rhoa*probest.a*sum(pia*rhoa*Wtot))/sum(pia*rhoa^2) #force of infection by age
    for (j in 1:na) { #'for' loop (first age group to total number (na))
      if (j==1) { #first age group
        for (i in 1:tot.comp) { #for worm compartment (i) in worm compartments W0:D (within the age loop)
          if (i==1) { #i.e. W0
            #Tranmisssion equation for W0: foi - worm progression from this compartment - aging from this age strata
            dymat[i,j] = lambda.a[j]-(nw*mu.w + da)*ymat[i,j] 
          } else if (i<dead.worms.id) { #i.e. W1..Wi
            #Transmission equation for W1..Wi worm progression from previous compartment - worm progression from this compartment - aging from this age strata
            dymat[i,j] = nw*mu.w*ymat[i-1,j] - (nw*mu.w + da)*ymat[i,j] 
          } else if (i==dead.worms.id) { #i.e. D
            #Transmission equation for D: gain from previous worm compartment - aging from this age strata
            dymat[i,j] = nw*mu.w*ymat[nw,j] - da*ymat[i,j]
          }
          else if (i==IgE.id) { #i.e. Imm
            #Transmission equation for immunity: product of dying worms and a constant 
            dymat[i,j] = gain.imm*ymat[nw,j]*nw*mu.w - ymat[i,j]*loss.imm - ymat[i,j] *da
          }
        }
      } 
      #previously = else if (j>1), changed to (j>1&j<na)
      else if (j>1&j<na) { #for age classes j+1..na-1
        for (i in 1:tot.comp) {
          if (i==1) { #i.e. W0
            #Tranmisssion equation for W0: foi - worm progression - aging from this age strata + aging from previous age strata
            dymat[i,j] = lambda.a[j]-(nw*mu.w +da)*ymat[i,j] + ymat[i,j-1]*da 
          } else if (i<dead.worms.id) { #i.e. W1..Wi
            #Transmission equation for W1..Wi: worm progression from previous compartment - worm progression from this compartment - aging from this compartment + aging from previous age strata
            dymat[i,j] = nw*mu.w*ymat[i-1,j] - nw*mu.w*ymat[i,j] + (ymat[i,j-1] - ymat[i,j])*da 
          } else if (i==dead.worms.id) { #i.e. D
            #Transmission equation for D: worm progression from previous compartment - aging from this compartment +               aging from previous age strata
            dymat[i,j] = nw*mu.w*ymat[nw,j] + (ymat[i,j-1] - ymat[i,j])*da
          }
          else if (i==IgE.id) { #i.e. Imm
            #Transmission equation for immunity: product of dead worm compartment and constant
            #this component is not used
            dymat[i,j] = gain.imm*ymat[nw,j]*nw*mu.w  - ymat[i,j]*loss.imm + (ymat[i,j-1] - ymat[i,j])*da
          }
        } 
        #previously = else if (j>na), changed to (j==na) - could also be just 'else'?
      } else if (j==na) { #i.e. eldest age class 
        for (i in 1:tot.comp) {
          if (i==1) { #i.e. W0
            #Tranmisssion equation for W0: foi - worm progression + aging from previous age strata
            dymat[i,j] = lambda.a[j]-(nw*mu.w)*ymat[i,j] + (ymat[i,j-1] - ymat[i,j])*da
          } else if (i<dead.worms.id) { #i.e. W1..Wi
            #Transmission equation for W1..Wi: worm progression from previous compartment - worm progression from this compartment + aging from previous age strata
            dymat[i,j] = nw*mu.w*ymat[i-1,j] - nw*mu.w*ymat[i,j] + (ymat[i,j-1] - ymat[i,j])*da
          } else if (i==dead.worms.id) { #i.e. D
            #Transmission equation for D: worm progression from previous compartment + aging from previous age strata
            dymat[i,j] = nw*mu.w*ymat[nw,j] + (ymat[i,j-1] - ymat[i,j])*da
          }
          else if (i==IgE.id) { #i.e. Imm
            #Transmission equation for immunity: product of dead worm compartment and constant 
            # this component is not used
            dymat[i,j] = gain.imm*ymat[nw,j]*nw*mu.w  - ymat[i,j]*loss.imm + (ymat[i,j-1] - ymat[i,j])*da
          }
        }
      }
    }   
    
    ## return relevant variables alondside matrix output
    worms.a  <- colSums(ymat[1:nw,])
    worms.mn <- weighted.mean(worms.a, pia)
    
    ## egg output (no density dependence  - this may need to be revised)
    epg.a <- worms.a*0.5*eggrate
    epg.mn <- weighted.mean(epg.a, pia)
    
    dead.worms.a <- ymat[dead.worms.id,]
    dead.worms.mn <- weighted.mean(dead.worms.a, pia)
    
    probest.mn <- weighted.mean(probest.a, pia)
    IgE.mn <- weighted.mean(IgE.a, pia)
    
    return(list(rbind(dymat), 
                worms.age=worms.a, worms.mn=worms.mn, IgE.a=IgE.a, IgE.mn=IgE.mn,
                dead.worms.mn = dead.worms.mn, epg.a = epg.a ,
                epg.mn = epg.mn, probest.mn=probest.mn))
  })
}


runmod <- function(parameters) {
  y <- rep(1, (parameters["nw"]+2)*parameters["na"])
  ## stepsize and times
  stepsize=60/365
  ## 50 years to equilibrium
  times<-seq(0,100, stepsize)
  out<-rk4(y=y, times=times, func=model, parms=parameters)
  out <- as.data.frame(out)
  out
}

## function to extract age dependency at equilibrium
extract <- function(out, parameters) {
  tmp <- t(out[nrow(out), grep("worms.age", colnames(out))])
  tmp <- cbind(seq(1, parameters["na"]), tmp[,1])
  tmp <- as.data.frame(cbind(tmp, t(out[nrow(out), grep("IgE.a", colnames(out))]), 
                             t(out[nrow(out), grep("epg.a", colnames(out))])))
  colnames(tmp) <- c("age", "worms", "IgE", "epg")
  tmp
}

# test <- extract(out)
