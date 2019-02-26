library(rpart)
source("MedTrees.R")


id = c()
x = c()
m = c()
y = c()

Z1 = c()
Z2 = c()
Z3 =c()

a = 0
b = 0
cP = 0

a_2 = 1
b_2 = 1
cP_2 = 1

ai = c()
bi = c()
cPi = c()

node = c()

n = 300

for(i in 1:n){
  id[i] = i
  Z1[i] = round(rnorm(1))
  Z2[i] = round(rnorm(1))
  Z3[i] = round(rnorm(1))
  
  if(Z1[i] < 0){
    node[i] = 1
    x[i] = round(runif(1))
    m[i] = x[i]*a + .1*rnorm(1)
    y[i] = x[i]*cP + m[i]*b + .1*rnorm(1)
    
    ai[i] = a
    bi[i] = b
    cPi[i] = cP
    
  }else{
    node[i] = 2
    x[i] = round(runif(1))
    m[i] = x[i]*a_2 + .1*rnorm(1)
    y[i] = x[i]*cP_2 + m[i]*b_2 + .1*rnorm(1)
    
    ai[i] = a_2
    bi[i] = b_2
    cPi[i] = cP_2
  }
}

SimDATA = cbind.data.frame(id,x,m,y,Z1,Z2,Z3)

TREE = MedTrees(totalXY = y ~ x, 
                    bc_effect = y ~ x + m, 
                    a_effect = m ~ x, 
                    data = SimDATA, 
                    group = ~ id,
                    rPartFormula= ~ Z1 + Z2 + Z3,
                    control=rpart.control(cp=.05))

TREE$rpart_out

TREE$summary_XM
TREE$summary_XY
TREE$summary_XMY

library(rpart.plot)
prp(TREE$rpart_out)
