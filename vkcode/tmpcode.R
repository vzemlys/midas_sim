####How exactly %:% works?
###The answer is the curly braces. The following statements are identical:
set.seed(1)
rez.01 <- foreach(n.s=1:4,.combine='rbind')%do% {
          foreach(i=1:2,.combine='rbind') %do% {
              rnorm(1)
          }}
set.seed(1)
rez.00 <- foreach(n.s=1:4,.combine='rbind')%:% 
          foreach(i=1:2,.combine='rbind') %do% {
              rnorm(1)
          }
