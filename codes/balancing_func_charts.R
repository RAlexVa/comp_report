library(tidyverse)
library(latex2exp)
##### unbounded BF #####
x <- seq(0,4,by=0.1)
data_u <- tibble(r=x,
                 min=sapply(x,min,1),
                 sq=sqrt(x),
                 max=sapply(x,max,1))

(plot1 <- data_u |> pivot_longer(-r,names_to = "b.fun",values_to = 'h(r)') |>  
  ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
  geom_line(size=1)+
  geom_segment(aes(x=2,y=0,xend=2,yend=2),color = "blue", linetype = "dashed", size = 1)+
  geom_segment(aes(x=3,y=0,xend=3,yend=3),color = "red", linetype = "dashed", size = 1)+
  annotate("text", x=2.1, y=0, label= TeX("$y_1$"),size=6, color='blue')+
  annotate("text", x=3.1, y=0, label= TeX("$y_2$"),size=6, color='red')+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=14)))
  
jpeg("C:/Users/ralex/Documents/src/comp_report/compare_balancing.jpeg", width = 850, height = 300)
plot1
dev.off()
##### Bounded BF #####
bounded_bf <- function(t,f,K=3){
  return(min(min(f(t),K),t*min(f(1/t),K))/K)
}
bounded_bf_log <- function(t,K=log(3)){
  return(min(min(t/2,K),t+min(-t/2,K)))
}

# exp(sapply(log(seq(0.1,2.5,by=0.1)),bounded_bf_log))
# sapply(seq(0.1,2.5,by=0.1),bounded_bf,f=sqrt)*3

f_c <- function(x,c=1.3){
  e1 <- min(1,x*exp(-c));
  e2 <- min(x,exp(-c));
  return(max(e1,e2));
}
sapply(x,bounded_bf,f=sqrt)
sapply(x,bounded_bf,f=function(x){max(x,1)})

data_u <- tibble(r=x,
                 min=sapply(x,min,1),
                 b_sq=sapply(x,bounded_bf,f=sqrt),
                 # b_max=sapply(x,bounded_bf,f=function(x){max(1,x)}),
                 # b_min=sapply(x,bounded_bf,f=function(x){min(3,x)}),
                 f_c=sapply(x,bf_c))

(plot2 <- data_u |> pivot_longer(-r,names_to = "b.fun",values_to = 'h(r)') |>  
    ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
    geom_line(size=1)+
    geom_segment(aes(x=2,y=0,xend=2,yend=1),color = "blue", linetype = "dashed", size = 1)+
    geom_segment(aes(x=3,y=0,xend=3,yend=1),color = "red", linetype = "dashed", size = 1)+
    annotate("text", x=2.1, y=0, label= TeX("$y_1$"),size=6, color='blue')+
    annotate("text", x=3.1, y=0, label= TeX("$y_2$"),size=6, color='red')+
    # labs(title='Bounded balancing functions')+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=14)))

jpeg("C:/Users/ralex/Documents/src/comp_report/compare_balancing_bounded.jpeg", width = 850, height = 300)
plot2
dev.off()

bounded_bf(5,3)
bounded_bf(10,3)

bounded_bf(1/5,3)*5
bounded_bf(1/10,3)*10
x <- seq(0,10,by=0.1)
y <- sapply(x,bounded_bf,K=3)/3
plot(x,y)
lines(x,y)


####### Otros ejemplos de bounded bf
b_f_min <- function(x,f,k){
  t1 <- min(f(x),k)
  t2 <- min(f(1/x),k)
  return (min(t1,x*t2)/k)
}
b_f_max <- function(x,f){
  return (max(f(x),x*f(1/x)))
}

x <- seq(0.1,25,by=0.1)

data_u <- tibble(r=x,
                 f1=sapply(x,b_f_min,f=sqrt,k=1),
                 f2=sapply(x,b_f_min,f=sqrt,k=1.5),
                 f3=sapply(x,b_f_min,f=sqrt,k=2),
                 f4=sapply(x,b_f_min,f=sqrt,k=3),
                 f5=sapply(x,b_f_min,f=sqrt,k=5),
                 f6=sapply(x,b_f_min,f=sqrt,k=10))

n <- 5
n.root <- function(a,n=5){a^(1/n)}
data_u <- tibble(r=x,
                 f1=sapply(x,b_f_min,f=n.root,k=1),
                 f2=sapply(x,b_f_min,f=n.root,k=1.5),
                 f3=sapply(x,b_f_min,f=n.root,k=2),
                 f4=sapply(x,b_f_min,f=n.root,k=3),
                 f5=sapply(x,b_f_min,f=n.root,k=5),
                 f6=sapply(x,b_f_min,f=n.root,k=10))

 data_u |> pivot_longer(-r,names_to = "b.fun",values_to = 'h(r)') |>  
    ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
    geom_line(size=1)

###### Bounded balancing functions #####
 bf <- function(x,f,bb){
   boundf <- min(f(x),bb)/bb
   boundf_1 <- min(f(1/x),bb)/bb
   
   return(min(boundf,x*boundf_1))
 }

 for(i in seq(0,15,by=0.1)){
   print(paste(round(bf(i,sqrt,4),4),round(bf(1/i,sqrt,4),4),round(i*bf(1/i,sqrt,4),10)==round(bf(i,sqrt,4),10)))
 }
 
 x <- seq(0.1,20,by=0.1)
 
 data_u <- tibble(r=x,
                  f1=sapply(x,bf,f=sqrt,bb=4),
                  f2=sapply(x,bf,f=sqrt,bb=2),
                  f3=sapply(x,bf,f=function(x){x^(2)},bb=2),
                  f4=sapply(x,bf,f=function(x){x^(1/4)},bb=2),
                  f5=sapply(x/16,min,1),
                  f6=sapply(x/5,min,1))
 
 data_u |> pivot_longer(-r,names_to = "b.fun",values_to = 'h(r)') |>  
   ggplot(aes(x=r,y=`h(r)`,color = `b.fun`))+
   geom_line(size=1)

 # Using a power function doesnt seems to give good results
 #But perhaps we can just use min function with a different bounding constant
 #Here I think we can choose how important we want to make each neighbors
 
 
 