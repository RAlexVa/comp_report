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
  geom_segment(aes(x=0.5,y=0,xend=0.5,yend=1),color = "blue", linetype = "dashed", size = 1)+
  geom_segment(aes(x=3.5,y=0,xend=3.5,yend=3.5),color = "red", linetype = "dashed", size = 1)+
  annotate("text", x=0.7, y=0, label= TeX("$y_1$"),size=6, color='blue')+
  annotate("text", x=3.7, y=0, label= TeX("$y_2$"),size=6, color='red')+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=14)))
  
jpeg("C:/Users/ralex/Documents/src/comp_report/compare_balancing.jpeg", width = 850, height = 300)
plot1
dev.off()
##### Bounded BF #####
bounded_bf <- function(t,f,K=3){
  return(min(min(f(t),K),t*min(f(1/t),K)))
}
sapply(x,bounded_bf,f=sqrt)
sapply(x,bounded_bf,f=function(x){max(x,1)})

data_u <- tibble(r=x,
                 min=sapply(x,min,1),
                 # sq=sqrt(x),
                 max=sapply(x,max,1),
                 # f_c=sapply(x,bf_c),
                 b_f=sapply(x,bounded_bf,f=sqrt))


bounded_bf(5,3)
bounded_bf(10,3)

bounded_bf(1/5,3)*5
bounded_bf(1/10,3)*10
x <- seq(0,10,by=0.1)
y <- sapply(x,bounded_bf,K=3)/3
plot(x,y)
lines(x,y)


