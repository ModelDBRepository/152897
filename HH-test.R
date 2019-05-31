# Hodgkin-Huxley kinetics
library(deSolve)

LVmod0D <- function(Time, State, Pars) {
 with(as.list(c(State, Pars)), {
if(Time>=5&&Time<=5.5){Amp}else{Amp<-0}
taum<-1/(0.1*(V+40)/(1-exp(-(V+40)/10))+4*exp(-(V+65)/18))
minf<-(0.1*(V+40)/(1-exp(-(V+40)/10)))*taum
tauh<-1/(0.07*exp(-(V+65)/20)+(1/(1+exp(-(V+35)/10))))
hinf<-0.07*exp(-(V+65)/20)*tauh
taun<-1/(0.01*(V+55)/(1-exp(-(V+55)/10))+0.125*exp(-(V+65)/80))
ninf<-(0.01*(V+55)/(1-exp(-(V+55)/10)))*taun

dV<-(-gna*(m^3)*h*(V-50)-gk*(n^4)*(V-(-77))-gl*(V-(-54.4))+Amp)/Cm
dm<--(m-minf)/taum
dh<--(h-hinf)/tauh
dn<--(n-ninf)/taun

return(list(c(dV, dm,dh,dn)))
 })
 }
pars<-c(gna=120,gk=36,gl=0.3, Cm=1, Amp=20)

yini<-c(V=-65, m=0.052, h=0.596, n=0.317)

times <- seq(0, 20, by = 0.1)

print(system.time(
out <- ode(func = LVmod0D, y = yini, parms = pars, times = times)))

## the model output is plotted, using R function matplot
## http://www.jstatsoft.org/v33/i09/paper
matplot(out[,"time"], out[,2:3], type = "l", xlab = "time", ylab = "Membrane potential",
main = "Hodgkin-Huxley model", lwd = 2, col="blue")
# legend("topright", c("Membrane potential", "m"), col = 1:2, lty = 1:2)

