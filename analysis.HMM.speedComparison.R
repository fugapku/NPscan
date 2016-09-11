# 
require('HMM') # for discrete
require('seqHMM') # for discrete 
require('HiddenMarkov')
require('mhsmm') # for continous 
require('depmixS4') #79
require('msm') # 260
#RHmm #2
#HMMmix

# Sequence of observation
n = 100
a = sample(c(rep("L",100*n),rep("R",300*n)))
b = sample(c(rep("L",300*n),rep("R",100*n)))
observation = c(a,b)

### initialize
hmm = initHMM(c("A","B"),c("L","R"),startProbs = c(0.5,0.5),
              transProbs=matrix(c(.9,.1,.1,.9),2),
              emissionProbs=matrix(c(.7,.3,.3,.7),2))
print(hmm)
seq.formated.s <- seqdef(t(observation), 1:length(observation),
                         labels = c('L','R'))
HMM.s = build_hmm(state_names = c("A","B"), 
                  observations = seq.formated.s,
                  initial_probs = c(0.5,0.5),
                  transition_probs = matrix(c(.9,.1,.1,.9),2),
                  emission_probs = matrix(c(.5,.51,.5,.49),2))


# Viterbi-training: HMM and seq HMM similar in speed
ptm <- proc.time()
vt = viterbiTraining(hmm,observation,1000)
proc.time()-ptm
# BW-training
ptm <- proc.time()
vt.bw = baumWelch(hmm,observation,1000)
proc.time()-ptm
# seq HMM
ptm <- proc.time()
fit.HMM <- fit_model(HMM.s,control_em = list(maxeval = 1000, restart = list(times = 0)),
                     global_step=T, control_global = list(maxtime=1000),
                     local_step=T)
proc.time()-ptm


## decoding: HMM is faster than seqHMM
ptm <- proc.time()
statePred = viterbi(vt.bw$hmm, observation)
postP = posterior(vt.bw$hmm, observation)
proc.time()-ptm

ptm <- proc.time()
statePred = hidden_paths(fit.HMM$model)
statePred1 = as.character(unlist(as.list(statePred)))
postP = posterior_probs(fit.HMM$model)
postP1 = postP[,,1]
proc.time()-ptm

print(vt$hmm)
print(vt.bw$hmm)
print(seqHMM2HMM(fit.HMM$model))

#### depmixS4
n = 10
a = sample(c(rep("L",100*n),rep("M",100*n),rep("R",300*n)))
b = sample(c(rep("L",300*n),rep("M",100*n),rep("R",100*n)))
observation = c(a,b)

hmm = initHMM(c("A","B"),c("L", 'M', "R"),startProbs = c(1,0),
              transProbs=t(matrix(c(.999,.001,0.1, 0.9),2)),
              emissionProbs=t(matrix(c(.20,.20,.6, .60,0.20,.20),3)));
mod <- HMM2depmix(hmm, seqs = observation)
ptm <- proc.time()
fm <- fit(mod)
proc.time()-ptm
summary(fm)

hmm = initHMM(c("A","B"),c("L", 'M', "R"),startProbs = c(0,1),
              transProbs=t(matrix(c(.95,.05,.1,.9),2)),
              emissionProbs=t(matrix(c(.7,.2,.1, .2,0.3,.5),3)));
observation1 = simHMM(hmm, 30000)
dat = data.frame(obs = observation1$observation)
mod <- HMM2depmix(hmm, seqs = observation1$observation)
set.seed(1)
ptm <- proc.time()
fm <- fit(mod)
proc.time()-ptm
summary(fm)
mod1 <- HMM2depmix(depmix2HMM(fm1), seqs = observation1$observation)
set.seed(1)
ptm <- proc.time()
fm1 <- fit(mod1)
proc.time()-ptm
summary(fm1)

set.seed(1)
fm1 <- fit(fm1, emc = em.control(rand=F))
summary(fm1)

if (0){
  rModels <- list(list(GLMresponse(obs~1, data=dat, family = multinomial('identity'), pstart = hmm$emissionProbs[1,])),
                  list(GLMresponse(obs~1, data=dat, family = multinomial('identity'), pstart = hmm$emissionProbs[2,])))
  transition <- list()
  transition[[1]] <- transInit(~1,nstates=2,data=dat, pstart = hmm$transProbs[1,])
  transition[[2]] <- transInit(~1,nstates=2,data=dat, pstart = hmm$transProbs[2,])
  inMod <- transInit(~1,ns=2,family=multinomial("identity"), pstart = hmm$startProbs)
  mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,
                    homogeneous=FALSE)
}

mod <- HMM2depmix(hmm, seqs = observation$observation)
depmix2HMM(mod)
mod1 <- depmix(list(obs~1), data=dat,nstates=2,
              respstart = hmm1$emissionProbs, trstart = hmm1$transProbs, instart = hmm1$startProbs,
              family=list(multinomial('identity')))

#set.seed(1)
ptm <- proc.time()
fm <- fit(mod)
proc.time()-ptm
summary(fm)
ptm <- proc.time()
fm1 <- fit(mod1)
proc.time()-ptm
summary(fm1)

fm@transition
fm@response
#posterior(fm)


data(speed)	
mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,
              family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
# print the model, formulae and parameter values
mod
set.seed(1)
# fit the model by calling fit
fm <- fit(mod)
posterior(fm)


# Volatility of S & P 500 returns
# (thanks to Chen Haibo for providing this example)

data(sp500)

# fit some models
msp <- depmix(logret~1,nstates=2,data=sp500)
set.seed(1)
fmsp <- fit(msp)
posterior(fmsp)
plot(ts(posterior(fmsp)[,2], start=c(1950,2),deltat=1/12),ylab="probability",
     main="Posterior probability of state 1 (volatile, negative markets).",
     frame=FALSE)


HMM
x <- rnorm(1000)
p <- plogis(x)
ss <- rbinom(1000,1,p)
mod <- GLMresponse(cbind(ss,1-ss)~x,family=binomial())
fitted = fit(mod)
posterior(fitted)
glm(cbind(ss,1-ss)~x, family=binomial)

