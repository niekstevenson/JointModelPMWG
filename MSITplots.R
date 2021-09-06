library(tidyr)

df <- aggregate(correct ~ resp*flank, dat, mean)
ggplot(df[df$correct > 0.6,], aes(x = as.factor(flank), y = correct, color = as.factor(resp))) + geom_point(size = 3) + theme_bw()


df <- aggregate(correct ~ resp*pos, dat[dat$flank == 0,], mean)
ggplot(df[df$correct > 0.6,], aes(x = as.factor(pos), y = correct, color = as.factor(resp))) + geom_point(size = 3) + theme_bw()

dat$resp[dat$pos == 3 & dat$correct == 0] %>% table()


df <- dat[!is.na(dat$rt),] %>% count(pos, resp, correct)

ggplot(df[df$correct == F,], aes(x = as.factor(pos), y = n, color = as.factor(resp))) + geom_point(size = 3) + theme_bw()



df <- dat[!is.na(dat$rt),] %>% count(flank, resp, correct)


ggplot(df[df$correct == F,], aes(x = as.factor(flank), y = n, color = as.factor(resp))) + geom_point(size = 3) + theme_bw()




df <- aggregate(correct ~ pos*resp*flank, dat[dat$condit == 2,], mean)
ggplot(df[df$correct > 0.1,], aes(x = as.factor(pos), y = correct, color = as.factor(resp), label = flank)) + 
         geom_text(aes(fontface = 2), size = 4, position = position_dodge(width = .2)) + 
         theme_bw()


df <- aggregate(correct ~ pos, dat, mean)
ggplot(df, aes(x = as.factor(pos), y = correct)) + geom_point(size = 3) + theme_bw()





# relevant shit -----------------------------------------------------------

#Flanker only condition
df <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 3,] %>% count(flank, resp)
n <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 3,] %>% count(flank)
df$sum <- rep(n[,2], each = 3)
df$errorsPerc <- df$n/df$sum
ggplot(df, aes(x = as.factor(flank), y = errorsPerc, color = as.factor(resp))) + 
  geom_point(size = 3) + 

  ggtitle(label = "Flanker only") + 
  theme_bw()


#Simon only condition
df <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 2,] %>% count(pos, resp)
n <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 2,] %>% count(pos)
df$sum <- rep(n[,2], each = 3)
df$errorsPerc <- df$n/df$sum
ggplot(df, aes(x = as.factor(pos), y = errorsPerc, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ylim(c(0,1)) + 
  ggtitle(label = "Simon only") + 
  theme_bw()

#Flanker for condition 4
df <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 4,] %>% count(flank, resp)
n <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 4,] %>% count(flank)
df$sum <- rep(n[,2], each = 3)
df$errorsPerc <- df$n/df$sum
ggplot(df, aes(x = as.factor(flank), y = errorsPerc, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ylim(c(0,1)) + 
  ggtitle(label = "Simon + Flanker") + 
  theme_bw()


#Simon for condition 4
df <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 4,] %>% count(pos, resp)
n <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 4,] %>% count(pos)
df$sum <- rep(n[,2], each = 3)
df$errorsPerc <- df$n/df$sum
ggplot(df, aes(x = as.factor(pos), y = errorsPerc, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ylim(c(0,1)) + 
  ggtitle(label = "Simon + Flanker") + 
  theme_bw()

#General flanker, regardless of condition
df <- dat[!is.na(dat$rt) & !dat$correct,] %>% count(flank, resp)
n <- dat[!is.na(dat$rt) & !dat$correct,] %>% count(flank)
df$sum <- rep(n[,2], each = 3)
df$errorsPerc <- df$n/df$sum
ggplot(df, aes(x = as.factor(flank), y = errorsPerc, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ylim(c(0,1)) + 
  theme_bw()

#General Simon, regardless of condition
df <- dat[!is.na(dat$rt) & !dat$correct,] %>% count(pos, resp)
n <- dat[!is.na(dat$rt) & !dat$correct,] %>% count(pos)
df$sum <- rep(n[,2], each = 3)
df$errorsPerc <- df$n/df$sum
ggplot(df, aes(x = as.factor(pos), y = errorsPerc, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ylim(c(0,1)) + 
  theme_bw()


#General flanker and Simon
df <- dat[!is.na(dat$rt) & !dat$correct,] %>% count(flank, pos, resp)
n <- dat[!is.na(dat$rt) & !dat$correct,] %>% count(flank, pos)
df$sum <- rep(n[,3], each = 3)
df$errorsPerc <- df$n/df$sum

ggplot(df, aes(x = as.factor(pos), y = errorsPerc, color = as.factor(resp), label = flank)) + 
  geom_text(aes(fontface = 2), size = 5, position = position_dodge(width = .25)) + 
  ylim(c(0,1)) + 
  theme_bw()


#Condition 4 only
df <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 4,] %>% count(flank, pos, resp)
n <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 4,] %>% count(flank, pos)
df$sum <- rep(n[,3], c(3,2,2,2,3,2,2,2,3))
df$errorsPerc <- df$n/df$sum

ggplot(df, aes(x = as.factor(pos), y = errorsPerc, color = as.factor(resp), label = flank)) + 
  geom_text(aes(fontface = 2), size = 5, position = position_dodge(width = .25)) + 
  ylim(c(0,1)) + 
  ggtitle(label = "Simon + Flanker") + 
  theme_bw()

ggplot(df, aes(x = as.factor(flank), y = errorsPerc, color = as.factor(resp), label = pos)) + 
  geom_text(aes(fontface = 2), size = 5, position = position_dodge(width = .25)) + 
  ylim(c(0,1)) + 
  ggtitle(label = "Simon + Flanker") + 
  theme_bw()



df <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 4,] %>% count(uniq, resp)
n <- dat[!is.na(dat$rt) & !dat$correct & dat$condit == 4,] %>% count(uniq)
df$sum <- rep(n[,2], each = 2)
df$errorsPerc <- df$n/df$sum

ggplot(df, aes(x = as.factor(uniq), y = errorsPerc, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ylim(c(0,1)) + 
  ggtitle(label = "Errors in condit 4") + 
  theme_bw()


df <- aggregate(rt ~ resp*condit, dat[dat$correct,], mean)
ggplot(df, aes(x = as.factor(condit), y = rt, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ggtitle(label = "correct RTs per condit") + 
  theme_bw()

df <- aggregate(rt ~ resp*condit, dat[!dat$correct,], mean)
ggplot(df, aes(x = as.factor(condit), y = rt, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ggtitle(label = "error RTs per condit") + 
  theme_bw()

df <- aggregate(correct ~ resp*condit, dat, mean)
ggplot(df, aes(x = as.factor(condit), y = correct, color = as.factor(resp))) + 
  geom_point(size = 3) + 
  ggtitle(label = "%correct per condit") + 
  theme_bw()

df <- aggregate(rt ~ pos*condit, dat[dat$correct,], mean)
ggplot(df, aes(x = as.factor(condit), y = rt, color = as.factor(pos))) + 
  geom_point(size = 3) + 
  ggtitle(label = "correct RTs per condit") + 
  theme_bw()

df <- aggregate(rt ~ pos*condit, dat[!dat$correct,], mean)
ggplot(df, aes(x = as.factor(condit), y = rt, color = as.factor(pos))) + 
  geom_point(size = 3) + 
  ggtitle(label = "error RTs per condit") + 
  theme_bw()

df <- aggregate(correct ~ pos*condit, dat, mean)
ggplot(df, aes(x = as.factor(condit), y = correct, color = as.factor(pos))) + 
  geom_point(size = 3) + 
  ggtitle(label = "%correct per condit") + 
  theme_bw()

df <- aggregate(rt ~ uniq*condit, dat[dat$correct,], mean)
ggplot(df, aes(x = as.factor(condit), y = rt, color = as.factor(uniq))) + 
  geom_point(size = 3) + 
  ggtitle(label = "correct RTs per condit") + 
  theme_bw()

df <- aggregate(rt ~ uniq*condit, dat[!dat$correct,], mean)
ggplot(df, aes(x = as.factor(condit), y = rt, color = as.factor(uniq))) + 
  geom_point(size = 3) + 
  ggtitle(label = "error RTs per condit") + 
  theme_bw()

df <- aggregate(correct ~ uniq*condit, dat, mean)
ggplot(df, aes(x = as.factor(condit), y = correct, color = as.factor(uniq))) + 
  geom_point(size = 3) + 
  ggtitle(label = "%correct per condit") + 
  theme_bw()