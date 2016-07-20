# -------------------------------------------------------------------------
# Rainfall and roof shape
# Mark Dingemanse 2016
# -------------------------------------------------------------------------

# Objectives --------------------------------------------------------------

# 1. Play with the new D-PLACE database 

# 2. Test a hunch about the relation of roof shape to rainfall (briefly, lots of
# rain negatively selects against flat roofs). Somewhat trivial, but easy to 
# check, lots of data, and fieldwork experience in Ghana provides relevant 
# ethnographic and historical background (including a possible account for an
# apparent exception).

# 3. Do this while controlling for language family and distance


# Preliminaries -----------------------------------------------------------

# Clear workspace
rm(list=ls())

# check for /in/ and /out/ directories (create them if needed)
add_working_dir <- function(x) { if(file.exists(x)) { cat(x,"dir:",paste0(getwd(),"/",x,"/")) } else { dir.create(paste0(getwd(),"/",x)) 
  cat("subdirectory",x,"created in",getwd()) } }
add_working_dir("in")
add_working_dir("out")

# Packages and useful functions
list.of.packages <- c("ggplot2","ggthemes","dplyr","lme4","extremevalues","mapproj","sp","rgeos")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
lapply(list.of.packages, require, character.only=T)
rm(list.of.packages,new.packages)

`%notin%` <- function(x,y) !(x %in% y) 

# "Research that uses data from D-PLACE should cite both the original source(s)
# of  the data and the paper by Kirby et al. in which D-PLACE was first
# presented  (e.g., research using cultural data from the Binford
# Hunter-Gatherer dataset: ""Binford (2001);  Binford and Johnson (2006); Kirby
# et al. 2016)."" The reference list should include the date data were 
# accessed and URL for D-PLACE (https://d-place.org), in addition to the full
# references for Binford (2001),  Binford and Johnson (2006), and Kirby et al.
# (2016)."



# Hypothesis -----------------------------------------------------------

# Main hunch: Roof type is related to precipitation such that round/sloping 
# roofs more likely when there is much rain. Note the directionality: we expect 
# few flat roofs in high rainfall societies, but we have no predictions for what
# happens in low rainfall societies. In causal terms, high rainfall selects 
# *against* flat roofs, but low rainfall does not select *for* them.

# Data -----------------------------------------------------------------


# Here's a direct link to the query on D-PLACE
# https://d-place.org/societies?c=%5B%7B%22variable%22:83,%22id%22:654%7D,%7B%22variable%22:83,%22id%22:655%7D,%7B%22variable%22:83,%22id%22:656%7D,%7B%22variable%22:83,%22id%22:657%7D,%7B%22variable%22:83,%22id%22:658%7D,%7B%22variable%22:83,%22id%22:659%7D,%7B%22variable%22:83,%22id%22:660%7D,%7B%22variable%22:83,%22id%22:661%7D,%7B%22variable%22:83,%22id%22:662%7D,%7B%22variable%22:88,%22id%22:700%7D,%7B%22variable%22:88,%22id%22:701%7D,%7B%22variable%22:88,%22id%22:702%7D,%7B%22variable%22:88,%22id%22:703%7D,%7B%22variable%22:88,%22id%22:704%7D,%7B%22variable%22:88,%22id%22:705%7D,%7B%22variable%22:88,%22id%22:706%7D,%7B%22variable%22:88,%22id%22:707%7D,%7B%22variable%22:88,%22id%22:708%7D,%7B%22variable%22:88,%22id%22:699%7D%5D&e=%5B%5B3,%22inrange%22,%5B%220.9970%22,%22547.9482%22%5D%5D%5D

# Load data (skipping line 1, which is a note about attribution)

d <- read.csv(file="in/dplace-export.csv",skip=1)

# Rename unwieldy variable names
names(d)[names(d)=="Description..EA082.House.construction..shape.of.roof"] <- "roofdesc"
names(d)[names(d)=="Code..EA082.House.construction..shape.of.roof"] <- "roofcode"
names(d)[names(d)=="Variable..Monthly.Mean.Precipitation..mm."] <- "rainfall"
names(d)[names(d)=="Code..EA087.House.construction..secondary.house.type...shape.of.roof"] <- "roofcode.secondary"
names(d)[names(d)=="Revised.longitude"] <- "longitude"
names(d)[names(d)=="Revised.latitude"] <- "latitude"

# What's the mapping exactly?
code <- unique(d$roofcode)
desc <- unique(d$roofdesc)
roofs <- data.frame(code,desc)

# make summary variable flatroof, 1 for flat, 0 for non-flat
# we have only 95 societies with flat roofs and 1035 with non-flat
d$flatroof <- ifelse(d$roofcode == 7,1,0) # i.e. flatroof = 1 if roofcode 7, 0 otherwise
d$flatroof <- as.factor(d$flatroof)
xtabs(~ flatroof,data=d)

# there is also a bit of data for roof shape of secondary buildings, though I 
# don't think we'll use it in this first analysis as it only covers 30% of the
# dataset.
d$flatroof.2nd <- ifelse(d$roofcode.secondary == 7,1,0) # i.e. flatroof = 1 if roofcode 7, 0 otherwise
d$flatroof.2nd <- as.factor(d$flatroof.2nd)
xtabs(~ flatroof.2nd,data=d)

xtabs(data=d, ~ flatroof + flatroof.2nd)

# So now: how does roof flatness relate to precipitation?

# Comparing means, we see that the average monthly rainfall for flat roof 
# societies is 48.1 (582 annually), while for non-flat roof societies it is
# twice as high at 105.8 (1270 annually). To put this in context, average annual
# rainfall over land on earth is 59.5 monthly (715 annually).
d %>%
  group_by(flatroof) %>%
  summarise(mean=mean(rainfall))

# What are the outliers? getOutliers() computes extreme values against a
# baseline expectation of a normal distribution

d.outliers <- getOutliers(d$rainfall)
outlierPlot(d$rainfall, d.outliers) 
# That QQ plot based on a normal distribution doesn't look optimal: two extreme 
# values are classified as outliers, together with the tail end of what looks like a
# continuous distribution. Can we do better?

# Oh wait, precipitation is probably better modelled using a Weibull 
# distribution: after all, 0 is a hard minimum whereas there is more room for 
# variation on the other end (cf. Wilks 1989). 
d.outliers.weibull <- getOutliers(d$rainfall,distribution="weibull")
outlierPlot(d$rainfall, d.outliers.weibull, title="Rainfall: QQ plot with outliers\nWeibull distribution, R² = 0.9988") 
# Whoah, rarely does one see such a beautiful fit — R² = 0.9988!

# Side note: we can also do this in ggplot2 of course
ggplot(d,aes(sample=rainfall)) +
  stat_qq(distribution=qweibull,dparams = list(shape=1.42, scale=110))
# And colouring the points already hints at different distributions for roof shape
ggplot(d,aes(sample=rainfall,colour=flatroof)) +
  stat_qq(distribution=qweibull,dparams = list(shape=1.8, scale=130))


# So what are those outliers?
d.outliers.weibull$nOut # both are on the right
d[c(d.outliers.weibull$iRight),c("Preferred.society.name","rainfall","latitude","longitude")]
# Monthly rainfall for both is recorded as 547.9482, probably taken or 
# interpolated from the same weather station(s). The lat/long places both on 
# Saipan Island, for which the annual rainfall is highly variable, but ranges 
# from 1900mm to 2300mm, i.e. 158-192 monthly, according to a 2004 report [1]. 
# It is hard to find sources with historical precipitation data for the focal 
# years of the Ethnographic Atlas (1930-1950), but I consider this enough reason
# to throw them out.

# [1] http://www.weriguam.org/reports/item/rainfall-climatology-for-saipan-distribution-return-periods-el-nino-tropical-cyclones-and-long-term-variations.html


d.no <- d[-d.outliers.weibull$iRight,] # we create a copy of the data with those societies excluded

d.no.outliers <- getOutliers(d.no$rainfall,distribution="weibull")
d.no.outliers$nOut # now there are no outliers anymore
outlierPlot(d.no$rainfall, d.no.outliers,title="Rainfall: QQ plot with outliers\nWeibull distribution, R² = 0.9989") 


# So what are the 3 societies with the highest rainfall in the flat roof 
# category? We take d.no, filter for flat roof societies, arrange to sort by
# rainfall, then slice to take the first few rows
top3 <- d.no %>%
  filter(flatroof == 1) %>%
  arrange(desc(rainfall)) %>%
  slice(1:3)




# Maps and plots
# --------------

# Mapping

m <- ggplot(d.no,aes(longitude,latitude)) + 
  theme_map() +
  borders("world", colour=NA, fill="#e7e8ea") +
  geom_point(data=filter(d.no,flatroof==0),aes(fill=rainfall,colour=flatroof),size=3,shape=21,stroke=1,alpha=0.7) +
  geom_point(data=filter(d.no,flatroof==1),aes(fill=rainfall,colour=flatroof),size=3,shape=22,stroke=1,alpha=0.7) +
  scale_colour_manual(values = c(NA,"brown"), guide = FALSE) +
  scale_fill_gradient("rainfall", low="#5edcff", high="#035280")
plot(m)
ggsave(file="out/map-rainfall-by-rooftype.png", width=12,height=8)

# map highest
m + geom_point(data=top3,aes(longitude,latitude),fill="red",colour=NA,size=3,shape=22,stroke=1)
ggsave(file="out/map-rainfall-by-rooftype-highest.png", width=12,height=8)

# zoom in on West-Africa
m + theme_bw() # use to target desired lat/long for zooming

zoom1 <- m + borders(colour="white") + coord_map(xlim = c(-18,7), ylim = c(20,2))

zoom1 + geom_point(data=top3,aes(longitude,latitude),fill="red",colour=NA,size=3,shape=22,stroke=1)
ggsave(file="out/map-rainfall-by-rooftype-west-africa.png", width=12,height=8)


# Plotting
myseed=101
set.seed(myseed)
p <- ggplot(d.no, aes(x=flatroof,y=rainfall,fill=rainfall,na.rm=T)) +
  theme_bw() +
  geom_jitter(alpha=0.7,shape=21,width=0.5,colour=NA,size=2,stroke=1) +
  scale_fill_gradient("rainfall", low="#5edcff", high="#035280") +
  stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size=0.3,width=0.33, geom = "crossbar")
plot(p)
ggsave(file="out/rainfall-by-rooftype.png", width=3,height=5)

# Highlight highest 5
d.no$top3 <- ifelse(d.no$Preferred.society.name %in% top3$Preferred.society.name, 1, NA)
d.no$top3 <- as.factor(d.no$top3)

set.seed(myseed)
pm <- ggplot(d.no, aes(x=flatroof,y=rainfall,fill=rainfall,colour=top3)) +
  theme_bw() +
  geom_jitter(alpha=0.7,width=0.5,shape=21,size=2,stroke=1) +
  scale_fill_gradient("rainfall", low="#5edcff", high="#035280") +
  scale_colour_manual(values=c("red",NA), guide = FALSE) +
  stat_summary(aes(x=flatroof,y=rainfall,fill=rainfall), inherit.aes = FALSE, fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", size=0.3,width=0.33, geom = "crossbar")
plot(pm)
ggsave(file="out/rainfall-by-rooftype-top3.png", width=3,height=5)


# Linear modelling --------------------------------------------------------------


# A simple linear model doesn't explain a lot of variation (R² 0.05) but comes 
# out as highly significant: if rain would not have any effect on roof type, the
# probability of this data is very low (F(1,1126)=59.5, p < 0.000001).

d.lm <- d.no
d.lm$flatroof <- as.integer(d.lm$flatroof)

model1 <- lm(flatroof ~ rainfall, d.lm)
summary(model1)
hist(residuals(model1))
qqnorm(residuals(model1))

# But not all observations are independent: we want to control for language 
# family, so we add that as a random intercept (meaning we allow the value of 
# roof to be 'dependent' to some degree on language family). We compare a model 
# with just language family with one where we also add rainfall as a fixed
# effect. The latter explains the spread of data significantly better (χ2(1) =
# 45.91, p < 0.00001, log likelihood difference 22.95).

model2 <- lmer(flatroof ~ (1|Language.family), data=d.lm)
summary(model2)
model3 <- lmer(flatroof ~ rainfall + (1|Language.family), data=d.lm)
summary(model3)

comparison <- anova(model2, model3)
comparison[1,4]-comparison[2,4] # log likelihood difference

# But the maps show that some of the flat roofs appear in 'clumps'. Perhaps 
# horizontal diffusion due to cultural contact can explain part of it. We need 
# to take into account distance to nearest neighbour. To do this we use the sp
# and rgeos packages (https://is.gd/FFLbKH).

# We create a spatial object
d.s <- d.no
coordinates(d.s) <- ~longitude+latitude # convert coordinates
distances <- gDistance(d.s, byid=T) # compute distances
min.dist <- apply(distances,1, function(x) order(x, decreasing =F)[2]) # find second shortest distance (the shortest is to itself)
d.no$distance <- min.dist # add distances to df

d.lm <- d.no # update the copy of the data used for linear modelling
d.lm$flatroof <- as.integer(d.lm$flatroof)

# We create a model where we try to predict flat roof by rainfall and distance. 
# We keep language family as a random effect. Essentially, we check whether 
# rainfall is still important when we control for historical contingency 
# (language family as a random effect) as well as direct social diffusion 
# (distance as a fixed effect). Turns out distance doesn't have a very big
# effect (t value 2.25, model estimate shows it's two orders of magnitude
# smaller than rainfall).
model4 <- lmer(flatroof ~ rainfall + distance + (1|Language.family),data=d.lm)
summary(model4)

# Comparing model3 to model4, the latter may explain the spread of data
# slightly better (χ²(1) = 5.086, p = 0.02412, log likelihood difference
# 2.54), but only if we use a loose criterion of p < 0.05

comparison <- anova(model3, model4)
comparison[1,4]-comparison[2,4] # log likelihood difference

# Conclusion: in this dataset, rainfall is the best predictor of roof shape, at 
# least in the sense that flat roofs are unlikely to occur in locales with a lot
# of precipitation.