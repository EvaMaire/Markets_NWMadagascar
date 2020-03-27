################################################################################################################              
#                                                                                                              #                          
# Code by Eva Maire (emg.maire@gmail.com). Last update on 2020-03-27                                           #
# This code provides results of analysis used in Maire, E., D'agata, S., Aliaume, C., Mouillot, D.,            #
# Darling, D., Ramahery, V., Ranaivoson, R., Randriamanantsoa, B., Tianarisoa, T., Santisy, A., &              #
# Cinner J. (2020). Disentangling the complex roles of markets on coral reefs in northwest Madagascar.         #
# Ecology and Society.                                                                                         #
#                                                                                                              #
################################################################################################################

#loading required libraries
require(visreg)
require(mgcv)
require(ggplot2)
require(dplyr)
require(MuMIn)
require(FactoMineR)
require(factoextra)
require(gridExtra)

# setting working directory
my_path<-"" #  <= folder with RData
setwd(my_path)

#############################################################################################################
# Figure 2 - Partial effects of each socioeconomic covariate predicting log fish biomass in the model 
# while considering the other predictor variables are held constant. Relationships between fish biomass and 
# travel time from the nearest community (a), 
# travel time from the nearest market (b), 
# human population size (c) and 
# management (d) for reefs where fishing is permitted (orange) and prohibited (green).

#Load fish biomass data
load("fish_biomass.RData")

# Step 1 - we first performe a Principal Components Analysis (PCA) using a set of six reef habitat and environmental variables 
# (depth, weekly average SST and primary productivity, reef complexity, percent cover of macroalgae and live hard coral) to 
# describe similarities between our ecological sites while dealing with multicollinearity.

# Six reef habitat and environmental variables
Pred <- c("Depth","Live_Hard_Coral","Macroalgae","Complexity","MeanChl","MeanTemp")
Ncol=vector(length=length(Pred))

for (i in 1:length(Pred))
{
  Ncol[i]=which(colnames(fish_biomass)==Pred[i])
}
Ncol

env=fish_biomass[,Ncol]
head(env)

PCA(env, scale.unit = TRUE, graph = F)
res.pca <- PCA(env, graph = FALSE)

# Appendix 1 - Associations between environmental and benthic conditions of reefs through a Principal Component 
# Analysis and corresponding loadings. 

pca_S1 <- fviz_pca_var(res.pca, col.var = "contrib",
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),midpoint=15,
                  repel = T,title="") + theme_minimal()

pca_S1

#Loadings
loadings <- sweep(res.pca$var$coord,2,sqrt(res.pca$eig[1:ncol(res.pca$var$coord),1]),FUN="/")
loadings_S1 <- tableGrob(round(loadings,2),theme = ttheme_minimal(base_size = 15))

# 2-panel Plot
grid.arrange(pca_S1,loadings_S1) # Figure_S1 

# We only retaine the first two components (representing 56% of the total variance, see Figure S1) 
# that mix abiotic and benthic conditions as environmental covariates for further analysis.  
pc1 <- res.pca$ind$coord[,1]
pc2 <- res.pca$ind$coord[,2]

fish_biomass$pc1 <- pc1
fish_biomass$pc2 <- pc2

# To explore how proximity to markets and communities affect reef fish biomass beyond ecological and human population size effects, 
# we build Generalized Additive Models (GAMs) considering the two environmental covariates provided by PCA (pc1 et pc2), 
# human population size, travel time from human settlements and markets, and management. 

mod<-gam(log_fish_biomass~s(pc1,bs="cr",k=3)+s(pc2,bs="cr",k=3)+s(human_pop,bs="cr",k=3)+
          management+s(tt_village,bs="cr",k=3)+s(tt_market,bs="cr",k=3)
        +s(tt_village, bs = "cr", by = management, k = 3, m = 1)
        ,data=fish_biomass,na.action="na.fail") 

AICc(mod) #AICc = 14.9
summary(mod) #R-sq (adj) = 0.83

# Step 2 - Each partial effect plot needs to be limited to the range of the data

#Travel time from the nearest village
visreg_village <- visreg(mod,"tt_village",by="management",overlay=T,ylim=c(2,4),band=T,plot=F) 
pred <- mod$fitted.values
ttv <- fish_biomass$tt_village
management <- fish_biomass$management

nb_perm <- length(which(fish_biomass$management=="fishing_permitted"))
nb_pro <- length(which(fish_biomass$management=="fishing_prohibited"))
nafperm <- rep(NA,length(which(visreg_village$fit$management=="fishing_permitted"))-nb_perm)
nafpro <- rep(NA,length(which(visreg_village$fit$management=="fishing_prohibited"))-nb_pro) 

pred_bm_fperm <- c(pred[which(management=="fishing_permitted")],nafperm)
pred_bm_fpro <- c(pred[which(management=="fishing_prohibited")],nafpro)
ttv_fperm <- c(ttv[which(management=="fishing_permitted")],nafperm)
ttv_fpro <- c(ttv[which(management=="fishing_prohibited")],nafpro)

datpro <- subset(visreg_village$fit, visreg_village$fit$management == "fishing_prohibited")
datpro <- datpro[,c("tt_village","visregFit","visregLwr","visregUpr")]
colnames(datpro) <- c("tt_village","fitpro","lowpro","highpro")
datperm <- subset(visreg_village$fit, visreg_village$fit$management == "fishing_permitted")
datperm <- datperm[,c("visregFit","visregLwr","visregUpr")]
colnames(datperm) <- c("fitperm","lowperm","highperm")

datvillage <- cbind(datpro,datperm,pred_bm_fperm,ttv_fperm,pred_bm_fpro,ttv_fpro) 
delv <- which(datvillage$tt_village > max(datvillage$ttv_fperm,na.rm=T) )
datvillage[delv,c("highperm","lowperm","fitperm")]<-NA

#Plot
white_theme<-theme(axis.text=element_text(colour="black",size=16),
                   axis.title=element_text(size=18),
                   axis.ticks=element_line(colour="black"),
                   panel.grid.minor=element_blank(),
                   panel.background=element_rect(fill="white",colour="black"),
                   legend.justification=c(1,0),legend.position=c(.95, .05),line= element_blank(),
                   plot.background=element_rect(fill="transparent",colour=NA))

pvillage <- ggplot(datvillage)+
  geom_line(aes(tt_village,fitpro),colour="#5ab4ac",size=1)+
  geom_ribbon(aes(x=tt_village,ymax=highpro,ymin=lowpro),fill="#5ab4ac",alpha=0.5)+
  geom_line(aes(tt_village,fitperm),colour="#d8b365",size=1)+
  geom_ribbon(aes(x=tt_village,ymax=highperm,ymin=lowperm),fill="#d8b365",alpha=0.5)+
  scale_x_continuous("Travel time from community (h)")+
  scale_y_continuous("Log fish biomass (kg/ha)",limits=c(1.72,3.66),breaks=c(1.88,2.2,2.48,2.7,2.9,3,3.18,3.3,3.48,3.6),
                     labels=c(75,150,300,500,800,1000,1500,2000,3000,4000))+
  white_theme+
  annotate("text", x = 0, y=3.66, label = "a", size=7,fontface =2)+
  #legend
  annotate("text", x = 3.29, y=1.75, label = "Fishing permitted",size=5)+
  annotate("text", x = 3.31, y=1.85, label = "Fishing prohibited",size=5)+ 
  geom_segment(aes(x = 2, y = 1.75, xend = 2.4, yend = 1.75), size=2, col = "#d8b365")+
  geom_segment(aes(x = 2, y = 1.85, xend = 2.4 , yend = 1.85), size=2, col = "#5ab4ac")


#Human population size 
visreg_pop <- visreg(mod,"human_pop",by="management",overlay=T,ylim=c(1,4),band=T,plot=F) 
pop <- fish_biomass$human_pop
pop_fperm <- c(pop[which(management=="fishing_permitted")],nafperm)
pop_fpro <- c(pop[which(management=="fishing_prohibited")],nafpro)

datpro <- subset(visreg_pop$fit, visreg_pop$fit$management == "fishing_prohibited")
datpro <- datpro[,c("human_pop","visregFit","visregLwr","visregUpr")]
colnames(datpro) <- c("human_pop","fitpro","lowpro","highpro")
datperm <- subset(visreg_pop$fit, visreg_pop$fit$management == "fishing_permitted")
datperm <- datperm[,c("visregFit","visregLwr","visregUpr")]
colnames(datperm) <- c("fitperm","lowperm","highperm")

datpop <- cbind(datpro,datperm,pred_bm_fperm,pop_fperm,pred_bm_fpro,pop_fpro) 
del <- which(datpop$human_pop > max(datpop$pop_fpro,na.rm=T) )
datpop[del,c("highpro","lowpro","fitpro")]<-NA

#Plot
ppop <- ggplot(datpop)+
  geom_line(aes(human_pop,fitpro),colour="#5ab4ac",size=1)+
  geom_ribbon(aes(x=human_pop,ymax=highpro,ymin=lowpro),fill="#5ab4ac",alpha=0.5)+
  geom_line(aes(human_pop,fitperm),colour="#d8b365",size=1)+
  geom_ribbon(aes(x=human_pop,ymax=highperm,ymin=lowperm),fill="#d8b365",alpha=0.5)+
  scale_x_continuous("Log human population size",limits=c(0,4))+
  scale_y_continuous("Log fish biomass (kg/ha)",limits=c(1.72,3.66),breaks=c(1.88,2.2,2.48,2.7,2.9,3,3.18,3.3,3.48,3.6),
                     labels=c(75,150,300,500,800,1000,1500,2000,3000,4000))+
  white_theme+
  annotate("text", x = 0, y=3.66, label = "c", size=7,fontface =2)


#Travel time from the nearest market
visreg_market <- visreg(mod,"tt_market",by="management",overlay=T,ylim=c(1,4),band=T,plot=F) 
ttm <- fish_biomass$tt_market
ttm_fperm <- c(ttm[which(management=="fishing_permitted")],nafperm)
ttm_fpro <- c(ttm[which(management=="fishing_prohibited")],nafpro)

datpro <- subset(visreg_market$fit, visreg_market$fit$management == "fishing_prohibited")
datpro <- datpro[,c("tt_market","visregFit","visregLwr","visregUpr")]
colnames(datpro) <- c("tt_market","fitpro","lowpro","highpro")
datperm <- subset(visreg_market$fit, visreg_market$fit$management == "fishing_permitted")
datperm <- datperm[,c("visregFit","visregLwr","visregUpr")]
colnames(datperm) <- c("fitperm","lowperm","highperm")

datmarket <- cbind(datpro,datperm,pred_bm_fperm,ttm_fperm,pred_bm_fpro,ttm_fpro) 
delm <- which(datmarket$tt_market > max(datmarket$ttm_fperm,na.rm=T) )
datmarket[delm,c("highperm","lowperm","fitperm")]<-NA

# Plot
pmarket <- ggplot(datmarket)+
  geom_line(aes(tt_market,fitpro),colour="#5ab4ac",size=1)+
  geom_ribbon(aes(x=tt_market,ymax=highpro,ymin=lowpro),fill="#5ab4ac",alpha=0.5)+
  geom_line(aes(tt_market,fitperm),colour="#d8b365",size=1)+
  geom_ribbon(aes(x=tt_market,ymax=highperm,ymin=lowperm),fill="#d8b365",alpha=0.5)+
  scale_x_continuous("Travel time from market (h)",limits=c(1,10),breaks=c(1,2,3,4,5,6,7,8,9,10),
                     labels=c(1,2,3,4,5,6,7,8,9,10))+
  scale_y_continuous("Log fish biomass (kg/ha)",limits=c(1.72,3.66),breaks=c(1.88,2.2,2.48,2.7,2.9,3,3.18,3.3,3.48,3.6),
                     labels=c(75,150,300,500,800,1000,1500,2000,3000,4000))+
  white_theme+
  annotate("text", x = 1, y=3.66, label = "b", size=7,fontface =2)


#Management 
visreg_management <- visreg(mod,"management",ylim=c(1,4),band=T,plot=F) 

ttm_fperm <- c(ttm[which(management=="fishing_permitted")],nafperm)
ttm_fpro <- c(ttm[which(management=="fishing_prohibited")],nafpro)

datpro <- subset(visreg_management$fit, visreg_management$fit$management == "fishing_prohibited")
datpro <- datpro[,c("visregFit","visregLwr","visregUpr")]
colnames(datpro) <- c("fitpro","lowpro","highpro")
datperm <- subset(visreg_management$fit, visreg_management$fit$management == "fishing_permitted")
datperm <- datperm[,c("visregFit","visregLwr","visregUpr")]
colnames(datperm) <- c("fitperm","lowperm","highperm")

x1 <- c(1,2.5) 
y1 <- c(datperm$lowperm,datpro$lowpro)
x2 <- c(2,3.5) 
y2 <- c(datperm$highperm,datpro$highpro)
type=c("Fishing permitted","Fishing prohibited")

d <- data.frame(y1,y2,x1,x2,type)

#Plot
white_theme2<-theme(axis.text=element_text(colour="black",size=16),
                    axis.title=element_text(size=18),
                    axis.ticks=element_line(colour="black"),
                    panel.grid.minor=element_blank(),
                    panel.background=element_rect(fill="white",colour="black"),
                    legend.position='none',
                    plot.background=element_rect(fill="transparent",colour=NA))

pmanagement <- ggplot()+
  geom_rect(data = d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2,fill=type),  alpha=0.5) +
  scale_fill_manual(values = alpha(c("Fishing permitted"= "#d8b365","Fishing prohibited" = "#5ab4ac"), 0.5))+
  geom_segment(aes(x = 1, y = fitperm, xend = 2, yend = fitperm, colour = "segment"),col="#d8b365",size=2, data = datperm)+
  geom_segment(aes(x = 2.5, y = fitpro, xend = 3.5, yend = fitpro, colour = "segment"),col="#5ab4ac",size=2, data = datpro)+
  scale_x_continuous("Management",limits=c(1,3.5),breaks=c(1.5,3),
                     labels=c("Fishing permitted","Fishing prohibited"))+
  scale_y_continuous("Log fish biomass (kg/ha)",limits=c(1.72,3.66),breaks=c(1.88,2.2,2.48,2.7,2.9,3,3.18,3.3,3.48,3.6),
                     labels=c(75,150,300,500,800,1000,1500,2000,3000,4000))+
  white_theme2+
  annotate("text", x = 1, y=3.66, label = "d", size=7,fontface =2)


#Multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} # end of the function

# 4-panel plot
#jpeg("Partial_effects_Maire_et_al_2020.jpeg",res=300, width=3000, height=3000)
multiplot(pvillage, ppop, pmarket, pmanagement,  cols=2)
#graphics.off()

rm(list=ls())
#############################################################################################################



#############################################################################################################
# Figure 3 - Associations between market proximity and (a) socioeconomic characteristics of communities 
# through a Principal Components Analysis and (b) selling strategies.

# (a) 
load("coastal_communities.RData")

#PCA 
res.pca <- PCA(coastal_communities[,-13], graph = F, quanti.sup=12) #Management as supplementary variable and remove the identity of Market for PCA

colvar <- c("#66c2a5","#66c2a5",
            "#fc8d62","#fc8d62","#fc8d62","#fc8d62",
            "#8da0cb","#8da0cb","#8da0cb","#8da0cb",
            "black")

var <- c("aComposition","aComposition",
         "cScale","cScale","cScale","cScale",
         "bTechnique","bTechnique","bTechnique","bTechnique",
         "dMarket access")

themepca <-theme(legend.title = element_text(size=15), 
                 legend.text=element_text(size=14),
                 axis.text.x = element_text(size=14),
                 axis.text.y = element_text(size=14),
                 axis.title.x =element_text(size=14),
                 axis.title.y =element_text(size=14))

p <- fviz_pca_biplot(res.pca, 
                     #Indivdiduals
                     pointshape = 21, pointsize=3, repel=T, col.ind = "black", fill.ind = coastal_communities$Market, mean.point = F,
                     #Variables
                     col.var = var, addEllipses=F, ellipse.level=0.90,fill.var = var, col.quanti.sup = "black",arrowsize=1,
                     labelsize=4,title="", legend.title = list(fill = "Nearest market", color = "Effects") ) + themepca 

pca <- p + scale_color_manual(labels = c("Composition","Technique","Scale","Market access"), values = c("#1b9e77", "#d95f02","#7570b3","#e7298a"))+
           scale_fill_manual(labels=c("Ambanja","Ambilobe","Hell Ville"),values=c("#f7f7f7","#252525","#969696"))

# (b)
load("selling_strategies.RData")
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

white_theme<-theme(axis.text=element_text(colour="black",size=14),
                   axis.title=element_text(size=15),
                   axis.ticks=element_line(colour="black"),
                   panel.grid.minor=element_blank(),
                   panel.background=element_rect(fill="white",colour="black"),
                   legend.justification=c(1,0),legend.position=c(.95, .05),line= element_blank(),
                   plot.background=element_rect(fill="transparent",colour=NA))

midd <- arrange(selling_strategies, factor(selling, levels = c("Middlemen & market","Middlemen","Community only")))

middlemen <- ggplot(midd, aes(x=selling, y=travel_time_market)) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill="#8C9091")+
  scale_y_continuous("Travel time from market (h)")+
  scale_x_discrete("") +
  coord_flip() +
  white_theme

# 2-panel plot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} # end of the function

multiplot(pca, middlemen, cols=2)

#############################################################################################################

# Appendices

#############################################################################################################
# Appendix 2 - Correlogram showing correlations between the ten socioeconomic indicators measured of coastal 
# communities and market access. 

require(corrgram)
require(corrplot)

ccor <- coastal_communities[,-c(12,13)] #Delete Management and Market
cor_mat <- cor(ccor)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(ccor)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(cor_mat, method="color", col=col(200),  
         type="upper",  
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 0.9, number.cex = .7, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 1, insig = "blank",  
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

#############################################################################################################
# Appendix 3 - Scores (cos2) of (a) each variable and (b) community integrated in the PCA linking 
# market access and social characteristics of coastal communities.
res.pca <- PCA(coastal_communities[,-13], graph = F, quanti.sup=12) #Management as supplementary variable and remove the identity of Market for PCA

var <- fviz_cos2(res.pca, choice = "var", labelsize=3, axes = 1:2,
               title="")+theme_minimal()+geom_hline(yintercept=0.4)

var2 <- ggpubr::ggpar(var,
                    title = "",
                    xlab = "Variables", ylab = "Cos2",
                    ggtheme = theme_minimal(), palette = "jco",
                    cex.lab=2,ylim=c(0,1),font.tickslab=c(10),font.xtickslab=8)

var3 <- var2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

comm <- fviz_cos2(res.pca, choice = "ind", labelsize=3, axes = 1:2,
               title="")+theme_minimal()+geom_hline(yintercept=0.4)

comm2 <- ggpubr::ggpar(comm,
                    title = "",
                    xlab = "Communities", ylab = "Cos2",
                    ggtheme = theme_minimal(),
                    palette = "jco",
                    cex.lab=2,ylim=c(0,1),font.tickslab=c(10),font.xtickslab=8)

comm3 <- comm2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

multiplot(var3, comm3, cols=2)

#############################################################################################################
# Appendix 4 - Associations between market proximity and selling fish catches
sold <- arrange(selling_strategies, factor(selling, levels = c("Community only","Middlemen & market","Middlemen")))

prop_fish_sold <- ggplot(aes(x=selling, y=fish_sold),data=sold) +
  stat_summary(fun.data = quantiles_95, geom="boxplot",fill="#8C9091")+
  scale_y_continuous(c("% of fish sold"),limits=c(65,95),breaks=c(65,70,75,80,85,90,95))+
  scale_x_discrete("") +
  coord_flip() +
  white_theme


plot(prop_fish_sold)

#############################################################################################################
#End of script
#
