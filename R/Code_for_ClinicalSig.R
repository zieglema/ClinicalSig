#setwd("/Users/Matthias/Seafile/Meine Bibliothek/R Projects/ClinicalSig")
#setwd("~/Meine Bibliothek/R/R Projects/ClinicalSig")

#Example data set####
#library(openxlsx)
#data <- read.xlsx("clsi.xlsx")
#dat <- data

#t test function for paired t test including Hedges g and power####

#' t test function for paired t test including Hedges g and power
#'
#' @param pre The variable containing the values from before the intervention
#' @param post The variable containing the values from after the intervention
#' @param data The data set containing the variables
#' @param alternative Direction of the alternative. Can have the values "less" or "greater".
#' @param sig.level The accepted probability for a type I error.
#'
#' @return a list containing results from t test for dependent variables, effect size and power
#' @export
#' @examples
#' tTest(pre=U1_GDS_G, post=U2_GDS_G, data=dat, alternative="less", sig.level=.05)
#' @author Matthias Ziegler

#'
tTest <- function(pre, post, data, alternative, sig.level){
  ttest <- t.test(eval(substitute(post), data),eval(substitute(pre), data),
                  alternative=alternative, paired=TRUE, conf.level=1-sig.level)
  sed <- sd(eval(substitute(pre),data)-eval(substitute(post),data), na.rm=T)
  diff <- ttest$estimate
  g <- diff/sed
  names(g) <- "Effect Size Hedge's g"
  power <- pwr::pwr.t.test(n=nrow(data), d=g, sig.level=sig.level,
                      power=NULL, type="paired", alternative="two.sided")
  list(ttest, g, power)
}
#testing the function
#tTest(pre=U1_GDS_G, post=U2_GDS_G, data=dat, alternative="less", sig.level=.05)

#Function for the reliable change index####

#' Function for the reliable change index
#'
#' @param pre The variable containing the values from before the intervention
#' @param post The variable containing the values from after the intervention
#' @param data Data set containing the variables
#' @param rtt Test-retest correlation for the used measure. Should be from the relevant population and with no intervention in between
#' @param sdNorm Standard deviation of the used measure in the norm group.
#'
#' @return Returns the RCI for each person in the data set
#' @export
#' @examples
#' dataRci <- rci(pre="U1_GDS_G", post="U2_GDS_G", data=dat, rtt=.83, sdNorm=3.32)
#' hist(dataRci$rci, main="Histogram of RCI", xlab="RCI", breaks=30, col="blue")
#'
rci <- function(pre, post, data, rtt, sdNorm){
  if ("rci" %in% names(data)) {warning("variable rci already exists in data frame and will be overwritten")}
  difference <- data[,post]-data[,pre]
  Serr <- sdNorm*sqrt(1-rtt)
  Sdiff <- sqrt((2*Serr^2))
  data$rci <- difference/Sdiff
  data
}


#Function for classifying the RCI####
#' Title
#'
#' @param data Data set
#' @param nameRci Name of the variable in data set representing the reliable change index
#' @param sigLevel The accepted probability for a type I error.
#'
#' @return RCI variable is classified as changed or not
#' @export
#' @examples
#' dataRciClass <- rciClass(data=dataRci, nameRci="rci", sigLevel=.05)
#' barplot(table(dataRciClass$rciClass))
#'
rciClass <- function(data, nameRci, sigLevel){
  if ("rciClass" %in% names(data)) {warning("variable rciClass already exists in data frame and will be overwritten")}
  data <- within(data, {
    rciClass <- NA
    rciClass[abs(rci)<abs(qnorm(sigLevel/2, lower.tail=T))] <- "no reliable change"
    rciClass[abs(rci)>abs(qnorm(sigLevel/2, lower.tail=T))] <- "reliable change"
  })
  data
}


#Function for categorizing clinical significance####

#' Function for categorizing clinical significance
#'
#' @param data Data set
#' @param pre The variable containing the values from before the intervention
#' @param post The variable containing the values from after the intervention
#' @param consistency Estimate for the internal consistency of the used measure
#' @param mPath Mean score in clinical norm group
#' @param sdPath Standard deviation in clinical norm group
#' @param sigLevel The accepted probability for a type I error.
#'
#' @return Each person is assigned to one of 7 classes
#' @export
#' @examples
#' dataFinal <- clinSig(data=dataRciClass, pre="U1_GDS_G", post="U2_GDS_G", consistency=.91,
#' mPath=14.78, sdPath=3.32, sigLevel=.05)
#' table(dataFinal$csClass)
#' barplot(table(dataFinal$csClass))
#'
clinSig <- function(data, pre, post, consistency, mPath, sdPath, sigLevel){
  #Defining the criterion
  healthy <- mPath-2*sdPath
  #Assessing the status before and after therapy
  Se <- abs(qnorm(sigLevel/2, lower.tail=T))*sdPath*sqrt(1-consistency)
  upperCI <- healthy+Se
  lowerCI <- healthy-Se
  max <- max(data[[pre]])
  preStatus <- cut(x=data[[pre]],breaks=c(0,upperCI, max),
                   labels=c("negative diagnosis", "positive diagnosis"))
  postStatus <- cut(x=data[[post]],breaks=c(0,upperCI, max),
                    labels=c("negative diagnosis", "positive diagnosis"))
  #Actual difference
  difference <- data[,post]-data[,pre]
  #Classifying according to criterion 1 from Jacobson & Truax
  data <- within(data, {
    csClass <- NA
    csClass[abs(rci)<=abs(qnorm(sigLevel/2, lower.tail=T))] <-"no reliable change"
    csClass[abs(rci)>abs(qnorm(sigLevel/2, lower.tail=T)) & preStatus=="positive diagnosis" & postStatus=="negative diagnosis" & difference<0] <-"Clinically significant improvement"
    csClass[abs(rci)>abs(qnorm(sigLevel/2, lower.tail=T)) & preStatus=="positive diagnosis" & postStatus=="positive diagnosis" & difference<0] <-"Clinically insignificant improvement"
    csClass[abs(rci)>abs(qnorm(sigLevel/2, lower.tail=T)) & preStatus=="positive diagnosis" & postStatus=="positive diagnosis" & difference>0] <-"Clinically insignificant deterioration"
    csClass[abs(rci)>abs(qnorm(sigLevel/2, lower.tail=T)) & preStatus=="negative diagnosis" & postStatus=="positive diagnosis" & difference>0] <-"Clinically significant deterioration"
    csClass[abs(rci)>abs(qnorm(sigLevel/2, lower.tail=T)) & preStatus=="negative diagnosis" & postStatus=="negative diagnosis" & difference>0] <-"Significant but clinically irrelevant deterioration"
    csClass[abs(rci)>abs(qnorm(sigLevel/2, lower.tail=T)) & preStatus=="negative diagnosis" & postStatus=="negative diagnosis" & difference<0] <-"Significant but clinically irrelevant improvement"
  })
  data
}


#Function for plotting rci and rci class####
#' Function for plotting rci and rci class
#'
#' @param data Data set
#' @param pre The variable containing the values from before the intervention
#' @param post The variable containing the values from after the intervention
#' @param sigLevel The accepted probability for a type I error.
#' @param rtt Test-retest correlation for the used measure. Should be from the relevant population and with no intervention in between
#' @param consistency Estimate for the internal consistency of the used measure.
#' @param mPath Mean score in clinical norm group
#' @param sdPath Standard deviation in clinical norm group
#' @param scaleName Name of the used instrument
#' @param min Theoretical score minimum in instrument used
#' @param max Theoretical score maximum in instrument used
#'
#' @return A plot depicting the 7 classes and all CIs
#' @export
#' @examples
#' plotCs(data=dat, pre="U1_GDS_G", post="U2_GDS_G", sigLevel=.05, consistency=.91, rtt=.83,
#' mPath=14.78, sdPath=3.32, scaleName="GDS Score", min=0, max=30)
#'
plotCs <- function(data, pre, post, sigLevel, rtt, consistency, mPath, sdPath, scaleName, min, max){
  #Estimating Standard error of difference
  difference <- data[,post]-data[,pre]
  Serr <- sdPath*sqrt(1-rtt)
  Sdiff <- sqrt((2*Serr^2))
  #Defining health status
  healthy <- mPath-2*sdPath
  #Defining the CIs
  Se <- abs(qnorm(sigLevel/2, lower.tail=TRUE))*sdPath*sqrt(1-consistency)
  upperCI <- healthy+Se
  lowerCI <- healthy-Se
  #plot
  plot(data[,pre]~data[,post], xlab=paste("Post", scaleName, sep=" "),
       ylab=paste("Pre", scaleName, sep=" "), ylim=c(min,max), xlim=c(min,max))
  #Adding the bisecting line
  abline(a=0, b=1, col="red")
  #Adding the CI
  abline(a=qnorm(sigLevel/2, lower.tail=TRUE)*Sdiff, b=1, col="red", lty="dashed")
  abline(a=abs(qnorm(sigLevel/2, lower.tail=TRUE))*Sdiff, b=1, col="red", lty="dashed")
  #Adding the cutoff for diagnosis
  abline(a=healthy, b=0, col="blue")
  abline(a=upperCI, b=0, col="blue", lty="dashed")
  abline(a=lowerCI, b=0, col="blue", lty="dashed")
  abline(v=healthy, col="green")
  abline(v=upperCI, col="green", lty="dashed")
  abline(v=lowerCI, col="green", lty="dashed")
}

