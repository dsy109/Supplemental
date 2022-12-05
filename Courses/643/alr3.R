#  pureErrorAnova function from alr3 package, which is now archived

pureErrorAnova <- function(mod){UseMethod("pureErrorAnova")}
pureErrorAnova.lm <- function(mod) {
 if (is.null(mod$model)) mod <- update(mod, model=TRUE)
 p <- dim(mod$model)[2] -1
 mod$model$Lack.of.Fit <-
   factor(randomLinComb(model.matrix(mod), 101319853))
 aov1 <- anova(mod)
 #set.seed(save.seed) # restore random number seed
 if (length(levels(mod$model$Lack.of.Fit)) == length(mod$model$Lack.of.Fit))
  aov1 else {
  aov2 <- anova(lm(mod$model[ , 1]~mod$model$Lack.of.Fit, weights=weights(mod)))
  rrow <- dim(aov1)[1]
  aov2[1, 1] <- aov1[rrow, 1]-aov2[2, 1]
  aov2[1, 2] <- aov1[rrow, 2]-aov2[2, 2]
  aov2[1, 3] <- aov2[1, 2]/aov2[1, 1]
  aov1[1:(rrow-1), 4] <- aov1[1:(rrow-1), 3]/aov2[2, 3]
  aov2[1, 4] <- aov2[1, 3]/aov2[2, 3]
  row.names(aov2) <- c(" Lack of fit", " Pure Error")
  aov <- rbind(aov1, aov2)
  aov[ , 5] <- pf(aov[ , 4], aov[ , 1], aov2[2, 1], lower.tail=FALSE)
  aov
  }}


randomLinComb <- function(X, seed=NULL) {UseMethod("randomLinComb")}

randomLinComb.default <- function(X, seed=NULL) {
 if(!is.null(seed)) set.seed(seed)
 std <- function(x){
    s <- sd(x)
    if( s > 0) (x-mean(x))/s else x}
 as.vector(apply(X, 2, std)%*% as.vector(2*rnorm(dim(X)[2])-1) )
 }

randomLinComb.lm <- function(X, ...) {
  randomLinComb(model.matrix(X), ...)}

randomLinComb.lm <- function(X, seed=NULL) {
 if(is.null(X$model)) X <- update(X, model=TRUE)
 randomLinComb(X$model[ , -1], seed=seed)}

