# last modified 21 December 2006 by J. Fox
# last modified 15 April 2009 by M. Friendly 
#  -- Fixed numerous warnings resulting from axes=FALSE
#  -- prepare to generalize diagonal panel
# last modified 2 Feb 2010 by M. Friendly
#  -- added code for repeated measures designs
# last modified 13 Oct 2011 by M. Friendly
#  -- added var.labels 




#' Pairwise HE Plots
#' 
#' The function (in the form of an \code{mlm} method for the generic
#' \code{\link[graphics]{pairs}} function) constructs a ``matrix'' of pairwise
#' HE plots (see \link{heplot}) for a multivariate linear model.
#' 
#' 
#' @param x an object of class \code{mlm}.
#' @param variables indices or names of the three of more response variables to
#'        be plotted; defaults to all of the responses.
#' @param var.labels labels for the variables plotted in the diagonal panels;
#'        defaults to names of the response variables.
#' @param var.cex character expansion for the variable labels.
#' @param type type of sum-of-squares-and-products matrices to compute; one
#'        of \code{"II"}, \code{"III"}, \code{"2"}, or \code{"3"}, where \code{"II"}
#'        is the default (and \code{"2"} is a synonym).
#' @param idata an optional data frame giving a factor or factors defining the
#'        intra-subject model for multivariate repeated-measures data.  See Details of
#' \code{\link[car]{Anova}} for an explanation of the intra-subject design and
#'        for further explanation of the other arguments relating to intra-subject factors.
#' @param idesign a one-sided model formula using the ``data'' in idata and
#'        specifying the intra-subject design for repeated measure models.
#' @param icontrasts names of contrast-generating functions to be applied by
#'        default to factors and ordered factors, respectively, in the within-subject
#'        ``data''; the contrasts must produce an intra-subject model matrix in which
#'        different terms are orthogonal. The default is c("contr.sum", "contr.poly").
#' @param imatrix In lieu of \code{idata} and \code{idesign}, you can specify
#'         the intra-subject design matrix directly via \code{imatrix}, in the form of
#'         list of named elements.  Each element gives the columns of the
#'         within-subject model matrix for an intra-subject term to be tested, and must
#'         have as many rows as there are responses; the columns of the within-subject
#'         model matrix for \emph{different} terms must be mutually orthogonal.
#' @param iterm For repeated measures designs, you must specify one
#'        intra-subject term (a character string) to select the SSPE (E) matrix used
#'        in the HE plot.  Hypothesis terms plotted include the \code{iterm} effect as
#'        well as all interactions of \code{iterm} with \code{terms}.
#' @param manova optional \code{Anova.mlm} object for the model; if absent a
#'        MANOVA is computed. Specifying the argument can therefore save computation
#'        in repeated calls.
#' @param offset.axes proportion to extend the axes in each direction; defaults to 0.05.
#' @param digits number of significant digits in axis end-labels; taken from
#'        the \code{"digits"} option.
#' @param fill A logical vector indicating whether each ellipse should be
#'        filled or not.  The first value is used for the error ellipse, the rest ---
#'        possibly recycled --- for the hypothesis ellipses; a single fill value can
#'        be given.  Defaults to FALSE for backward compatibility. See Details of
#' \code{\link{heplot}}
#' @param fill.alpha Alpha transparency for filled ellipses, a numeric scalar
#'        or vector of values within \code{[0,1]}, where 0 means fully transparent and
#'        1 means fully opaque. Defaults to 0.3.
#' @param \dots arguments to pass down to \code{heplot}, which is used to draw
#'        each panel of the display.
#' @author Michael Friendly
#' @seealso \code{\link{heplot}}, \code{\link{heplot3d}}
#' @references 
#' Friendly, M. (2006).  Data Ellipses, HE Plots and Reduced-Rank
#' Displays for Multivariate Linear Models: SAS Software and Examples
#' \emph{Journal of Statistical Software}, 17(6), 1-42.
#' \url{https://www.jstatsoft.org/v17/i06/}
#' 
#' Friendly, M. (2007).  HE plots for Multivariate General Linear Models.
#' \emph{Journal of Computational and Graphical Statistics}, 16(2) 421-444.
#' \url{http://datavis.ca/papers/jcgs-heplots.pdf}
#' @keywords hplot multivariate
#' @examples
#' 
#' # ANCOVA, assuming equal slopes
#' rohwer.mod <- lm(cbind(SAT, PPVT, Raven) ~ SES + n + s + ns + na + ss, data=Rohwer)
#' 
#' # View all pairs, with ellipse for all 5 regressors
#' pairs(rohwer.mod, hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss")))
#' 
#' 
#' @exportS3Method pairs mlm
#' @importFrom car Anova
#' @importFrom stats model.frame
pairs.mlm <-
  function(x, variables, var.labels, var.cex = 2,
           type=c("II", "III", "2", "3"),
           idata=NULL,
           idesign=NULL,
           icontrasts=NULL,
           imatrix=NULL,
           iterm=NULL,
           manova,        # an optional Anova.mlm object
           offset.axes=0.05, 
           digits=getOption("digits") - 1,
           fill=FALSE,         ## whether to draw filled ellipses (vectorized)
           fill.alpha=0.3,     ## alpha transparency for filled ellipses
           ...){
    
    #	manova <- Anova(x, type)
    if (missing(manova)) {
      type <- match.arg(type)
      if (is.null(imatrix)) {
        manova <- car::Anova(x, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts)
      }
      else {
        #			if (packageDescription("car")[["Version"]] >= 2)
        manova <- car::Anova(x, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix)
        #			else stop("imatrix argument requires car 2.0-0 or later")
      } 
    }   
    
    data <- model.frame(x)
    #	Y <- model.response(model.frame(x))
    if (is.null(idata) && is.null(imatrix)) {
      Y <- model.response(data) 
      #		SSPE <- manova$SSPE
    } 
    else {
      if (is.null(iterm)) stop("Must specify a within-S iterm for repeated measures designs" )
      ### FIXME::car -- workaround for car::Anova.mlm bug: no names assigned to $P component
      if (is.null(names(manova$P))) names(manova$P) <- names(manova$SSPE)
      Y <- model.response(data) %*% manova$P[[iterm]]
      #		SSPE <- manova$SSPE[[iterm]]
    }   
    
    vars <- colnames(Y)
    if (!missing(variables)){
      if (is.numeric(variables)) {
        vars <- vars[variables]
        if (any(is.na(vars))) stop("Bad response variable selection.")
      }
      else {
        check <- !(variables %in% vars)
        if (any(check)) stop(paste("The following", 
                                   if (sum(check) > 1) "variables are" else "variable is",
                                   "not in the model:", paste(variables[check], collapse=", ")))
        vars <- variables
      }
    }
    
    if(missing(var.labels)) var.labels <- vars
    else {
      if (length(var.labels) < length(vars)) stop("Too few var.labels supplied")
    }
    
    n.resp <- length(vars)
    if (n.resp < 3) stop("Fewer than 3 response variables.")
    range <- apply(Y, 2, range)
    min <- - offset.axes
    max <- 1 + offset.axes
    old.par <- par(mfrow=c(n.resp, n.resp), mar=rep(0,4))
    on.exit(par(old.par))
    
    panel.label <- function(x, ...) {
      plot(c(min, max),c(min, max), type="n", axes=FALSE)
      text(0.5, 0.5, var.labels[i], cex=var.cex)
      text(1, 0, signif(range[1, i], digits=digits), adj=c(1, 0))
      text(0, 1, signif(range[2, i], digits=digits), adj=c(0, 1)) 
      box()
    }	
    for (i in 1:n.resp){
      for (j in 1:n.resp){
        if (i == j){
          panel.label()
        }
        else {
          heplot(x, variables=c(vars[j], vars[i]), manova=manova, axes=FALSE,
                 idata=idata, idesign=idesign, imatrix=imatrix, iterm=iterm,
                 offset.axes=offset.axes, fill=fill, fill.alpha=fill.alpha, ...)
          box()
        }
      }
    }
  }

# Change log
# last modified 23 January 2007 by J. Fox
# last modified 14 May 2007 by M. Friendly -- return xlim, ylim
# last modified 18 May 2007 by M. Friendly -- fix xlim, ylim return when !add
# last modified 20 May 2007 by M. Friendly -- pass ... to text
# last modified 23 May 2007 by J. Fox -- add ... to call to points()
# last modified 22 Oct 2007 by M. Friendly
# -- moved lambda.crit to utility.R
# -- added he.rep to handle common task of repeating HE argument values
#  13 Apr 2009 by M. Friendly -- fix label.ellipse
#  15 Apr 2009 by M. Friendly -- added axes= to fix warnings from pairs.mlm
#  24 Dec 2009 by M. Friendly -- added idate=, idesign=, icontrasts, iterm for repeated measures
#  26 Dec 2009 by M. Friendly -- workaround for car::Anova buglet
#  27 Dec 2009 by M. Friendly -- made it work for designs with no between effects
#  28 Dec 2009 by M. Friendly -- made it work with car 2.0 for doubly multivariate
#  10 Jan 2010 by M. Friendly -- merged with heplot.mlm.R
#  23 Jul 2010 by M. Friendly -- return radius
#  05 Nov 2010 by M. Friendly 
# -- added fill= and fill.alpha for filled ellipses
# -- replaced lines() with polygon() for H and E ellipses
# -- calculate H.rank to distinguish degenerate ellipses
# -- added last() to utility.R
# -- added err.label to allow changing label for Error ellipse
# -- changed default colors from palette()[-1] to a better collection, also allowing options("heplot.colors")
#  15 Jan 2013 by M. Friendly
# -- replaced internal label.ellipse with separate function; added label.pos= argument
#  22 Feb 2013
# -- added ... to label.ellipse to be able to pass cex=



#' Two-Dimensional HE Plots
#' 
#' This function plots ellipses representing the hypothesis and error
#' sums-of-squares-and-products matrices for terms and linear hypotheses in a
#' multivariate linear model.  These include MANOVA models (all explanatory
#' variables are factors), multivariate regression (all quantitative
#' predictors), MANCOVA models, homogeneity of regression, as well as repeated
#' measures designs treated from a multivariate perspective.
#' 
#' The \code{heplot} function plots a representation of the covariance ellipses
#' for hypothesized model terms and linear hypotheses (H) and the corresponding
#' error (E) matrices for two response variables in a multivariate linear model
#' (mlm).
#' 
#' The plot helps to visualize the nature and dimensionality response variation
#' on the two variables jointly in relation to error variation that is
#' summarized in the various multivariate test statistics (Wilks' Lambda,
#' Pillai trace, Hotelling-Lawley trace, Roy maximum root). Roy's maximum root
#' test has a particularly simple visual interpretation, exploited in the
#' \code{size="evidence"} version of the plot. See the description of argument
#' \code{alpha}.
#' 
#' For a 1 df hypothesis term (a quantitative regressor, a single contrast or
#' parameter test), the H matrix has rank 1 (one non-zero latent root of \eqn{H
#' E^{-1}}) and the H "ellipse" collapses to a degenerate line.
#' 
#' Typically, you fit a mlm with \code{mymlm <- lm(cbind(y1, y2, y3, ...) ~
#' modelterms)}, and plot some or all of the \code{modelterms} with
#' \code{heplot(mymlm, ...)}.  Arbitrary linear hypotheses related to the terms
#' in the model (e.g., contrasts of an effect) can be included in the plot
#' using the \code{hypotheses} argument.  See
#' \code{\link[car]{linearHypothesis}} for details.
#' 
#' For repeated measure designs, where the response variables correspond to one
#' or more variates observed under a within-subject design, between-subject
#' effects and within-subject effects must be plotted separately, because the
#' error terms (E matrices) differ.  When you specify an intra-subject term
#' (\code{iterm}), the analysis and HE plots amount to analysis of the matrix
#' \bold{Y} of responses post-multiplied by a matrix \bold{M} determined by the
#' intra-subject design for that term.  See Friendly (2010) or the
#' \code{vignette("repeated")} in this package for an extended discussion and
#' examples.
#' 
#' The related \code{\link[candisc]{candisc}} package provides functions for
#' visualizing a multivariate linear model in a low-dimensional view via a
#' generalized canonical discriminant analyses.
#' \code{\link[candisc]{heplot.candisc}} and
#' \code{\link[candisc]{heplot3d.candisc}} provide a low-rank 2D (or 3D) view
#' of the effects for a given term in the space of maximum discrimination.
#' 
#' When an element of \code{fill} is \code{TRUE}, the ellipse outline is drawn
#' using the corresponding color in \code{col}, and the interior is filled with
#' a transparent version of this color specified in \code{fill.alpha}.  To
#' produce filled (non-degenerate) ellipses without the bounding outline, use a
#' value of \code{lty=0} in the corresponding position.
#' 
#' @aliases heplot heplot.mlm
#' @param mod a model object of class \code{"mlm"}.
#' @param terms a logical value or character vector of terms in the model for
#'              which to plot hypothesis matrices; if missing or \code{TRUE}, defaults to
#'              all terms; if \code{FALSE}, no terms are plotted.
#' @param hypotheses optional list of linear hypotheses for which to plot
#'              hypothesis matrices; hypotheses are specified as for the
#'              \code{\link[car]{linearHypothesis}} function in the \code{car} package; the
#'              list elements can be named, in which case the names are used.
#' @param term.labels logical value or character vector of names for the terms
#'              to be plotted. If \code{TRUE} (the default) the names of the terms are used;
#'              if \code{FALSE}, term labels are not plotted.
#' @param hyp.labels logical value or character vector of names for the
#'              hypotheses to be plotted. If \code{TRUE} (the default) the names of
#'              components of the list of hypotheses are used; if \code{FALSE}, hypothesis
#'              labels are not plotted.
#' @param err.label Label for the error ellipse
#' @param label.pos Label position, a vector of integers (in \code{0:4}) or
#'              character strings (in \code{c("center", "bottom", "left", "top", "right")},
#'              or in \code{c("C", "S", "W", "N", "E")} use in labeling ellipses, recycled
#'              as necessary.  Values of 1, 2, 3 and 4, respectively indicate positions
#'              below, to the left of, above and to the right of the max/min coordinates of
#'              the ellipse; the value 0 specifies the centroid of the \code{ellipse}
#'              object.  The default, \code{label.pos=NULL} uses the correlation of the
#' \code{ellipse} to determine "top" (r>=0) or "bottom" (r<0).  Even more
#'              flexible options are described in \code{\link{label.ellipse}}
#' @param variables indices or names of the two response variables to be
#'              plotted; defaults to \code{1:2}.
#' @param error.ellipse if \code{TRUE}, plot the error ellipse; defaults to
#'              \code{TRUE}, if the argument \code{add} is \code{FALSE} (see below).
#' @param factor.means logical value or character vector of names of factors
#'              for which the means are to be plotted, or \code{TRUE} or \code{FALSE};
#'              defaults to \code{TRUE}, if the argument \code{add} is \code{FALSE} (see
#'              below).
#' @param grand.mean if \code{TRUE}, plot the centroid for all of the data;
#'              defaults to \code{TRUE}, if the argument \code{add} is \code{FALSE} (see
#'              below).
#' @param remove.intercept if \code{TRUE} (the default), do not plot the
#'              ellipse for the intercept even if it is in the MANOVA table.
#' @param type ``type'' of sum-of-squares-and-products matrices to compute; one
#'              of \code{"II"}, \code{"III"}, \code{"2"}, or \code{"3"}, where \code{"II"}
#'              is the default (and \code{"2"} is a synonym).
#' @param idata an optional data frame giving a factor or factors defining the
#'              intra-subject model for multivariate repeated-measures data.  See Friendly
#'              (2010) and Details of \code{\link[car]{Anova}} for an explanation of the
#'              intra-subject design and for further explanation of the other arguments
#'              relating to intra-subject factors.
#' @param idesign a one-sided model formula using the ``data'' in idata and
#'              specifying the intra-subject design for repeated measure models.
#' @param icontrasts names of contrast-generating functions to be applied by
#'              default to factors and ordered factors, respectively, in the within-subject
#'              ``data''; the contrasts must produce an intra-subject model matrix in which
#'              different terms are orthogonal. The default is c("contr.sum", "contr.poly").
#' @param imatrix In lieu of \code{idata} and \code{idesign}, you can specify
#'              the intra-subject design matrix directly via \code{imatrix}, in the form of
#'              list of named elements.  Each element gives the columns of the
#'              within-subject model matrix for an intra-subject term to be tested, and must
#'              have as many rows as there are responses; the columns of the within-subject
#'              model matrix for \emph{different} terms must be mutually orthogonal.
#' @param iterm For repeated measures designs, you must specify one
#'              intra-subject term (a character string) to select the SSPE (E) matrix used
#'              in the HE plot.  Hypothesis terms plotted include the \code{iterm} effect as
#'              well as all interactions of \code{iterm} with \code{terms}.
#' @param markH0 A logical value (or else a list of arguments to
#'              \code{\link{mark.H0}}) used to draw cross-hairs and a point indicating the
#'              value of a point null hypothesis.  The default is TRUE if \code{iterm} is
#'              non-NULL.
#' @param manova optional \code{Anova.mlm} object for the model; if absent a
#'              MANOVA is computed. Specifying the argument can therefore save computation
#'              in repeated calls.
#' @param size how to scale the hypothesis ellipse relative to the error
#'              ellipse; if \code{"evidence"}, the default, the scaling is done so that a
#'              ``significant'' hypothesis ellipse at level \code{alpha} extends outside of
#'              the error ellipse; if \code{"effect.size"}, the hypothesis ellipse is on the
#'              same scale as the error ellipse.
#' @param level equivalent coverage of ellipse for normally-distributed errors,
#'              defaults to \code{0.68}, giving a standard 1 SD bivariate ellipse.
#' @param alpha significance level for Roy's greatest-root test statistic; if
#'              \code{size="evidence"}, then the hypothesis ellipse is scaled so that it
#'              just touches the error ellipse at the specified alpha level; a larger
#'              hypothesis ellipse \emph{somewhere} in the space of the response variables
#'              therefore indicates statistical significance; defaults to \code{0.05}.
#' @param segments number of line segments composing each ellipse; defaults to \code{60}.
#' @param center.pch character to use in plotting the centroid of the data;
#'              defaults to \code{"+"}.
#' @param center.cex size of character to use in plotting the centroid of the data; defaults to \code{2}.
#' @param col a color or vector of colors to use in plotting ellipses; the
#'              first color is used for the error ellipse; the remaining colors --- recycled
#'              as necessary --- are used for the hypothesis ellipses.  A single color can
#'              be given, in which case it is used for all ellipses.  For convenience, the
#'              default colors for all heplots produced in a given session can be changed by
#'              assigning a color vector via \code{options(heplot.colors =c(...)}.
#'              Otherwise, the default colors are \code{c("red", "blue", "black",
#'              "darkgreen", "darkcyan", "magenta", "brown", "darkgray")}.
#' @param lty vector of line types to use for plotting the ellipses; the first
#'              is used for the error ellipse, the rest --- possibly recycled --- for the
#'              hypothesis ellipses; a single line type can be given. Defaults to \code{2:1}.
#' @param lwd vector of line widths to use for plotting the ellipses; the first
#'              is used for the error ellipse, the rest --- possibly recycled --- for the
#'              hypothesis ellipses; a single line width can be given. Defaults to
#'              \code{1:2}.
#' @param fill A logical vector indicating whether each ellipse should be
#'              filled or not.  The first value is used for the error ellipse, the rest ---
#'              possibly recycled --- for the hypothesis ellipses; a single fill value can
#'              be given.  Defaults to FALSE for backward compatibility. See Details below.
#' @param fill.alpha Alpha transparency for filled ellipses, a numeric scalar
#'              or vector of values within \code{[0,1]}, where 0 means fully transparent and 1 means fully opaque.
#' @param xlab x-axis label; defaults to name of the x variable.
#' @param ylab y-axis label; defaults to name of the y variable.
#' @param main main plot label; defaults to \code{""}.
#' @param xlim x-axis limits; if absent, will be computed from the data.
#' @param ylim y-axis limits; if absent, will be computed from the data.
#' @param axes Whether to draw the x, y axes; defaults to \code{TRUE}
#' @param offset.axes proportion to extend the axes in each direction if
#'              computed from the data; optional.
#' @param add if \code{TRUE}, add to the current plot; the default is
#'              \code{FALSE}.  If \code{TRUE}, the error ellipse is not plotted.
#' @param verbose if \code{TRUE}, print the MANOVA table and details of
#'              hypothesis tests; the default is \code{FALSE}.
#' @param warn.rank if \code{TRUE}, do not suppress warnings about the rank of
#'              the hypothesis matrix when the ellipse collapses to a line; the default is \code{FALSE}.
#' @param \dots arguments to pass down to \code{plot}, \code{text}, and \code{points}.
#' @return The function invisibly returns an object of class \code{"heplot"},
#' with coordinates for the various hypothesis ellipses and the error ellipse,
#' and the limits of the horizontal and vertical axes.  These may be useful for
#' adding additional annotations to the plot, using standard plotting
#' functions.  (No methods for manipulating these objects are currently
#' available.)
#' 
#' The components are:
#' \describe{
#' \item{H}{a list containing the coordinates of each ellipse for the hypothesis terms} 
#' \item{E}{a matrix containing the coordinates for the error ellipse} 
#' \item{center}{x,y coordinates of the centroid} 
#' \item{xlim}{x-axis limits} 
#' \item{ylim}{y-axis limits}
#' \item{radius}{the radius for the unit circles used to generate the ellipses}
#' }
#' 
#' @seealso \code{\link[car]{Anova}}, \code{\link[car]{linearHypothesis}} for
#' details on testing MLMs.
#' 
#' \code{\link{heplot1d}}, \code{\link{heplot3d}}, \code{\link{pairs.mlm}},
#' \code{\link{mark.H0}} for other HE plot functions.
#' \code{\link{coefplot.mlm}} for plotting confidence ellipses for parameters
#' in MLMs.
#' 
#' \code{\link{trans.colors}} for calculation of transparent colors.
#' \code{\link{label.ellipse}} for labeling positions in plotting H and E
#' ellipses.
#' 
#' \code{\link[candisc]{candisc}}, \code{\link[candisc]{heplot.candisc}} for
#' reduced-rank views of \code{mlm}s in canonical space.
#' 
#' @references 
#' Friendly, M. (2006).  Data Ellipses, HE Plots and Reduced-Rank
#' Displays for Multivariate Linear Models: SAS Software and Examples
#' \emph{Journal of Statistical Software}, \bold{17}(6), 1--42. %
#' \url{https://www.jstatsoft.org/v17/i06/},
#' DOI: 10.18637/jss.v017.i06
#' 
#' Friendly, M. (2007).  HE plots for Multivariate General Linear Models.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{16}(2)
#' 421--444.  \url{http://datavis.ca/papers/jcgs-heplots.pdf}
#' 
#' Friendly, Michael (2010). HE Plots for Repeated Measures Designs.
#' \emph{Journal of Statistical Software}, 37(4), 1-40.  
#' DOI: 10.18637/jss.v037.i04.
#' 
#' Fox, J., Friendly, M. & Weisberg, S. (2013). Hypothesis Tests for
#' Multivariate Linear Models Using the car Package. \emph{The R Journal},
#' \bold{5}(1),
#' \url{https://journal.r-project.org/archive/2013-1/fox-friendly-weisberg.pdf}.
#' 
#' Friendly, M. & Sigal, M. (2014) Recent Advances in Visualizing Multivariate
#' Linear Models. \emph{Revista Colombiana de Estadistica}, \bold{37}, 261-283.
#' %\url{http://ref.scielo.org/6gq33g}.
#' @keywords hplot aplot multivariate
#' @examples
#' 
#' ## iris data
#' contrasts(iris$Species) <- matrix(c(0,-1,1, 2, -1, -1), 3,2)
#' contrasts(iris$Species)
#' 
#' iris.mod <- lm(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~
#' Species, data=iris)
#' 
#' hyp <- list("V:V"="Species1","S:VV"="Species2")
#' heplot(iris.mod, hypotheses=hyp)
#' # compare with effect-size scaling
#' heplot(iris.mod, hypotheses=hyp, size="effect", add=TRUE)
#' 
#' # try filled ellipses; include contrasts
#' heplot(iris.mod, hypotheses=hyp, fill=TRUE, 
#'        fill.alpha=0.2, col=c("red", "blue"))
#' heplot(iris.mod, hypotheses=hyp, fill=TRUE, 
#'        col=c("red", "blue"), lty=c(0,0,1,1))
#'
#' # vary label position and fill.alpha
#' heplot(iris.mod, hypotheses=hyp, fill=TRUE, fill.alpha=c(0.3,0.1), col=c("red", "blue"), 
#'        lty=c(0,0,1,1), label.pos=0:3)
#' 
#' # what is returned?
#' hep <-heplot(iris.mod, variables=c(1,3),  hypotheses=hyp)
#' str(hep)
#' 
#' # all pairs
#' pairs(iris.mod, hypotheses=hyp, hyp.labels=FALSE)
#' 
#' 
#' ## Pottery data, from car package
#' data(Pottery, package = "carData")
#' pottery.mod <- lm(cbind(Al, Fe, Mg, Ca, Na) ~ Site, data=Pottery)
#' heplot(pottery.mod)
#' heplot(pottery.mod, terms=FALSE, add=TRUE, col="blue", 
#'   hypotheses=list(c("SiteCaldicot = 0", "SiteIsleThorns=0")),
#'   hyp.labels="Sites Caldicot and Isle Thorns")
#' 
#' ## Rohwer data, multivariate multiple regression/ANCOVA
#' #-- ANCOVA, assuming equal slopes
#' rohwer.mod <- lm(cbind(SAT, PPVT, Raven) ~ SES + n + s + ns + na + ss, data=Rohwer)
#' car::Anova(rohwer.mod)
#' col <- c("red", "black", "blue", "cyan", "magenta", "brown", "gray")
#' heplot(rohwer.mod, col=col)
#' 
#' # Add ellipse to test all 5 regressors
#' heplot(rohwer.mod, hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss")), 
#'        col=col, fill=TRUE)
#' # View all pairs
#' pairs(rohwer.mod, hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss")))
#' # or 3D plot
#' 
#' if(requireNamespace("rgl")){
#' col <- c("pink", "black", "blue", "cyan", "magenta", "brown", "gray")
#' heplot3d(rohwer.mod, hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss")), col=col)
#' }
#' 
#' @export heplot
heplot <-
  function(mod, ...) UseMethod("heplot")

#' @rdname heplot
#' @exportS3Method  heplot mlm
#' @importFrom car linearHypothesis Anova
heplot.mlm <-
  function ( 
    mod,           # an mlm object
    terms,         # vector of terms to plot H ellipses
    hypotheses,    # list of linear hypotheses for which to plot H ellipses
    term.labels=TRUE,  # TRUE, FALSE or a vector of term labels of length(terms)
    hyp.labels=TRUE,   # as above for term.labels
    err.label="Error",
    label.pos=NULL,    # label positions: NULL or 0:4
    variables=1:2,     # x,y variables for the plot [variable names or numbers]
    error.ellipse=!add,
    factor.means=!add,
    grand.mean=!add,
    remove.intercept=TRUE,
    type=c("II", "III", "2", "3"),
    idata=NULL,
    idesign=NULL,
    icontrasts=c("contr.sum", "contr.poly"),
    imatrix=NULL,
    iterm=NULL,
    markH0=!is.null(iterm),
    manova,        # an optional Anova.mlm object
    size=c("evidence", "effect.size"),
    level=0.68,
    alpha=0.05,
    segments=60,   # line segments in each ellipse
    center.pch="+",   # doesn't have to be an argument
    center.cex=2,
    col=getOption("heplot.colors", c("red", "blue", "black", "darkgreen", "darkcyan","magenta", 
                                     "brown","darkgray")),
    # colors for H matrices, E matrix
    lty=2:1,
    lwd=1:2,
    fill=FALSE,         ## whether to draw filled ellipses (vectorized)
    fill.alpha=0.3,     ## alpha transparency for filled ellipses
    xlab,
    ylab,
    main="",
    xlim,           # min/max for X (override internal min/max calc) 
    ylim,
    axes=TRUE,      # whether to draw the axes
    offset.axes,    # if specified, the proportion by which to expand the axes on each end (e.g., .05)
    add=FALSE,      # add to existing plot?
    verbose=FALSE,
    warn.rank=FALSE,  
    ...) {
    ell <- function(center, shape, radius) {
      angles <- (0:segments)*2*pi/segments
      circle <- radius * cbind( cos(angles), sin(angles))
      if (!warn.rank){
        warn <- options(warn=-1)
        on.exit(options(warn))
      }
      Q <- chol(shape, pivot=TRUE)
      order <- order(attr(Q, "pivot"))
      t( c(center) + t( circle %*% Q[,order]))
    }
    #	label.ellipse <- function(ellipse, label, col){
    #		if (cor(ellipse)[1,2] >= 0){
    #			index <- which.max(ellipse[,2])
    #			x <- ellipse[index, 1] + 0.5 * strwidth(label)  # was: "A"
    #			y <- ellipse[index, 2] + 0.5 *strheight("A")
    #			adj <- c(1, 0) 
    #		}
    #		else {
    #			index <- which.min(ellipse[,2])
    #			x <- ellipse[index, 1] - 0.5 * strwidth(label)  # was: "A"
    #			y <- ellipse[index, 2] - 0.5 * strheight("A")
    #			adj <- c(0, 1) 
    #		}
    #		text(x, y, label, adj=adj, xpd=TRUE, col=col, ...)
    #	}
    #	last <- function(x) {x[length(x)]}
    
    #if (!require(car)) stop("car package is required.")
    # avoid deprecated warnings from car
    #	if (packageDescription("car")[["Version"]] >= 2) linear.hypothesis <- linearHypothesis
    type <- match.arg(type)
    size <- match.arg(size)
    data <- model.frame(mod)
    
    if (missing(manova)) {
      if (is.null(imatrix)) {
        manova <- car::Anova(mod, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts)
      }
      else {
        manova <- car::Anova(mod, type=type, idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix)
      } 
    }   
    if (verbose) print(manova)
    
    if (is.null(idata) && is.null(imatrix)) {
      Y <- model.response(data) 
      SSPE <- manova$SSPE
    } 
    else {
      if (is.null(iterm)) stop("Must specify a within-S iterm for repeated measures designs" )
      ### DONE::car -- workaround for car::Anova.mlm bug: no names assigned to $P component
      if (is.null(names(manova$P))) names(manova$P) <- names(manova$SSPE)
      Y <- model.response(data) %*% manova$P[[iterm]]
      SSPE <- manova$SSPE[[iterm]]
    }   
    
    if (!is.null(rownames(SSPE))) {response.names <- rownames(SSPE)}
    else {response.names <- paste("V.", 1:nrow(SSPE), sep="")}
    p <- length(response.names)
    
    if (!is.numeric(variables)) {
      vars <- variables
      variables <- match(vars, response.names)
      check <- is.na(variables)
      if (any(check)) stop(paste(vars[check], collapse=", "), 
                           " not among response variables.") 
    }
    else {
      if (any (variables > length(response.names))) stop("There are only ", 
                                                         length(response.names), " response variables.")
      vars <- response.names[variables]
    }
    if (length(variables) != 2) {
      extra <- if (length(variables) == 1) 'heplot1d()' else 
        if (length(variables) == 3) 'heplot3d()' else 'pairs()'
      stop(paste("You may only plot 2 response variables. Use", extra))
    }
    
    if (missing(terms) || (is.logical(terms) && terms)) {
      terms <- manova$terms
      # FIXME:  This does mot work if the between-S design includes only an intercept 
      # FIXME: && terms="(Intercept)" is specified
      if (!is.null(iterm)) {
        #			if (terms=="(Intercept)")  terms <- iterm else 
        terms <- terms[grep(iterm, terms)]   ## only include those involving iterm
      }
      if (remove.intercept) terms <- terms[terms != "(Intercept)"]
    }
    n.terms <- if (!is.logical(terms)) length(terms) else 0 
    # note: if logical here, necessarily FALSE
    n.hyp <- if (missing(hypotheses)) 0 else length(hypotheses)
    n.ell <- n.terms + n.hyp
    if (n.ell == 0) stop("Nothing to plot.")
    
    Y <- Y[,vars] 
    gmean <- if (missing(data))  c(0,0) 
    else colMeans(Y)
    if (missing(xlab)) xlab <- vars[1]
    if (missing(ylab)) ylab <- vars[2] 
    dfe <- manova$error.df
    scale <- 1/dfe
    radius <- sqrt(2 * qf(level, 2, dfe))
    
    # assign colors and line styles
    col <- he.rep(col, n.ell) 
    lty <- he.rep(lty, n.ell)
    lwd <- he.rep(lwd, n.ell)
    # handle filled ellipses
    fill <- he.rep(fill, n.ell)
    fill.alpha <- he.rep(fill.alpha, n.ell)
    fill.col <- trans.colors(col, fill.alpha)
    label.pos <- he.rep(label.pos, n.ell)
    # TODO:  take account of rank=1?
    fill.col <- ifelse(fill, fill.col, NA)
    E.col<- last(col)
    
    H.ellipse <- as.list(rep(0, n.ell))
    # keep track of ranks to distinguish degenerate ellipses
    H.rank <- rep(0, n.ell)
    
    if (n.terms > 0) for (term in 1:n.terms){
      term.name <- terms[term]
      H <- manova$SSP[[term.name]]
      if (!(all(variables %in% 1:nrow(H)))) {
        warning(paste("Skipping H term ", term.name, "(size: ", nrow(H), ")", sep=""))
        next
      }
      H <- H[variables, variables]
      dfh <- manova$df[term.name]
      factor <- if (size == "evidence") lambda.crit(alpha, p, dfh, dfe) else 1
      H <- H * scale/factor
      if (verbose){
        cat(term.name, " H matrix (", dfh, " df):\n")
        print(H)
      }
      H.ellipse[[term]] <- ell(gmean, H, radius)
      H.rank[term] <- qr(H)$rank
    }
    if (n.hyp > 0) for (hyp in 1:n.hyp){
      lh <- car::linearHypothesis(mod, hypotheses[[hyp]])
      H <- lh$SSPH[variables, variables]
      dfh <- lh$df
      factor <- if (size == "evidence") lambda.crit(alpha, p, dfh, dfe) else 1
      H <- H * scale/factor
      if (verbose){
        cat("\n\n Linear hypothesis: ", names(hypotheses)[[hyp]], "\n") 
        print(lh)
      }
      H.ellipse[[n.terms + hyp]] <- ell(gmean, H, radius)
    }
    E <- SSPE
    E <- E[variables, variables]
    E <- E * scale[1]
    E.ellipse <- ell(gmean, E, radius)
    H.ellipse$E <- E.ellipse     
    if (!add){
      max <- apply(sapply(H.ellipse, function(X) apply(X, 2, max)), 1, max)
      min <- apply(sapply(H.ellipse, function(X) apply(X, 2, min)), 1, min)
      factors <- data[, sapply(data, is.factor), drop=FALSE]
      if (!is.logical(factor.means)){
        factor.names <- colnames(factors) 
        which <- match(factor.means, factor.names)
        check <- is.na(which)
        if (any(check)) stop(paste(factor.means[check], collapse=", "), 
                             " not among factors.")
        factors <- factors[, which, drop=FALSE]
      }
      if (!is.logical(factor.means) || factor.means){   
        for (fac in factors){
          means <- aggregate(Y, list(fac), mean)
          min[1] <- min(min[1], means[,2])
          max[1] <- max(max[1], means[,2])
          min[2] <- min(min[2], means[,3])
          max[2] <- max(max[2], means[,3])
        }
      }
      if (!missing(offset.axes)){
        range <- max - min
        min <- min - offset.axes*range
        max <- max + offset.axes*range
      }
      xlim <- if(missing(xlim)) c(min[1], max[1]) else xlim
      ylim <- if(missing(ylim)) c(min[2], max[2]) else ylim
      plot(xlim, ylim,  type = "n", xlab=xlab, ylab=ylab, main=main, axes=axes, ...)
    }
    # no longer need H.ellipse$E, since we return it separately
    H.ellipse$E <- NULL
    if (grand.mean) 
      points(gmean[1], gmean[2], pch=center.pch, cex=center.cex, col="black", xpd=TRUE)
    if (error.ellipse){
      #		lines(E.ellipse, col=E.col, lty=lty[length(lty)], lwd=lwd[length(lwd)])
      polygon(E.ellipse, col=last(fill.col), border=last(col), lty=last(lty), lwd=last(lwd))
      label.ellipse(E.ellipse, err.label, col=last(col), label.pos=last(label.pos), ...)
    }
    term.labels <- if (n.terms == 0) NULL
    else if (!is.logical(term.labels)) term.labels
    else if (term.labels) terms else rep("", n.terms)  
    if (n.terms > 0) for (term in 1:n.terms){
      #			lines(H.ellipse[[term]], col=col[term], lty=lty[term], lwd=lwd[term])
      # TODO: avoid polygon if rank=1 ???
      polygon(H.ellipse[[term]], col=fill.col[term], border=col[term],  lty=lty[term], lwd=lwd[term])
      label.ellipse(H.ellipse[[term]], term.labels[term], col=col[term], label.pos=label.pos[term], ...) 
    }   
    hyp.labels <- if (n.hyp == 0) NULL
    else if (!is.logical(hyp.labels)) hyp.labels
    else if (hyp.labels) names(hypotheses) else rep("", n.hyp)  
    if (n.hyp > 0) for (hyp in 1:n.hyp){
      ell <- n.terms + hyp
      #			lines(H.ellipse[[ell]], col=col[ell], lty=lty[ell], lwd=lwd[ell])
      polygon(H.ellipse[[ell]], col=fill.col[ell], border=col[ell],  lty=lty[ell], lwd=lwd[ell])
      label.ellipse(H.ellipse[[ell]], hyp.labels[hyp], col=col[ell], label.pos=label.pos[ell], ...)
    }
    if (!add && (!is.logical(factor.means) || factor.means)){
      for (fac in factors){
        means <- aggregate(Y, list(fac), mean)
        points(means[,2], means[,3], pch=16, xpd=TRUE, ...)
        text(means[,2], means[,3], labels=as.character(means[,1]), pos=3, xpd=TRUE, ...)
      }
    }
    
    if(is.logical(markH0) && markH0) mark.H0()
    else if (is.list(markH0)) do.call(mark.H0, markH0)
    
    names(H.ellipse) <- c(if (n.terms > 0) term.labels, if (n.hyp > 0) hyp.labels)
    result <- if (!add) list(H=H.ellipse, E=E.ellipse, center=gmean, xlim=xlim, ylim=ylim, radius=radius)
    else list(H=H.ellipse, E=E.ellipse, center=gmean, radius=radius)
    class(result) <- "heplot"
    invisible(result)
  }


#' Internal heplots functions
#' 
#' Internal functions for the heplots package
#' 
#' These functions calculate critical values of multivariate test statistics (Wilks' Lambda, Hotelling-Lawley
#' trace, Roy's maximum root test) used in setting the size of H ellipses relative to E.
#' They are not intended to be called by the user.
#' 
#' 
#' @name heplots-internal
#' @aliases lambda.crit HLT.crit Roy.crit he.rep termInfo last
#' @param alpha significance level for critical values of multivariate
#' statistics
#' @param p Number of variables
#' @param dfh degrees of freedom for hypothesis
#' @param dfe degrees of freedom for error
#' @param test.statistic Test statistic used for the multivariate test
#' @param x An argument to \code{\link{heplot}} or \code{\link{heplot3d}} that
#' is to be repeated for Error and all hypothesis terms
#' @param n Number of hypothesis terms
#' @author Michael Friendly \email{friendly@yorku.ca}
#' @keywords internal
#' @return The critical value of the test statistic
#'
lambda.crit <- function(alpha, p, dfh, dfe, 
                        test.statistic=c("Roy", "HLT", "Hotelling-Lawley")){
  test.statistic <- match.arg(test.statistic)
  switch(test.statistic,
         Roy = Roy.crit(alpha, p, dfh, dfe),
         HLT = HLT.crit(alpha, p, dfh, dfe),
         "Hotelling-Lawley" = HLT.crit(alpha, p, dfh, dfe)
  )
}
# see: http://wiki.math.yorku.ca/index.php/Statistics:_Ellipses
## Critical value for \lambda_1 in Roy test

#' @rdname heplots-internal
Roy.crit <- function(alpha, p, dfh, dfe){
  df1 <- max(p, dfh)
  df2 <- dfe - df1 + dfh
  (df1/df2) * qf(alpha, df1, df2, lower.tail=FALSE)
}

## Critical value for \bar{\lambda_i} in HLT test
#' @rdname heplots-internal
HLT.crit <- function ( alpha, p, dfh, dfe) {
  s <- min(p, dfh)
  m <- (abs(p-dfh)-1)/2
  n <- (dfe-p-1)/2
  df1 <- 2*m + s + 1
  df2 <- 2*(s*n +1)
  s * (df1/df2) * qf(alpha, df1, df2, lower.tail=FALSE)	
}


# extend HE parameters for given number of terms
#   return vector in the form H1, H2, ..., E
#' @rdname heplots-internal
he.rep <- function (x, n) {
  if (length(x) < 2) x <- rep(x, 2)
  x <- c(rep(x[-1], n)[1:n], x[1])
  return(x)
}

last <- function(x) {x[length(x)]}

# copied from stats::: to avoid using :::
#' @rdname heplots-internal
Pillai <- function (eig, q, df.res) 
{
  test <- sum(eig/(1 + eig))
  p <- length(eig)
  s <- min(p, q)
  n <- 0.5 * (df.res - p - 1)
  m <- 0.5 * (abs(p - q) - 1)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * n + s + 1
  c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s * tmp2)
}

#' @rdname heplots-internal
Wilks <- function (eig, q, df.res) 
{
  test <- prod(1/(1 + eig))
  p <- length(eig)
  tmp1 <- df.res - 0.5 * (p - q + 1)
  tmp2 <- (p * q - 2)/4
  tmp3 <- p^2 + q^2 - 5
  tmp3 <- if (tmp3 > 0) 
    sqrt(((p * q)^2 - 4)/tmp3)
  else 1
  c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q, 
    p * q, tmp1 * tmp3 - 2 * tmp2)
}

#' @rdname heplots-internal
HL <- function (eig, q, df.res) 
{
  test <- sum(eig)
  p <- length(eig)
  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (df.res - p - 1)
  s <- min(p, q)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * (s * n + 1)
  c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
}

#' @rdname heplots-internal
Roy <- function (eig, q, df.res) 
{
  p <- length(eig)
  test <- max(eig)
  tmp1 <- max(p, q)
  tmp2 <- df.res - tmp1 + q
  c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
}

# make transparent colors, suitable for filling areas
# alpha transparency: (0 means fully transparent and 1 means opaque).
# names: optional character vector of names for the colors





#' Make Colors Transparent
#' 
#' Takes a vector of colors (as color names or rgb hex values) and adds a
#' specified alpha transparency to each.
#' 
#' Colors (\code{col}) and \code{alpha} need not be of the same length. The
#' shorter one is replicated to make them of the same length.
#' 
#' @param col A character vector of colors, either as color names or rgb hex
#' values
#' @param alpha alpha transparency value(s) to apply to each color (0 means
#' fully transparent and 1 means opaque)
#' @param names optional character vector of names for the colors
#' @return A vector of color values of the form \code{"#rrggbbaa"}
#' @author Michael Friendly
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' \code{\link[grDevices]{col2rgb}}, \code{\link[grDevices]{rgb}},
#' \code{\link[grDevices]{adjustcolor}},
#' @keywords color
#' @examples
#' 
#' trans.colors(palette(), alpha=0.5)
#' 
#' # alpha can be vectorized
#' trans.colors(palette(), alpha=seq(0, 1, length=length(palette())))
#' 
#' # lengths need not match: shorter one is repeated as necessary
#' trans.colors(palette(), alpha=c(.1, .2))
#' 
#' trans.colors(colors()[1:20])
#' 
#' # single color, with various alphas
#' trans.colors("red", alpha=seq(0,1, length=5))
#' # assign names
#' trans.colors("red", alpha=seq(0,1, length=5), names=paste("red", 1:5, sep=""))
#' 
#' 
#' @export trans.colors
trans.colors <- function(col, alpha=0.5, names=NULL) {
  nc <- length(col)
  na <- length(alpha)
  # make lengths conform, filling out to the longest
  if (nc != na) {
    col <- rep(col, length.out=max(nc,na))
    alpha <- rep(alpha, length.out=max(nc,na))
  }
  clr <-rbind(col2rgb(col)/255, alpha=alpha)
  col <- rgb(clr[1,], clr[2,], clr[3,], clr[4,], names=names)
  col
}

#'  Label an ellipse
#'
#' @description 
#'  \code{label.ellipse} is used to a draw text label on an ellipse at its center or
#'  somewhere around the periphery in a very flexible way.
#'
#' @details 
#' If \code{label.pos=NULL}, the function uses the sign of the correlation
#' represented by the ellipse to determine a position
#' at the top (\eqn{r>=0}) or bottom (\eqn{r<0}) of the ellipse.

#' Integer values of 0, 1, 2, 3 and 4, respectively indicate positions 
#' at the center, below, to the left of, above 
#' and to the right of the max/min coordinates of the ellipse.
#' Label positions can also be specified as the corresponding character strings
#' \code{c("center", "bottom", "left", "top", "right")}, or compass directions, 
#' \code{c("C", "S", "W", "N", "E")}, or  
#
#' Other integer \code{label.pos} values, \code{5:nrow(ellipse)} are taken as indices of the row coordinates
#' to be used for the ellipse label. 
#' Equivalently, \code{label.pos} can also be a \emph{fraction} in (0,1), interpreted
#' as the fraction of the way around the unit circle, counterclockwise from the point (1,0).
#'
#' @param ellipse A two-column matrix of coordinates for the ellipse boundary
#' @param label   Character string to be used as the ellipse label 
#' @param col     Label color
#' @param label.pos  Label position relative to the ellipse.  See details 
#' @param xpd     Should the label be allowed to extend beyond the plot limits?
#' @param tweak   A vector of two lengths used to tweak label positions
#' @param ...     Other parameters passed to \code{text}, e.g., \code{cex}, \dots
#' 
#' @author Michael Friendly
#' @export
#' @seealso \code{\link{heplot}}
#' @examples 
#' circle <- function(center=c(0,0), radius=1, segments=60) {
#'    angles <- (0:segments)*2*pi/segments
#'    circle <- radius * cbind( cos(angles), sin(angles))
#'    t( c(center) + t( circle ))
#' }
#' 
#' label_demo <- function(ell) {
#'   plot(-2:2, -2:2, type="n", asp=1, main="label.pos values and points (0:60)")
#'   lines(ell, col="gray")
#'   points(0, 0, pch="+", cex=2)
#'   
#'   labs <- c("center", "bot", "left", "top", "right")
#'   for (i in 0:4) {
#'     label.ellipse(ell, label=paste(i, ":", labs[i+1]), label.pos = i)
#'   }
#'   for( i in 5*c(1,2, 4,5, 7,8, 10,11)) {
#'     points(ell[i,1], ell[i,2], pch=16)
#'     label.ellipse(ell, label=i, label.pos=i)
#'   }
#' }
#' 
#' circ <- circle(radius=1.8)
#' label_demo(circ)
#' 
#' ell <-circ %*% chol(matrix( c(1, .5, .5, 1), 2, 2)) 
#' label_demo(ell)

#' @export label.ellipse
label.ellipse <- function(ellipse, label, col="black", 
                          label.pos=NULL, xpd=TRUE, 
                          tweak=0.5*c(strwidth("M"), strheight("M")), ...){
  
  ellipse <- as.matrix(ellipse)
  if (ncol(ellipse) < 2) stop("ellipse must be a 2-column matrix")
  
  if (is.null(label.pos)) {
    r = cor(ellipse, use="complete.obs")[1,2]
    label.pos <- if (r>0) 3 else 1
  }
  else if(length(label.pos) > 1) {
    warning("label.pos = ", paste(label.pos, collapse=", "), " has length ", length(label.pos), " only 1st used." )
    label.pos <- label.pos[1]    # only use 1st if > 1
  }
  
  #		index <- if (1:4 %% 2) ... 
  
  posn <- c("center", "bottom", "left", "top", "right")
  poss <- c("C",      "S",      "W",    "N",   "E")
  if (is.character(label.pos)) {
    if (label.pos %in% posn) label.pos <- pmatch(label.pos, posn, nomatch=3) - 1
    if (label.pos %in% poss) label.pos <- pmatch(label.pos, poss, nomatch=3) - 1
  }
  pos <- label.pos
  
  if (label.pos==1) {   # bottom
    index <- which.min(ellipse[,2])
    x <- ellipse[index, 1]
    y <- ellipse[index, 2] + tweak[2]
  }
  else if (label.pos==2) {   # left
    index <- which.min(ellipse[,1])
    x <- ellipse[index, 1] + tweak[1]
    y <- ellipse[index, 2]
  }
  else if (label.pos==3) {   # top
    index <- which.max(ellipse[,2])
    x <- ellipse[index, 1] 
    y <- ellipse[index, 2] - tweak[2]
  }
  else if (label.pos==4) {   # right
    index <- which.max(ellipse[,1])
    x <- ellipse[index, 1] - tweak[1]
    y <- ellipse[index, 2]
  }
  else if (label.pos==0) {   # center
    x <- mean(ellipse[, 1])
    y <- mean(ellipse[, 2]) - tweak[2]
    pos <-3
  }
  else  {  # use as index into ellipse coords
    if (0 < label.pos & label.pos < 1) 
      label.pos <- floor(label.pos * nrow(ellipse))
    index <- max(1, min(label.pos, nrow(ellipse)))
    x <- ellipse[index, 1]
    y <- ellipse[index, 2]
    pos <- 4 - floor(4*(index-1)/nrow(ellipse))
  }
  
  text(x, y, label, pos=pos, xpd=xpd, col=col, ...)
}


