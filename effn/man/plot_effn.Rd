\name{plot_effn}
\alias{plot_effn}
\title{Plot an effn report}
\description{
Draws three panels: the distribution of the effective number of donors across
genes, the minimum detectable effect as a function of the number of donors per
group, and the loss ladder from cells to donors to effective donors.
}
\usage{
plot_effn(report, which = c("neff", "mde", "loss"), target_mde = NULL,
          col = c("#2B6CB0", "#C05621", "#2F855A"), ...)
}
\arguments{
\item{report}{An \code{effn_report} produced by \code{\link{effn}}.}

\item{which}{Character vector of panels to draw; any of \code{"neff"},
\code{"mde"} and \code{"loss"}.}

\item{target_mde}{Optional target effect on the log2 scale, drawn as a
horizontal line in the MDE panel.}

\item{col}{Three colours, one per panel.}

\item{...}{Further arguments passed to \code{\link[graphics]{par}}.}
}
\value{
Invisibly, a list of the data frames used for each panel.
}
\examples{
set.seed(1)
G <- 60; d <- rep(paste0("d", 1:8), each = 20)
g <- rep(rep(c("ctrl", "case"), each = 4), each = 20)
mu <- exp(rnorm(G, 2, 0.6)); b <- matrix(rnorm(G * 8, 0, 0.5), G)
cnt <- sapply(seq_along(d), function(j)
  rnbinom(G, mu = mu * exp(b[, match(d[j], paste0("d", 1:8))]), size = 10))
rownames(cnt) <- paste0("g", 1:G)
r <- effn(cnt, g, d, n_boot = 20, n_tech = 5)
plot_effn(r, target_mde = 1)
}
\seealso{
\code{\link{effn}}, \code{\link{prescribe}}
}
