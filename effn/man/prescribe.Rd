\name{prescribe}
\alias{prescribe}
\title{Prescribe the number of donors needed for a target effect}
\description{
Inverts the minimum detectable effect formula: given the biological standard
deviation and the depth penalty estimated from the data, returns how many donors
per group are required to detect \code{target_mde} at the requested power, and
how many are missing.
}
\usage{
prescribe(report, target_mde = 1, gene = NULL,
          alpha = report$params$alpha, power = report$params$power,
          nmax = 10000L)
}
\arguments{
\item{report}{An \code{effn_report} produced by \code{\link{effn}}.}

\item{target_mde}{Target effect on the log2 scale; \code{1} is a two-fold
change and \code{log2(1.5)} a 1.5-fold change.}

\item{gene}{Optional gene name. When \code{NULL} the median across estimable
genes is used.}

\item{alpha}{Two-sided type I error.}

\item{power}{Target power.}

\item{nmax}{Largest donor count the bisection will consider; when the target is
not reachable within \code{nmax} the result is \code{NA} and the sentence says
so.}
}
\value{
A list with \code{n_required}, \code{n_required_ci}, \code{n_current},
\code{shortfall}, \code{target_mde}, \code{target_fold} and \code{sentence}.
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
prescribe(r, target_mde = 1)
}
\seealso{
\code{\link{effn}}
}
