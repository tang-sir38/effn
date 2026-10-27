\name{effn}
\alias{effn}
\title{Effective sample size, minimum detectable effect and prescription}
\description{
Takes a single-cell count matrix with donor labels and returns, per gene, the
effective number of donors after the depth penalty, the minimum detectable
effect at the requested power, and how many additional donors a target effect
would require.
}
\usage{
effn(object, group, donor = NULL, assay = NULL, layer = "counts",
     n_boot = 200L, n_tech = 40L,
     tech = c("bootstrap", "splithalf"), seed = 20261021L,
     alpha = 0.05, power = 0.80, min_donor = 3L, paired = "auto",
     jack = TRUE, verbose = FALSE)
}
\arguments{
\item{object}{A genes x cells count matrix (\code{matrix}, \code{Matrix},
\code{data.frame}), a \code{SingleCellExperiment}, or a Seurat object (v4 or
v5).}

\item{group}{Vector of length \code{ncol(object)} giving the group of each
cell. Exactly two groups are allowed. For a \code{SingleCellExperiment} or a
\code{Seurat} object the name of a metadata column may be given instead, so
\code{effn(sce, "group")} and \code{effn(sce, sce$group)} are equivalent.}

\item{donor}{Optional vector of donor labels. When \code{NULL}, labels are
looked up in sample metadata; if none are found the function stops rather than
treating the matrix as a single replicate.}

\item{assay}{Assay name for \code{SingleCellExperiment} or Seurat input.}

\item{layer}{Layer (Seurat v5) or slot (Seurat v4) holding raw counts.}

\item{n_boot}{Number of donor-level bootstrap replicates for the confidence
interval of the median effective sample size and median MDE.}

\item{n_tech}{Number of within-donor resamples used to measure the technical
variance of the pseudobulk mean.}

\item{tech}{Estimator for the technical variance: \code{"bootstrap"} resamples
the donor's cells with replacement at the observed depth, \code{"splithalf"}
uses the cheaper split-half estimator whose \code{1/n} divisor only holds on a
linear scale.}

\item{seed}{Integer seed; identical input and seed give identical output.}

\item{alpha}{Two-sided type I error.}

\item{power}{Target power.}

\item{min_donor}{Minimum donors per group; below this the gene is reported as
unestimable rather than given a number.}

\item{paired}{\code{"auto"} (default) treats the design as paired when every
donor contributes to both groups; \code{TRUE} / \code{FALSE} force the choice.}

\item{jack}{Use a leave-one-donor-out jackknife for the reported intervals of
the median effective count and the median MDE (paired designs). The donor
bootstrap duplicates donors and shrinks the between-donor variance, so its
interval sits below the point estimate; the jackknife does not.}

\item{verbose}{Print stage timings to the message stream.}
}
\value{
An object of class \code{effn_report}, a list with elements
\code{per_gene} (one row per gene: \code{mu_cpm}, \code{mean_log2cpm},
\code{sigma_obs}, \code{sigma_tech}, \code{sigma_bio}, \code{rho_tech},
\code{n_eff}, \code{n_eff_group1}, \code{n_eff_group2}, \code{mde_log2},
\code{mde_fold}, \code{status}), \code{summary} (design-level quantities,
the two loss terms, \code{n_required_2fold} with its interquartile range,
and a one-sentence prescription), \code{boot} (bootstrap replicates) and
\code{params}.
}
\details{
Donor-level values are formed by summing raw counts within a donor first and
normalising once, giving \code{log2(CPM + 1)}. The observed pooled variance
across donors is the sum of a biological and a technical component. The
technical component is measured by resampling each donor's own cells with
replacement at the observed depth and taking the variance of the resulting
pseudobulk values, so no assumption about how that variance scales with depth
enters the estimate. The split-half estimator is retained as
\code{tech = "splithalf"} for comparison and for large matrices where the
cheaper estimator is preferred, but its \code{1/n} divisor holds only on a
linear scale. The effective number of donors is then
\code{n / (1 + rho_tech)} with \code{rho_tech = sigma_tech^2 / sigma_bio^2}.
The minimum detectable effect uses exact t quantiles with
\code{df = n1 + n2 - 2} for an unpaired design and \code{df = n - 1} on the
within-donor contrast for a paired one, where \code{n} counts donors rather
than effective donors; the normal approximation is not used. The effective
counts enter the standard error of the contrast, not its quantile, and using
them for the quantile as well would shrink the degrees of freedom by
\code{1 + rho} and inflate every reported MDE. The prescribed donor count is
the median over estimable genes of the number each gene would need, so that it
refers to the same typical gene as the reported median MDE.

Every gene carries a \code{status}. \code{all_zero} means the gene is absent
from every cell. \code{tech_dominated} means the observed donor-level variance
is no larger than the split-half technical variance, so no biological component
can be separated. \code{low_resolution} means fewer than two effective donors
remain after the depth penalty, where the t reference distribution has under one
degree of freedom and any interval would be an artefact of a vanishing
denominator. \code{low_expression} means a cell-level mean below 0.1 CPM.
\code{uninformative} covers genes with no usable degrees of freedom. Only
\code{ok} genes enter the medians reported in \code{summary}.
}
\examples{
set.seed(1)
G <- 60; d <- rep(paste0("d", 1:8), each = 20)
g <- rep(rep(c("ctrl", "case"), each = 4), each = 20)
mu <- exp(rnorm(G, 2, 0.6)); b <- matrix(rnorm(G * 8, 0, 0.5), G)
cnt <- sapply(seq_along(d), function(j)
  rnbinom(G, mu = mu * exp(b[, match(d[j], paste0("d", 1:8))]), size = 10))
rownames(cnt) <- paste0("g", 1:G)
r <- effn(cnt, g, d, n_boot = 20, n_tech = 20)
summary(r)
}
\seealso{
\code{\link{prescribe}}, \code{\link{plot_effn}}
}
