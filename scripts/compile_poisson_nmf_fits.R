# Compile the fitted Poisson non-negative factorizations into a single
# .RData file.
library(tools)
library(stringr)

# Combine results from all files of the form fit-*.rds in this
# directory.
out.dir <- "../output/droplet/rds"
  
# List all the RDS files containing the model fits.
files <- Sys.glob(file.path(out.dir,"fit-*.rds"))
n     <- length(files)

# Set up two data structures: "fits", a list used to store all the
# results; and "dat", a data frame summarizing the model parameters
# and optimization settings used to produce these fits.
fits   <- vector("list",n)
labels <- files
labels <- str_remove(labels,paste(out.dir,"/",sep = ""))
labels <- str_remove(labels,".rds")
names(fits) <- labels
dat <- data.frame(label       = labels,
                  k           = 0,
                  method      = "em",
                  extrapolate = FALSE,
                  stringsAsFactors = FALSE)

# Load the results from the RDS files.
for (i in 1:n) {
  out                  <- readRDS(files[i])
  fits[[i]]            <- out$fit
  dat[i,"k"]           <- ncol(out$fit$F)
  dat[i,"method"]      <- out$method
  dat[i,"extrapolate"] <- out$control$extrapolate
}

# Reorder the results in "fits" and "dat".
dat  <- transform(dat,method = factor(method,c("em","ccd","scd")))
i    <- with(dat,order(k,extrapolate,method))
dat  <- dat[i,]
fits <- fits[i]
rownames(dat) <- NULL

# Convert the "k" column to a factor.
dat <- transform(dat,k = factor(k))

# Save the combined results to an .RData file.
save(list = c("dat","fits"),
     file = "fits.RData")
resaveRdaFiles("fits.RData")
