# ------------------------------------------------------------
# Bounds for CDE and NDE under an IV-based potential outcomes model
# ------------------------------------------------------------
# IV (Z): 0 = Female, 1 = Male
# X: treatment / mediator (binary)
# Y: outcome (binary)
#
# The script:
# 1. stores observed joint counts p(Z, X, Y),
# 2. converts them to conditional distributions P(X, Y | Z),
# 3. enumerates latent response types for X and Y,
# 4. builds the linear constraints implied by the observed data,
# 5. solves LP problems for sharp bounds on CDE(0), CDE(1), NDE(0), NDE(1).
# ------------------------------------------------------------

library(lpSolve)

# ------------------------------------------------------------
# 1. Input data
# ------------------------------------------------------------

make_case_array <- function(case = c("case1", "case2")) {
  case <- match.arg(case)
  
  p <- array(
    0,
    dim = c(2, 2, 2),
    dimnames = list(
      Z = c("0", "1"),
      X = c("0", "1"),
      Y = c("0", "1")
    )
  )
  
  if (case == "case1") {
    p["0", "1", "1"] <- 19
    p["0", "1", "0"] <- 151
    p["0", "0", "1"] <- 7
    p["0", "0", "0"] <- 0
    p["1", "1", "1"] <- 182
    p["1", "1", "0"] <- 20
    p["1", "0", "1"] <- 205
    p["1", "0", "0"] <- 7
  }
  
  if (case == "case2") {
    p["0", "1", "1"] <- 527
    p["0", "1", "0"] <- 527
    p["0", "0", "1"] <- 55
    p["0", "0", "0"] <- 303
    p["1", "1", "1"] <- 880
    p["1", "1", "0"] <- 195
    p["1", "0", "1"] <- 5
    p["1", "0", "0"] <- 59
  }
  
  p
}

# Choose the case here
counts <- make_case_array("case1")

# Standardize to probabilities
p <- counts / sum(counts)

# Marginal P(Z)
pz <- c(sum(p["0", , ]), sum(p["1", , ]))
names(pz) <- c("0", "1")

# Conditional P(X, Y | Z = z)
p_cond_z0 <- p["0", , ] / pz["0"]
p_cond_z1 <- p["1", , ] / pz["1"]

# ------------------------------------------------------------
# 2. Latent response types
# ------------------------------------------------------------

# X potential types: (X(0), X(1))
X_TYPES <- rbind(
  c(0, 0), # never-taker
  c(1, 1), # always-taker
  c(0, 1), # complier
  c(1, 0)  # defier
)
rownames(X_TYPES) <- c("never", "always", "complier", "defier")

# Y potential types:
# columns correspond to:
# y00 = Y(X=0, Z=0)
# y01 = Y(X=0, Z=1)
# y10 = Y(X=1, Z=0)
# y11 = Y(X=1, Z=1)
Y_TYPES <- as.matrix(expand.grid(
  y00 = 0:1,
  y01 = 0:1,
  y10 = 0:1,
  y11 = 0:1
))

nX <- nrow(X_TYPES)   # 4
nY <- nrow(Y_TYPES)   # 16
m  <- nX * nY         # 64 latent classes

# Latent class indexing helpers
ix_of_k <- rep(seq_len(nX), each = nY)
iy_of_k <- rep(seq_len(nY), times = nX)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

# Map (x, z) to the appropriate Y-potential index in Y_TYPES
y_index <- function(x, z) {
  stopifnot(x %in% c(0, 1), z %in% c(0, 1))
  
  if (x == 0 && z == 0) return(1L)
  if (x == 0 && z == 1) return(2L)
  if (x == 1 && z == 0) return(3L)
  if (x == 1 && z == 1) return(4L)
  
  stop("Invalid (x, z) pair.")
}

# Realized X under latent X-type ix and instrument value z
X_of <- function(ix, z) {
  stopifnot(ix >= 1, ix <= nrow(X_TYPES), z %in% c(0, 1))
  X_TYPES[ix, z + 1]
}

# Realized Y under latent Y-type iy, treatment x, and instrument z
Y_of <- function(iy, x, z) {
  stopifnot(iy >= 1, iy <= nrow(Y_TYPES), x %in% c(0, 1), z %in% c(0, 1))
  Y_TYPES[iy, y_index(x, z)]
}

# ------------------------------------------------------------
# 4. Build equality constraints
# ------------------------------------------------------------

build_constraints <- function(p_cond_z0, p_cond_z1, ix_of_k, iy_of_k) {
  Aeq <- matrix(0, nrow = 9, ncol = m)
  beq <- numeric(9)
  
  row_id <- 0L
  
  for (z in 0:1) {
    p_cond <- if (z == 0) p_cond_z0 else p_cond_z1
    
    for (x in 0:1) {
      for (y in 0:1) {
        row_id <- row_id + 1L
        
        for (k in seq_len(m)) {
          ix <- ix_of_k[k]
          iy <- iy_of_k[k]
          
          Aeq[row_id, k] <- as.numeric(
            X_of(ix, z) == x && Y_of(iy, x, z) == y
          )
        }
        
        beq[row_id] <- p_cond[as.character(x), as.character(y)]
      }
    }
  }
  
  # Sum of latent probabilities = 1
  Aeq[9, ] <- 1
  beq[9] <- 1
  
  list(Aeq = Aeq, beq = beq)
}

constraints <- build_constraints(
  p_cond_z0 = p_cond_z0,
  p_cond_z1 = p_cond_z1,
  ix_of_k = ix_of_k,
  iy_of_k = iy_of_k
)

Aeq <- constraints$Aeq
beq <- constraints$beq

# ------------------------------------------------------------
# 5. Objective vectors
# ------------------------------------------------------------

# Controlled Direct Effect:
# CDE(x_fixed) = E[Y(x_fixed, Z=1) - Y(x_fixed, Z=0)]
make_cvec_CDE <- function(x_fixed, iy_of_k) {
  stopifnot(x_fixed %in% c(0, 1))
  
  cvec <- numeric(m)
  
  for (k in seq_len(m)) {
    iy <- iy_of_k[k]
    cvec[k] <- Y_of(iy, x_fixed, 1) - Y_of(iy, x_fixed, 0)
  }
  
  cvec
}

# Natural Direct Effect:
# NDE(0) = E[Y(X(0), Z=1) - Y(X(0), Z=0)]
make_cvec_NDE0 <- function(ix_of_k, iy_of_k) {
  cvec <- numeric(m)
  
  for (k in seq_len(m)) {
    ix <- ix_of_k[k]
    iy <- iy_of_k[k]
    
    x0 <- X_TYPES[ix, 1]  # X(0)
    y1 <- Y_of(iy, x0, 1)
    y0 <- Y_of(iy, x0, 0)
    
    cvec[k] <- y1 - y0
  }
  
  cvec
}

# NDE(1) = E[Y(X(1), Z=1) - Y(X(1), Z=0)]
make_cvec_NDE1 <- function(ix_of_k, iy_of_k) {
  cvec <- numeric(m)
  
  for (k in seq_len(m)) {
    ix <- ix_of_k[k]
    iy <- iy_of_k[k]
    
    x1 <- X_TYPES[ix, 2]  # X(1)
    y1 <- Y_of(iy, x1, 1)
    y0 <- Y_of(iy, x1, 0)
    
    cvec[k] <- y1 - y0
  }
  
  cvec
}

cvec_CDE0 <- make_cvec_CDE(0, iy_of_k)
cvec_CDE1 <- make_cvec_CDE(1, iy_of_k)
cvec_NDE0 <- make_cvec_NDE0(ix_of_k, iy_of_k)
cvec_NDE1 <- make_cvec_NDE1(ix_of_k, iy_of_k)

# ------------------------------------------------------------
# 6. LP solver
# ------------------------------------------------------------

solve_bounds <- function(cvec, Aeq, beq) {
  dir <- rep("=", nrow(Aeq))
  
  sol_min <- lp(
    direction = "min",
    objective.in = cvec,
    const.mat = Aeq,
    const.dir = dir,
    const.rhs = beq
  )
  
  sol_max <- lp(
    direction = "max",
    objective.in = cvec,
    const.mat = Aeq,
    const.dir = dir,
    const.rhs = beq
  )
  
  list(
    lower = sol_min$objval,
    upper = sol_max$objval,
    status_min = sol_min$status,
    status_max = sol_max$status,
    solution_min = sol_min$solution,
    solution_max = sol_max$solution
  )
}

# ------------------------------------------------------------
# 7. Estimate all bounds
# ------------------------------------------------------------

results <- list(
  CDE0 = solve_bounds(cvec_CDE0, Aeq, beq),
  CDE1 = solve_bounds(cvec_CDE1, Aeq, beq),
  NDE0 = solve_bounds(cvec_NDE0, Aeq, beq),
  NDE1 = solve_bounds(cvec_NDE1, Aeq, beq)
)

# ------------------------------------------------------------
# 8. Print a compact summary
# ------------------------------------------------------------

print_bounds <- function(results) {
  out <- data.frame(
    estimand = names(results),
    lower = sapply(results, `[[`, "lower"),
    upper = sapply(results, `[[`, "upper"),
    status_min = sapply(results, `[[`, "status_min"),
    status_max = sapply(results, `[[`, "status_max"),
    row.names = NULL
  )
  
  print(out, digits = 4)
  invisible(out)
}

print_bounds(results)