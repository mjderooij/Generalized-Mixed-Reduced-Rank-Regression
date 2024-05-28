read_data = function(){
library(haven) 
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(nnet)

ZA7953_v1_0_0 <- read_sav("ZA7953_v1-0-0.sav") # 37793 x 521
# View(ZA7953_v1_0_0)
# dim(ZA7953_v1_0_0) 


# FILTERING FOR NETHERLANDS 
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 26] == 1, ]
#dim(ZA7953_v1_0_0) 

# DATA FEATURING: creation of a new ordinal predictor 
education_columns <- ZA7953_v1_0_0[, 413:421]
ZA7953_v1_0_0$EDUCATION <- apply(education_columns, 1, function(row) which(row == 1)) 
ZA7953_v1_0_0$EDUCATION <- as.numeric(ZA7953_v1_0_0$EDUCATION)
#View(ZA7953_v1_0_0)

# DATA CLEANING FOR PREDICTOR VARIABLES

# POLITICAL ALIGNMENT
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 406] != 9, ]

# URBANIZATION
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 432] != 8, ]


# DATA CLEANING FOR ORDINAL RESPONSE VARIABLES 

ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 194] != 5, ]
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 211] != 5, ]
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 323] != 5, ]
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 328] != 5, ]
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 335] != 5, ]


# DATA CLEANING FOR BINARY RESPONSE VARIABLES 
# TRUST IN INSTITUTIONS 
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 150] != 3, ]
ZA7953_v1_0_0[, 150] = ZA7953_v1_0_0[, 150] - 1


# FUTURE ENLARGEMENT 
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 201] != 3, ]
ZA7953_v1_0_0 <- ZA7953_v1_0_0[ZA7953_v1_0_0[, 201] != 4, ]
ZA7953_v1_0_0[, 201] = ZA7953_v1_0_0[, 201] - 1



#dim(ZA7953_v1_0_0) # 847 x 522

############

X <- ZA7953_v1_0_0[, c(14, 406, 426, 432,522)]
Yo <- ZA7953_v1_0_0[, c(194, 211, 323, 328, 335)]
Yn = NULL
Yb = ZA7953_v1_0_0[,c(150, 201)]

# DELETING ALL THE NULL VALUES
complete_rows <- complete.cases(Yo, Yb, X)

Yo <- Yo[complete_rows, ]
Yo <- as.matrix(Yo)
colnames(Yo) <- c("Country Interests Respected in EU", "EU Minimum Wage", 
                  "EU Fiancial Support to Ukraine", "Defence Investment", "Renewable Energy Investment")


Yb <- Yb[complete_rows, ]
Yb <- as.matrix(Yb)
colnames(Yb) <- c("Trust in EU institutions", "Future Enlargement")


X <- X[complete_rows, ]
colnames(X) <- c("AGE", "POLITICAL ALIGNMENT", "GENDER", "URBANIZATION", "EDUCATION")
X = as.matrix(X)



# Factorization for predictor variables

X[ , 2] = as.factor(X[ , 2])
X[ , 3] = X[ , 3] - 1 
X[ , 3] = as.factor(X[ , 3])
X[ , 4] = as.factor(X[ , 4])
X[ , 5] = as.factor(X[ , 5])


# Xscale <- c("N", "O","C", "O", "O")
# 
# S = 2
# penalties = c(0,0,0)
# maxiter = 65536
# dcrit = 1e-6
# trace = TRUE
# 
# dim(Yo)
# dim(Yb)
# dim(Yn)
# dim(X)
# 
# 
# summary(X)
output = list(X = X, Yo = Yo, Yb = Yb, Yn = Yn)
return(output)
}