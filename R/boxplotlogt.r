boxplotlogt <- function (df, s, p = 5000)
{
    abetas <- abs(rt (5000, df = df)) * exp( 0.5 * s)  
    boxplot (abetas, ylab = "")
    title (ylab = bquote(paste (abs(beta))) ) 
    title (xlab = sprintf("t (df=%g, log(scale)= %.1f)", df, 0.5*s))
}

#boxplotlogt (1, exp (-10))

